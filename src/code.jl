# This script aim to run DETECT as a standalone (not inside ClimaLSM.jl)
# I use ClimaLSM.jl tutorial and Katherine ppt 

# Required deps, from tutorial
using OrdinaryDiffEq: ODEProblem, solve, RK4
using SciMLBase
using Plots
using ClimaCore
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
using ClimaLSM
using ClimaLSM.Domains
import ClimaLSM: name, make_rhs, prognostic_vars, prognostic_types
import ClimaLSM.Domains: coordinates


# 1. Parameters  

"""
    DETECTParameters{FT <: AbstractFloat}
A struct for storing parameters of the `DETECTModel`.

$(DocStringExtensions.FIELDS)
"""
struct DETECTParameters{FT <: AbstractFloat}
  # root submodel parameters
    "Total root biomass C in a 1 m deep by 1 cm² soil column (mg C cm⁻²)"
    Rᶜ::FT
    "Root mass-base respiration rate at 10°C and mean environmental conditions (mg C cm⁻³ h⁻¹)"
    Rᵦ::FT
    "The effect of soil water content (θ) on root respiration (unitless)"
    α₁ᵣ::FT
    "The effect of antecedent θ on root respiration (unitless)"
    α₂ᵣ::FT
    "The interactive effect of θ and antecedent θ on root respiration (unitless)"
    α₃ᵣ::FT

  # microbial submodel parameters
    "Total soil organic C in a 1 m deep by 1 cm² soil column (mg C cm⁻²)"
    Sᶜ::FT
    "Total microbial biomass C in a 1 m deep by 1 cm² column of soil (mg C cm⁻²)"
    Mᶜ::FT
    "Value of Vₘₐₓ at 10°C and mean environmental conditions (mg C cm⁻³ h⁻¹)"
    Vᵦ::FT
    "The effect of soil water content (θ) on microbial respiration (unitless)"
    α₁ₘ::FT
    "The effect of antecedent θ on microbial respiration (unitless)"
    α₂ₘ::FT
    "The interactive effect of θ and antecedent θ on microbial respiration (unitless)"
    α₃ₘ::FT
    "Michaelis-Menten half-saturation constant (mg C cm⁻³ h⁻¹)"
    Kₘ::FT
    "Microbial carbon-use efficiency (mg C mg⁻¹ C⁻¹)"
    CUE::FT
    "Fraction of soil organic C that is soluble (-)"
    p::FT
    "Diffusivity of soil C substrate in liquid (unitless)"
    Dₗᵢ::FT

  # shared parameters between root/microbial submodels
    "Temperature sensitivity parameter, somewhat analogous to an energy activation (Kelvin)"
    E₀::FT 
    "Temperature sensitivity-related parameter (Kelvin)"
    T₀::FT
    "The effect of antecedent soil temperature on root and microbial respiration (unitless)"
    α₄::FT

  # soil CO₂ diffusivity submodel parameters
    "Absolute value of the slope of the line relating log(ψ) versus log(θ) (unitless)"
    α₅::FT
    "Soil bulk density (g cm⁻³)"
    BD::FT
    "Air-filled porosity at soil water potential of -100 cm H₂O (~ 10 kPa) (%)"
    ϕ₁₀₀::FT
    "Particle density"
    PD::FT
end

function DETECTParameters(;
    Rᶜ::FT,
    Rᵦ::FT,
    α₁ᵣ::FT,
    α₂ᵣ::FT,
    α₃ᵣ::FT,
    Sᶜ::FT,
    Mᶜ::FT,
    Vᵦ::FT,
    α₁ₘ::FT,
    α₂ₘ::FT,
    α₃ₘ::FT,
    Kₘ::FT,
    CUE::FT,
    p::FT,
    Dₗᵢ::FT,
    E₀::FT, 
    T₀::FT,
    α₄::FT,
    α₅::FT,
    BD::FT,
    ϕ₁₀₀::FT,
    PD::FT,
) where {FT}
    return DETECTParameters{FT}(Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀, T₀, α₄, α₅, BD, ϕ₁₀₀, PD)
end


# 2. Model 

"""
    DETECTModel

A model for simulating the production and transport of CO₂ in the soil with dynamic
source and diffusion terms.

$(DocStringExtensions.FIELDS)
"""
struct DETECTModel{FT, PS, D, BC, S} <: AbstractSoilModel{FT}
    "the parameter set"
    parameters::PS
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BC
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
end

ClimaLSM.name(model::DETECTModel) = :DETECT;


# 3. Prognostic and Auxiliary variables 

ClimaLSM.prognostic_vars(::DETECTModel) = (:C,) # pCO2 in soil, [ppm]
ClimaLSM.prognostic_types(::DETECTModel{FT}) where {FT} = (FT,)
ClimaLSM.auxiliary_vars(::DETECTModel) = (:D, :S) # Diffusivity, Source (microbe + root)
ClimaLSM.auxiliary_types(::DETECTModel{FT}) where {FT} = (FT, FT)


# 4. RHS

"""
    make_rhs(model::DETECTModel)

An extension of the function `make_rhs`, for the DETECT equation.
This function creates and returns a function which computes the entire
right hand side of the PDE for `C`, and updates `dY.soil.C` in place
with that value.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_rhs(model::DETECTModel)
    function rhs!(dY, Y, p, t)
        @unpack Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀, T₀, α₄, α₅, BD, ϕ₁₀₀, PD = model.parameters
        
	top_flux_bc, bot_flux_bc =
            boundary_fluxes(model.boundary_conditions, p, t)

        interpc2f = Operators.InterpolateC2F()

        gradc2f_C = Operators.GradientC2F()

        # We are setting a boundary value on a flux, which is a gradient of a scalar
        # Therefore, we should set boundary conditions in terms of a covariant vector
        # We set the third component first - supply a Covariant3Vector

        # Without topography only
        # In Cartesian coordinates, W (z^) = Cov3 (z^)= Contra3 (n^ = z^)
        # In spherical coordinates, W (r^) = Cov3 (r^) = Contra3 (n^ = r^)

        # It appears that the WVector is converted internally to a Covariant3Vector for the gradient value
        # at the boundary. Offline tests indicate that you get the same thing if
        # the bc is WVector(F) or Covariant3Vector(F*Δr) or Contravariant3Vector(F/Δr)

        divf2c_C = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(bot_flux_bc)),
        )

        # GradC2F returns a Covariant3Vector, so no need to convert.
        @. dY.DETECT.C =
	                 divf2c_C(interpc2f(Y.DETECT.D)*gradc2f_C(Y.DETECT.C)) + Y.DETECT.S # !!!!! not sure this is correct

        # Horizontal contributions
        horizontal_components!(dY, model.domain, model, p, z)

        # Source terms
        for src in model.sources
            source!(dY, src, Y, p)
        end

        # This has to come last
        dss!(dY, model.domain)
    end
    return rhs!
end


# 5. Auxiliary variables

"""
    make_update_aux(model::DETECTModel)

An extension of the function `make_update_aux`, for the DETECT equation. 
This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::DETECTModel)
    function update_aux!(p, Y, t)
        @unpack Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀, T₀, α₄, α₅, BD, ϕ₁₀₀, PD = model.parameters
	@. p.DETECT.D = CO2_diffusivity(
					ϕ₁₀₀,
					diffusion_coefficient(),
					soil_porosity(BD, PD),
					) 

	# FUNCTION FOR EFFECTIVE DIFFUSIVITY OF CO2 in Soil/soil_respiration_parameterizations.jl
			# hydraulic_conductivity(
            # K_sat,
            # vg_m,
            # effective_saturation(ν, Y.soil.ϑ_l, θ_r),
        # )
	@. p.DETECT.S = CO2_source(
				   Root_source(
					       Rᵦ,
					       Rᶜ,
					       root_temp_θ_adj1(α₁ᵣ, α₂ᵣ, α₃ᵣ, α₄),
					       root_temp_θ_adj2(
								energy_act(E₀),
								T₀,
								),
					       ),
				   Microbe_source(
						  Kₘ, 
						  CUE,
						  Vmax(Vᵦ, E₀, T₀),
						  Csol(Dₗᵢ, p),
						  Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ,
						  ),
				   )

	# FUNCTION FOR SOURCE (i.e., production) OF CO2 in Soil/soil_respiration_parameterizations.jl
	# pressure_head(vg_α, vg_n, vg_m, θ_r, Y.soil.ϑ_l, ν, S_s)
    end
    return update_aux!
end



















