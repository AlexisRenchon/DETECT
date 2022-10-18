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

include("./DETECTModel_auxiliary.jl")

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
    E₀ₛ::FT 
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
    E₀ₛ::FT, 
    T₀::FT,
    α₄::FT,
    α₅::FT,
    BD::FT,
    ϕ₁₀₀::FT,
    PD::FT,
) where {FT}
    return DETECTParameters{FT}(Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀ₛ, T₀, α₄, α₅, BD, ϕ₁₀₀, PD)
end


# 2. Model 

"""
    DETECTModel

A model for simulating the production and transport of CO₂ in the soil with dynamic
source and diffusion terms.

$(DocStringExtensions.FIELDS)
"""
struct DETECTModel{FT, PS, D, BC, S, DT} <: AbstractModel{FT}
    "the parameter set"
    parameters::PS # in constructor you could enforce that it is ::DETECTParameters{FT}
    "the soil domain, using ClimaCore.Domains"
    domain::D
    "the boundary conditions, of type AbstractSoilBoundaryConditions"
    boundary_conditions::BC # maybe also have an FT
    "A tuple of sources, each of type AbstractSoilSource"
    sources::S
    " Drivers"
    driver::DT
end

ClimaLSM.name(model::DETECTModel) = :DETECT;


# 3. Prognostic and Auxiliary variables 

ClimaLSM.prognostic_vars(::DETECTModel) = (:C,) # pCO2 in soil, [ppm] # stored in Y.DETECT.C
ClimaLSM.prognostic_types(::DETECTModel{FT}) where {FT} = (FT,)
ClimaLSM.auxiliary_vars(::DETECTModel) = (:D, :Sₘ, :Sᵣ) # Diffusivity, Source (microbe + root) # stored in p.DETECT.D
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
        @unpack Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀ₛ, T₀, α₄, α₅, BD, ϕ₁₀₀, PD = model.parameters
        
	top_flux_bc, bot_flux_bc =
            boundary_fluxes(model.boundary_conditions, p, t)

        interpc2f = Operators.InterpolateC2F()
        gradc2f_C = Operators.GradientC2F()
        divf2c_C = Operators.DivergenceF2C(
            top = Operators.SetValue(Geometry.WVector.(top_flux_bc)),
            bottom = Operators.SetValue(Geometry.WVector.(bot_flux_bc)),
        ) # -∇ ⋅ (-D∇C), where -D∇C is a flux of C02. ∇C point in direction of increasing C, so the flux is - this.
        @. dY.DETECT.C =
	                 -divf2c_C(-interpc2f(p.DETECT.D)*gradc2f_C(Y.DETECT.C))

        # Source terms are added in here
        for src in model.sources
            source!(dY, src, Y, p, model.parameters)
        end

    end
    return rhs!
end

abstract type AbstractCarbonSource end
struct RootProduction end
struct MicrobeProduction end

function source!(dY, src::RootProduction, Y, p, params)

    dY .+= p.DETECT.Sᵣ
end

function source!(dY, src::MicrobeProduction, Y, p, params)
   
    dY .+= p.DETECT.Sₘ
end

# 5. Auxiliary variables
abstract type AbstractSoilDriver end
struct PrescribedSoil <: AbstractSoilDriver
    temperature::Function # (t,z) -> exp(-z)*sin(t), or e.g. a spline fit to data
    volumetric_liquid_fraction::Function
end


function soil_temperature(driver::PrescribedSoil, p, Y, t, z)
    return driver.temperature(t, z)
end
function soil_moisture(driver::PrescribedSoil, p, Y, t, z)
    return driver.volumetric_liquid_fraction(t, z)
end
#=
function soil_moisture(driver::PrognosticSoil, p, Y, t, z)
    return Y.soil.ϑ_l
end

function soil_temperature(driver::PrognosticSoil, p, Y, t, z)
    return p.soil.T
end
=#
    """
    make_update_aux(model::DETECTModel)

An extension of the function `make_update_aux`, for the DETECT equation. 
This function creates and returns a function which updates the auxiliary
variables `p.soil.variable` in place.

This has been written so as to work with Differential Equations.jl.
"""
function ClimaLSM.make_update_aux(model::DETECTModel)
    function update_aux!(p, Y, t) 
        @unpack Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, P, Dₗᵢ, E₀ₛ, T₀, α₄, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀ = model.parameters
        Ts = soil_temperature(model.drivers.soil.temperature,p, Y, t, z)
	@. p.DETECT.D = CO2_diffusivity(
					ϕ₁₀₀,
					diffusion_coefficient(Dstp, Ts, T₀, P₀, P),
					soil_porosity(BD, PD),
					) 

	@. p.DETECT.Sᵣ = root_source(
					       Rᵦ,
					       Rᶜ,
					       root_θ_adj(α₁ᵣ, α₂ᵣ, α₃ᵣ, α₄),
					       temp_adj(
							energy_act(E₀ₛ),
							T₀,
							),
    )
   @. p.DETECT.Sₘ = microbe_source(
						  Kₘ, 
						  CUE,
						  Vmax(
						       Vᵦ,
						       energy_act(E₀ₛ),
						       T₀
						       ),
						  Csol(Dₗᵢ, p),
						  Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ,
						  )

    end
    return update_aux!
end


struct FluxBC
    top::Function
    bottom::Function
end

function boundary_fluxes(f::FluxBC, p, t)
    return top(t), bottom(t)
end



