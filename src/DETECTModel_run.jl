# code to run DETECTModel
include("./DETECTModel.jl")

FT = Float32
nelems = 50 # number of layers in the vertical
zmin = FT(-5)
zmax = FT(0.0)
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
top_bc = FT(765.0) # atm CO2 in mg m-3
bot_bc = FT(765.0) # let's just assume starting condition vertically homogeneous for now 
sources = (RootProduction(),MicrobeProduction()) # is S a source? Do I need to code differently?
boundary_fluxes = FluxBC(top_flux_bc, bot_flux_bc)
params = DETECTParameters{FT}(θ, Tₛ, θₐᵣ, θₐₘ, Tₛₐ, p, Cᵣ, Csom, Cmic, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀ₛ, T₀, α₄, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b)

θ_soil(t, z) = exp(-z)*sin(t)*100 + 273 # for now
soil_temp(t, z) = exp(-z)*sin(t) #

DETECT = DETECTModel{FT}(;
	parameters = params,
	domain = soil_domain,
	boundary_conditions = boundary_fluxes,
	sources = sources,
	driver = PrescribedSoil(;temperature = soil_temp, volumetric_liquid_fraction = θ_soil)
	)

Y, p, coords = initialize(DETECT)

# Initial conditions
function init_DETECT!(YDETECT, z, params)
	function CO2_profile(
		z::FT,
		params::DETECTParameters{FT},
	) where {FT}
		@unpack θ, Tₛ, θₐᵣ, θₐₘ, Tₛₐ, p, Cᵣ, Csom, Cmic, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀ₛ, T₀, α₄, α₅, BD, ϕ₁₀₀, PD, Dstp, P₀, b = params
		C = FT(765.0) # should be f(z)	
		return FT(C)
	end
	YDETECT.C .= CO2_profile.(z, Ref(params))
end

init_DETECT!(Y, coords.z, DETECT.parameters)

DETECT_ode! = make_ode_function(DETECT)

t0 = FT(0)
tf = FT(60) # why not
dt = FT(1)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values) #?
prob = ODEProblem(DETECT_ode!, Y, (t0, tf), p)
sol = solve(prob, Euler(); dt = dt, callback = cb) # do we want Euler or another algorithm for this?


