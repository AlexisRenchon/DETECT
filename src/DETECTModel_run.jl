# code to run DETECTModel

nelems = 50
soil_domain = Column(; zlim = (zmin, zmax), nelements = nelems)
top_flux_bc = FT(0.0) # should be catm
bot_flux_bc = FT(0.0)
sources = () # is S a source? Do I need to code differently?
boundary_fluxes = FluxBC{FT}(top_flux_bc, bot_flux_bc)
params = DETECTParameters{FT}(Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀, T₀, α₄, α₅, BD, ϕ₁₀₀, PD)

DETECT = DETECTModel{FT}(;
	parameters = params,
	domain = soil_domain,
	boundary_conditions = boundary_fluxes,
	sources = sources,
	)

Y, p, coords = initialize(DETECT)

# Initial conditions
function init_DETECT!(YDETECT, z, params)
	function CO2_profile(
		z::FT,
		params::DETECTParameters{FT},
	) where {FT}
		@unpack Rᶜ, Rᵦ, α₁ᵣ, α₂ᵣ, α₃ᵣ, Sᶜ, Mᶜ, Vᵦ, α₁ₘ, α₂ₘ, α₃ₘ, Kₘ, CUE, p, Dₗᵢ, E₀, T₀, α₄, α₅, BD, ϕ₁₀₀, PD = params
		z_∇ = FT(-10) # zmin
		S = #?
		C = #?
		return FT(C)
	end
	YDETECT.C .= CO2_profile.(z, Ref(params))
end

init_DETECT!(Y, coords.z, soil.parameters)

DETECT_ode! = make_ode_function(DETECT)

t0 = FT(0)
tf = FT(60) # why not
dt = FT(1)
cb = SavingCallback((u, t, integrator) -> copy(integrator.p), saved_values) #?
prob = ODEProblem(DETECT_ode!, Y, (t0, tf), p)
sol = solve(prob, Euler(); dt = dt, callback = cb) # do we want Euler or another algorithm for this?


