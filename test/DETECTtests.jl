FT = Float32
include("src/DETECTModel_auxiliary.jl")
include("src/DETECTModel_params.jl")

fTₛ(z,t) = exp(-z)*sin(t)*100
fθ(z,t) = exp(-z)*sin(t)
fθₐᵣ(z,t) = exp(-z)*sin(t) # weighted mean θ past 4 days
fθₐₘ(z,t) = exp(-z)*sin(t) # weighted mean θ past 4 days
fTₛₐ(z,t) = exp(-z)*sin(t)*100 # weighted mean Ts past 4 weeks

t = 1
z = 1
Tₛ = FT(fTₛ(z,t))
θ = FT(fθ(z,t))
θₐᵣ = FT(fθₐᵣ(z,t))
θₐₘ = FT(fθₐᵣ(z,t))
Tₛₐ = FT(fTₛₐ(z,t))

P = P₀

Dg₀ = diffusion_coefficient(Dstp, Tₛ, T₀, P₀, P)
ϕ = air_filled_soil_porosity(BD, PD, θ)
Dgs = CO₂_diffusivity(Dg₀, ϕ₁₀₀, ϕ, b) 
fᵣ = root_θ_adj(α₁ᵣ, α₂ᵣ, α₃ᵣ, θ, θₐᵣ) 
E₀ = energy_act(E₀ₛ, α₄, Tₛₐ) 
g = temp_adj(E₀, Tᵣₑ, T₀, Tₛ) 
Sᵣ = root_source(Rᵦ, Cᵣ, fᵣ, g) 
fₘ = microbe_θ_adj(α₁ₘ, α₂ₘ, α₃ₘ, θ, θₐₘ) 
Vmax = fVmax(Vᵦ, fₘ, g) 
Csol = fCsol(Csom, p, θ, Dₗᵢ) 
Sₘ = microbe_source(Vmax, Csol, Kₘ, Cmic, CUE) 
S = fS(Sᵣ, Sₘ) 

# TO DO: 
# list things that vary in depth f(z)
# list things that vary in depth and time f(z,t)

