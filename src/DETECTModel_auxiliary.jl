# define functions for DETECT auxiliary variables

# 1. CO2 Diffusivity

function diffusion_coefficient(Dstp, Ts, T₀, P₀, P)
	Dg₀ = Dstp * (Ts/T₀)^1.75 * (P₀/P)
end

function soil_porosity(BD, PD, θ)
	ϕₜ = 1 + BD/PD
	ϕ = ϕₜ - θ 
end

function CO2_diffusivity(Dg₀, ϕ₁₀₀, ϕ, b)
	Dgs = Dg₀ * (2ϕ₁₀₀^3 + 0.04ϕ₁₀₀) * (ϕ/ϕ₁₀₀)^(2 + 3/b)
end

# 2.1 CO2 source: root

function root_temp_θ_adj1(α₁ᵣ, α₂ᵣ, α₃ᵣ, θ, θₐᵣ)
	fᵣ = exp(α₁ᵣθ + α₂ᵣθₐᵣ + α₃ᵣθθₐᵣ)
end

function temp_θ_adj2(E₀, Tᵣₑ, T₀, Tₛ, T₀)
	g = exp(E₀ * (1/(Tᵣₑ - T₀) - 1/(Tₛ - T₀)))
end

function energy_act(E₀ₛ, α₄, Tₛₐ)
	E₀ = E₀ₛ + α₄Tₛₐ
end

function Root_source(Rᵦ, Cᵣ, fᵣ, gᵣ)
	Sᵣ = Rᵦ * Cᵣ * f * g
end

# 2.2 CO2 source: microbe

function microbe_temp_θ_adj1(α₁ₘ, α₂ₘ, α₃ₘ, θ, θₐₘ)
	fₘ = exp(α₁ₘθ + α₂ₘθₐₘ + α₃ₘθθₐₘ)
end

function Vmax(Vᵦ, fₘ, g)
	Vmax = Vᵦ * fₘ * g # g defined in 2.1
end

function Csol(Csom, p, θ, Dₗᵢ)
	Csol = Csom * p * θ^3 * Dₗᵢ
end

function Microbe_source(Vmax, Csol, Km, Cmic, CUE)
	Sₘ = Vmax * Csol/(Kₘ + Csol) * Cmic * (1 - CUE)
end

function S(Sᵣ, Sₘ)
	S = Sᵣ + Sₘ
end



