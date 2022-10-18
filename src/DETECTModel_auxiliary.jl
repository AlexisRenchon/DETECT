# define functions for DETECT auxiliary variables

# 1. CO2 Diffusivity

function diffusion_coefficient(Dstp::FT, Ts::FT, T₀::FT, P₀::FT, P::FT) where {FT}
	Dg₀ = Dstp * (Ts/T₀)^FT(1.75) * (P₀/P)
	return Dg₀
end

function soil_porosity(BD::FT, PD::FT, θ::FT) where {FT}
	ϕₜ = 1 + BD/PD
	ϕ = ϕₜ - θ
	return ϕ
end

function CO2_diffusivity(Dg₀::FT, ϕ₁₀₀::FT, ϕ::FT, b::FT) where {FT}
	Dgs = Dg₀ * (FT(2)ϕ₁₀₀^FT(3) + FT(0.04)ϕ₁₀₀) * (ϕ/ϕ₁₀₀)^(FT(2) + FT(3)/b)
	return Dgs
end

# 2.1 CO2 source: root

function root_temp_θ_adj1(α₁ᵣ::FT, α₂ᵣ::FT, α₃ᵣ::FT, θ::FT, θₐᵣ::FT) where {FT}
	fᵣ = exp(α₁ᵣθ + α₂ᵣθₐᵣ + α₃ᵣθθₐᵣ)
	return fᵣ
end

function temp_θ_adj2(E₀::FT, Tᵣₑ::FT, T₀::FT, Tₛ::FT, T₀::FT) where {FT}
	g = exp(E₀ * (FT(1)/(Tᵣₑ - T₀) - FT(1)/(Tₛ - T₀)))
	return g
end

function energy_act(E₀ₛ::FT, α₄::FT, Tₛₐ::FT) where {FT}
	E₀ = E₀ₛ + α₄Tₛₐ
	return E₀
end

function Root_source(Rᵦ::FT, Cᵣ::FT, fᵣ::FT, gᵣ::FT) where {FT}
	Sᵣ = Rᵦ * Cᵣ * f * g
	return Sᵣ
end

# 2.2 CO2 source: microbe

function microbe_temp_θ_adj1(α₁ₘ::FT, α₂ₘ::FT, α₃ₘ::FT, θ::FT, θₐₘ::FT) where {FT}
	fₘ = exp(α₁ₘθ + α₂ₘθₐₘ + α₃ₘθθₐₘ) where {FT}
	return fₘ
end

function Vmax(Vᵦ::FT, fₘ::FT, g::FT) where {FT}
	Vmax = Vᵦ * fₘ * g # g defined in 2.1
	return Vmax
end

function Csol(Csom::FT, p::FT, θ::FT, Dₗᵢ::FT) where {FT}
	Csol = Csom * p * θ^FT(3) * Dₗᵢ
	return Csol
end

function Microbe_source(Vmax::FT, Csol::FT, Km::FT, Cmic::FT, CUE::FT) where {FT}
	Sₘ = Vmax * Csol/(Kₘ + Csol) * Cmic * (FT(1) - CUE)
	return Sₘ
end

function S(Sᵣ::FT, Sₘ::FT) where {FT}
	S = Sᵣ + Sₘ
	return S
end

