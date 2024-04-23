    ϕ₀ = [0,:E₊S]

# setup simulation parameters
    k_pars = Pars( # correct substrate
        k₁ = 0.1,
        k₋₁ = 1,
        τ = 4.6,
        k₋₁⁺ = 1,
        kₚ = 1,
        T = 1e6
    )
    q_pars = Pars( # incorrect substrate
        k₁ = 0.1,
        k₋₁ = 2,
        τ = 4.6,
        k₋₁⁺ = 2,
        kₚ = 1,
        T = 1e6
    )
    times = [10.0^i for i in 0:0.2:6]
    