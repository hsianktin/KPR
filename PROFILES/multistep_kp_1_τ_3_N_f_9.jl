ϕ₀ = [0,0]

# already large kₚ setting

# setup simulation parameters
k_pars = Pars( # correct substrate
    k₁ = 0.1,
    k₋₁ = 1,
    τ = 3,
    N_f = 9, # number of proofreading steps
    N_s = 9, # threshold for the stabilization step 
    k₋₁⁺ = 1,
    kₚ = 1,
    T = 5e3
)
q_pars = Pars( # incorrect substrate
    k₁ = 0.1,
    k₋₁ = 2,
    τ = 3,
    N_f = 9, # number of proofreading steps
    N_s = 9, # threshold for the stabilization step 
    k₋₁⁺ = 2,
    kₚ = 1,
    T = 5e3
)
times = [10.0^i for i in 0:0.2:6]    
