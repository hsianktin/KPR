# write the following template to /PROFILES/profile_name.jl
# ϕ₀ = [0,:E₊S]

# # setup simulation parameters
# k_pars = Pars( # correct substrate
#     k₁ = 0.1,
#     k₋₁ = 1,
#     τ = 3,
#     k₋₁⁺ = 1,
#     kₚ = 0.01,
#     T = 1e6
# )
# q_pars = Pars( # incorrect substrate
#     k₁ = 0.1,
#     k₋₁ = 2,
#     τ = 3,
#     k₋₁⁺ = 2,
#     kₚ = 0.01,
#     T = 1e6
# )
# times = [10.0^i for i in 0:0.2:6]

function write_profile_τ(profile_name, τ)
    string = """
    ϕ₀ = [0,:E₊S]

# setup simulation parameters
    k_pars = Pars( # correct substrate
        k₁ = 0.1,
        k₋₁ = 1,
        τ = $(τ),
        k₋₁⁺ = 1,
        kₚ = 1,
        T = 1e6
    )
    q_pars = Pars( # incorrect substrate
        k₁ = 0.1,
        k₋₁ = 2,
        τ = $(τ),
        k₋₁⁺ = 2,
        kₚ = 1,
        T = 1e6
    )
    times = [10.0^i for i in 0:0.2:6]
    """
    open("PROFILES/$profile_name.jl", "w") do f
        write(f, string)
    end
end

τs = [i for i in 0:0.2:5] ∪ [i for i in 0:0.1:1]
profile_names = ["profile_kp_1_τ_$τ" for τ in τs]
for (profile_name, τ) in zip(profile_names, τs)
    write_profile_τ(profile_name, τ)
end

#### verifying the N_f -> + ∞ limit
τs = [3]
profile_names = ["multistep_kp_1_τ_$(τ)_N_f_+∞" for  τ in τs]
for (profile_name, τ) in zip(profile_names, τs)
    write_profile_τ(profile_name, τ)
end
