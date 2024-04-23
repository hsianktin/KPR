function write_profile_N_f(profile_name, N_f)
    string = """
    ϕ₀ = [0,0]

    # already large kₚ setting
    
    # setup simulation parameters
    k_pars = Pars( # correct substrate
        k₁ = 0.1,
        k₋₁ = 1,
        τ = 3,
        N_f = $N_f, # number of proofreading steps
        N_s = $N_f, # threshold for the stabilization step 
        k₋₁⁺ = 1,
        kₚ = 1,
        T = 5e3
    )
    q_pars = Pars( # incorrect substrate
        k₁ = 0.1,
        k₋₁ = 2,
        τ = 3,
        N_f = $N_f, # number of proofreading steps
        N_s = $N_f, # threshold for the stabilization step 
        k₋₁⁺ = 2,
        kₚ = 1,
        T = 5e3
    )
    times = [10.0^i for i in 0:0.2:6]    
    """
    open("PROFILES/multistep_$profile_name.jl", "w") do f
        write(f, string)
    end
end

function write_profile_N_f_τ(profile_name, N_f, τ)
    string = """
    ϕ₀ = [0,0]

    # already large kₚ setting
    
    # setup simulation parameters
    k_pars = Pars( # correct substrate
        k₁ = 0.1,
        k₋₁ = 1,
        τ = $τ,
        N_f = $N_f, # number of proofreading steps
        N_s = $N_f, # threshold for the stabilization step 
        k₋₁⁺ = 1,
        kₚ = 1,
        T = 5e3
    )
    q_pars = Pars( # incorrect substrate
        k₁ = 0.1,
        k₋₁ = 2,
        τ = $τ,
        N_f = $N_f, # number of proofreading steps
        N_s = $N_f, # threshold for the stabilization step 
        k₋₁⁺ = 2,
        kₚ = 1,
        T = 5e3
    )
    times = [10.0^i for i in 0:0.2:6]    
    """
    open("PROFILES/multistep_$profile_name.jl", "w") do f
        write(f, string)
    end
end

#### verifying the N_f -> + ∞ limit
N_fs = [i for i in [1,2,3,4,5,6,7,8,9,10,15,20,30,40,50,60,70,80,90,100,200]]
τs = [i for i in [3]]
profile_names = ["kp_1_τ_$(τ)_N_f_$(N_f)" for N_f in N_fs, τ in τs]
paras = [(N_f, τ) for N_f in N_fs, τ in τs]
for (profile_name, (N_f, τ)) in zip(profile_names, paras)
    write_profile_N_f_τ(profile_name, N_f, τ)
end

#### Fix N_f = 6, change τ
N_f = 6
τs = [i for i in 0.2:0.2:5.0] ∪ [i for i in 0.0:0.1:1.0]
profile_names = ["kp_1_τ_$(τ)_N_f_$(N_f)" for τ in τs]
for (profile_name, τ) in zip(profile_names, τs)
    write_profile_N_f_τ(profile_name, N_f, τ)
end