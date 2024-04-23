using ProgressMeter

#### CC_FLD with respect to τ
τs = [t for t in 0.0:0.2:5.0] ∪ [t for t in 0.0:0.1:1.0]
profile_names = ["profile_kp_1_τ_$τ" for τ in τs]

# when down, send email
try
    # for profile in profiles
    #     @show profile
    #     run_profile(profile)
    # end
    # permute the profiles
    @showprogress 1 for (profile_name, τ) in zip(profile_names, τs)
        @show cmd = `julia multiround_simulation_task_dispatcher.jl --profile PROFILES/$profile_name.jl`
        # run cmd and write stdout and stderr to file
        run(`julia multiround_simulation_task_dispatcher.jl --profile PROFILES/$profile_name.jl`)
    end
catch e
    
end
