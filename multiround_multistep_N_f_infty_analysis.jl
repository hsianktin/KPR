using BSON
using StatsBase, DataFrames,CSV

using DataFrames, CSV
using PGFPlotsX
using LaTeXStrings
using Optim


using ProgressMeter
group_label = "kp_1_τ_3"
profiles_∞ = ["PROFILES/multistep_$(group_label)_N_f_+∞.jl"]
tasks = [(substrate, id) for substrate in ["k", "q"] for id in 1:100]
labels_∞ = [split(profile, "/")[end] |> x -> split(x, ".jl")[1] for profile in profiles_∞]
Results_k_∞ = []
Results_q_∞ = []
for label in labels_∞
    results_k = []
    results_q = []
    for task in tasks
        substrate, id = task
        result = BSON.load(
            "DATA/multiround_simulation_results_$(substrate)_$(id)_$(label).bson"
            )[:results]
        if substrate == "k"
            push!(results_k, result)
        elseif substrate == "q"
            push!(results_q, result)
        end
    end
    results_k = vcat(results_k...)
    results_q = vcat(results_q...)
    push!(Results_k_∞, results_k)
    push!(Results_q_∞, results_q)
end

include("multiround_simulation_product_analysis_utils.jl")
function handle_results_∞(results_k,results_q,profile)
    include("multiround_simulation_utils.jl");
    include(profile);
    
    𝐓ₚs_k =  [result[4] for result in results_k]
    𝐓ₚs_q =  [result[4] for result in results_q]
    𝐓ₚs = [𝐓ₚs_k, 𝐓ₚs_q]
    function P_τₚᵀ(𝐓ₚs; T=1e6, ξ = 0, nₜₕ = 1) # the probability of τₚᵀ < T
        if ξ == 0
            # incorrect substrate
            𝐓ₚ = 𝐓ₚs[2]
        else
            # correct substrate
            𝐓ₚ = 𝐓ₚs[1]
        end
        τₚ = [Tₚ[nₜₕ] for Tₚ in 𝐓ₚ]
        return length(findall(τₚ .< T)) / length(τₚ)
    end
    
    P_τₚᵀ_q = [
        P_τₚᵀ(𝐓ₚs; T=t, ξ = 0, nₜₕ = 1)
        for t in times
    ]
    
    P_τₚᵀ_k = [
        P_τₚᵀ(𝐓ₚs; T=t, ξ = 1, nₜₕ = 1)
        for t in times
    ]
    
    ∑ = sum  # Alias for sum function
    
    function ℙ_Xₚᵀ❘ξ(Xₚᵀ, ξ; t=T, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0, nₜₕ = 1)
        if ξ == 0
            # incorrect substrate
            𝐓ₚ = 𝐓ₚs[2]
        else
            # correct substrate
            𝐓ₚ = 𝐓ₚs[1]
        end
        τₚ = [Tₚ[nₜₕ] for Tₚ in 𝐓ₚ]
        if Xₚᵀ == 1 
            return length(findall(τₚ .< t)) / length(τₚ)
        else
            return 1 - length(findall(τₚ .< t)) / length(τₚ)
        end
    end
    
    function ℙ_ξ(ξ; P_ξ₀ = 0)
        if ξ == 0
            return P_ξ₀
        else
            return 1 - P_ξ₀
        end
    end
    Ξ = [0,1]
    𝐗 = [0,1]
    
    function ℙ_Xₚᵀξ(Xₚᵀ, ξ; t=T, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0, nₜₕ = 1)
        return ℙ_Xₚᵀ❘ξ(Xₚᵀ, ξ; t=t, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ = nₜₕ) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)
    end
    
    function ℙ_Xₚᵀ(Xₚᵀ; t=T, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0, nₜₕ = 1)
        return ∑(ℙ_Xₚᵀξ(Xₚᵀ, ξ; t=t, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ=nₜₕ) for ξ ∈ Ξ)
    end
    function log₂(x)
        if x == 0 || isnan(x)
            return 0
        else
            return log2(x)
        end
    end
    
    # mutual information
    function 𝕀ₚᵀ(t,𝐓ₚs, P_ξ₀, nₜₕ)
        I = ∑(ℙ_Xₚᵀξ(X,ξ;t,𝐓ₚs,P_ξ₀,nₜₕ)* log₂(
                        ℙ_Xₚᵀξ(X,ξ;t,𝐓ₚs,P_ξ₀,nₜₕ) / (ℙ_Xₚᵀ(X;t,𝐓ₚs,P_ξ₀,nₜₕ) * ℙ_ξ(ξ;P_ξ₀))
                    ) 
            for X ∈ 𝐗, ξ ∈ Ξ)
        return I 
    end
    # channel capacity
    function CCₚᵀ(;t=T, 𝐓ₚs=𝐓ₚs, nₜₕ=1)
        # Initial guess for p₀
        initial_guess = [0.5]
        lower_bounds = [0.0]  # Lower bound for p₀
        upper_bounds = [1.0]  # Upper bound for p₀
        loss = p -> -𝕀ₚᵀ(t,𝐓ₚs,p[1], nₜₕ)
        # Perform the optimization
        opt_result = optimize(loss, lower_bounds, upper_bounds ,initial_guess, Fminbox(LBFGS()); inplace = false)
    
        # Extract optimal input distribution and channel capacity
        optimal_p₀ = opt_result.minimizer[1]
        channel_capacity = -opt_result.minimum  # Negate to get the maximum
    
        return (channel_capacity, optimal_p₀)
    end
    times_refined = [10.0^i for i in 0:0.05:6]
    𝐂ₚᵀs = [[CCₚᵀ(t=t, 𝐓ₚs=𝐓ₚs, nₜₕ=n)[1] for t in times_refined] for n in 1:100]
    for (n, 𝐂ₚᵀ) in enumerate(𝐂ₚᵀs)
        for (m, 𝐂ₚᵀ) in enumerate(𝐂ₚᵀ)
            push!(df, (Inf,k_pars.τ, times_refined[m], 𝐂ₚᵀ, n, "τₚ"))
        end
        # push!(df, (k_pars.τ, times_refined[n], 𝐂ₚᵀ, "τₚ"))
    end
    Φs_k = [result[1] for result in results_k]
    # collection of times
    Ts_k = [result[2] for result in results_k]
    # collection of first activation times
    τₐs_k = [result[3] for result in results_k]
    Ps_k = [[ϕ[i][1] for ϕ in Φs_k] for i in 1:length(times)]

    # collection of states, times and first activation times for incorrect substrate
    Φs_q = [result[1] for result in results_q]
    Ts_q = [result[2] for result in results_q]
    τₐs_q = [result[3] for result in results_q]
    η_FLD = get_η_FLD(results_k, results_q)
    Ps_q = [[ϕ[i][1] for ϕ in Φs_q] for i in 1:length(times)]
    Pₜ_k = Ps_k[end]
    Pₜ_q = Ps_q[end]
    Pₜ = [Pₜ_k, Pₜ_q]
    Ωₚ(Pₜ) = Pₜ[1] ∪ Pₜ[2] |> unique
    ℙₜ_k(Pₜ) = [length(findall(Pₜ[1] .== ω)) / length(Pₜ[1]) for ω in Ωₚ(Pₜ)]
    ℙₜ_q(Pₜ) = [length(findall(Pₜ[2] .== ω)) / length(Pₜ[2]) for ω in Ωₚ(Pₜ)]
    ℙₜ(Pₜ) = [ℙₜ_k(Pₜ), ℙₜ_q(Pₜ)]

    function ℙ_n❘ξ(n;ξ = 0, Pₜ = Pₜ, P_ξ₀ = 0)
        if ξ == 0
            return length(findall(Pₜ[2] .== n)) / length(Pₜ[2])
        else
            return length(findall(Pₜ[1] .== n)) / length(Pₜ[1])
        end
    end
    function ℙ_nξ(n,ξ; Pₜ = Pₜ, P_ξ₀ = 0)
        return ℙ_n❘ξ(n;ξ = ξ, Pₜ = Pₜ, P_ξ₀ = P_ξ₀) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)
    end

    function ℙ_n(n; Pₜ = Pₜ, P_ξ₀ = 0)
        return ∑(ℙ_nξ(n,ξ; Pₜ = Pₜ, P_ξ₀ = P_ξ₀) for ξ ∈ Ξ)
    end

    function 𝕀ₚ(Pₜ, P_ξ₀)
        I = ∑(ℙ_nξ(n,ξ; Pₜ = Pₜ, P_ξ₀ = P_ξ₀) * log₂(
                        ℙ_nξ(n,ξ; Pₜ = Pₜ, P_ξ₀ = P_ξ₀) 
                        / (ℙ_n(n; Pₜ = Pₜ, P_ξ₀ = P_ξ₀) * ℙ_ξ(ξ;P_ξ₀))
                    ) 
            for n ∈ Ωₚ(Pₜ), ξ ∈ Ξ)
        return I 
    end

    function CCₚ(;Pₜ = Pₜ)
        # Initial guess for p₀
        initial_guess = [0.5]
        lower_bounds = [0.0]  # Lower bound for p₀
        upper_bounds = [1.0]  # Upper bound for p₀
        loss = p -> -𝕀ₚ(Pₜ,p[1])
        # Perform the optimization
        opt_result = optimize(loss, lower_bounds, upper_bounds ,initial_guess, Fminbox(LBFGS()); inplace = false)

        # Extract optimal input distribution and channel capacity
        optimal_p₀ = opt_result.minimizer[1]
        channel_capacity = -opt_result.minimum  # Negate to get the maximum

        return (channel_capacity, optimal_p₀)
    end


    𝐈ₚ = [𝕀ₚ([Ps_k[i], Ps_q[i]], 0.5) for i ∈ eachindex(times)]
    𝐂ₚ = [CCₚ(Pₜ = [Ps_k[i], Ps_q[i]])[1] for i ∈ eachindex(times)]
    for (n, 𝐂ₚ) in enumerate(𝐂ₚ)
        push!(df, (Inf, k_pars.τ, times[n], 𝐂ₚ, 0, "P"))
    end

    # collection of states
    Φs_k = [result[1] for result in results_k]
    # collection of times
    Ts_k = [result[2] for result in results_k]
    # collection of first activation times
    τₐs_k = [result[3] for result in results_k]
    Ps_k = [[ϕ[i][1] for ϕ in Φs_k] for i in 1:length(times)]

    # collection of states, times and first activation times for incorrect substrate
    Φs_q = [result[1] for result in results_q]
    Ts_q = [result[2] for result in results_q]
    τₐs_q = [result[3] for result in results_q]
    Ps_q = [[ϕ[i][1] for ϕ in Φs_q] for i in 1:length(times)]

    τₐs = [τₐs_k, τₐs_q]
    function P_τₐ(τₐs; T=1e6, ξ = 0) # the probability of τₐ < T
        if ξ == 0
            # incorrect substrate
            τₐ = τₐs[2]
        else
            # correct substrate
            τₐ = τₐs[1]
        end
        return length(findall(τₐ .< T)) / length(τₐ)
    end
    P_τₐs_q = [
        P_τₐ(τₐs; T=t, ξ = 0)
        for t in times
    ]

    P_τₐs_k = [
        P_τₐ(τₐs; T=t, ξ = 1)
        for t in times
    ]

    # For fixed t, and probability of P_ξ₀ (Probability of incorrect substrate),
    # compute the mutual information I(ξ; Xₐ ≡ τₐ < t | t)
    ∑ = sum  # Alias for sum function

    function ℙ_Xₐ❘ξ(Xₐ, ξ; t=T, τₐs=τₐs, P_ξ₀ = 0)
        if ξ == 0
            # incorrect substrate
            τₐ = τₐs[2]
        else
            # correct substrate
            τₐ = τₐs[1]
        end
        if Xₐ == 1 
            return length(findall(τₐ .< t)) / length(τₐ)
        else
            return 1 - length(findall(τₐ .< t)) / length(τₐ)
        end
    end


    Ξ = [0,1]
    𝐗 = [0,1]
    function ℙ_Xₐξ(Xₐ, ξ; t=T, τₐs=τₐs, P_ξ₀ = 0)
        return ℙ_Xₐ❘ξ(Xₐ, ξ; t=t, τₐs=τₐs, P_ξ₀ = P_ξ₀) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)
    end

    function ℙ_Xₐ(Xₐ; t=T, τₐs=τₐs, P_ξ₀ = 0)
        return ∑(ℙ_Xₐξ(Xₐ, ξ; t=t, τₐs=τₐs, P_ξ₀ = P_ξ₀) for ξ ∈ Ξ)
    end

    # mutual information
    function 𝕀ₐ(t,τₐs, P_ξ₀)
        I = ∑(ℙ_Xₐξ(X,ξ;t,τₐs,P_ξ₀)* log₂(
                        ℙ_Xₐξ(X,ξ;t,τₐs,P_ξ₀) / (ℙ_Xₐ(X;t,τₐs,P_ξ₀) * ℙ_ξ(ξ;P_ξ₀))
                    ) 
            for X ∈ 𝐗, ξ ∈ Ξ)
        return I 
    end
    # channel capacity
    function CCₐ(;t=T, τₐs=τₐs)
        # Initial guess for p₀
        initial_guess = [0.5]
        lower_bounds = [0.0]  # Lower bound for p₀
        upper_bounds = [1.0]  # Upper bound for p₀
        loss = p -> -𝕀ₐ(t,τₐs,p[1])
        # Perform the optimization
        opt_result = optimize(loss, lower_bounds, upper_bounds ,initial_guess, Fminbox(LBFGS()); inplace = false)

        # Extract optimal input distribution and channel capacity
        optimal_p₀ = opt_result.minimizer[1]
        channel_capacity = -opt_result.minimum  # Negate to get the maximum

        return (channel_capacity, optimal_p₀)
    end

    CCₐs = [CCₐ(t=t, τₐs=τₐs) for t in times]
    @show CCₐs
    for (n, CCₐ) in enumerate(CCₐs)
        push!(df, (Inf, k_pars.τ, times[n], CCₐs[n][1], 0,"τₐ"))
    end
    for (t, η) in zip(times, η_FLD)
        push!(df, (Inf, k_pars.τ, t, η, 0, "η_FLD"))
    end
end

df = DataFrame(
    N_f = Float64[],
    τ = Float64[],
    T = Float64[],
    CC = Float64[],
    nₜₕ = Int[],
    type = String[],
)
for (results_k, results_q, profile) in zip(Results_k_∞, Results_q_∞, profiles_∞)
    handle_results_∞(results_k, results_q, profile)
end

df

# append the results to the dataframe 
CSV.write("DATA/multiround_multistep_simulation_product_CC_summary_$(group_label).csv", df, append = true)