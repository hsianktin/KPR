@info "This script assumes simulation results are already available in DATA/"
using BSON
using StatsBase

using ProgressMeter
profile = "PROFILES/base.jl"
tasks = [(substrate, id) for substrate in ["k", "q"] for id in 1:100]
label = split(profile, "/")[end] |> x -> split(x, ".jl")[1]

# collect results
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

include("multiround_simulation_utils.jl");
include(profile);
𝐓ₐs_k = [result[3] for result in results_k]
𝐓ₐs_q = [result[3] for result in results_q]
𝐓ₚs_k =  [result[4] for result in results_k]
𝐓ₚs_q =  [result[4] for result in results_q]
𝐓ₐs = [𝐓ₐs_k, 𝐓ₐs_q]
𝐓ₚs = [𝐓ₚs_k, 𝐓ₚs_q]

function P_τₐᵀ(𝐓ₐs; T=1e6, ξ = 0)
    if ξ == 0
        # incorrect substrate
        𝐓ₐ = 𝐓ₐs[2]
    else
        # correct substrate
        𝐓ₐ = 𝐓ₐs[1]
    end
    τₐ = [Tₐ for Tₐ in 𝐓ₐ]
    return length(findall(τₐ .< T)) / length(τₐ)
end

P_τₐᵀ_q = [
    P_τₐᵀ(𝐓ₐs; T=t, ξ = 0)
    for t in times
]

P_τₐᵀ_k = [
    P_τₐᵀ(𝐓ₐs; T=t, ξ = 1)
    for t in times
]

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

function ℙ_Xₐᵀ❘ξ(Xₐᵀ, ξ; t=T, 𝐓ₐs=𝐓ₐs, P_ξ₀ = 0)
    if ξ == 0
        # incorrect substrate
        𝐓ₐ = 𝐓ₐs[2]
    else
        # correct substrate
        𝐓ₐ = 𝐓ₐs[1]
    end
    τₐ = [Tₐ for Tₐ in 𝐓ₐ]
    if Xₐᵀ == 1 
        return length(findall(τₐ .< t)) / length(τₐ)
    else
        return 1 - length(findall(τₐ .< t)) / length(τₐ)
    end
end

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

function ℙ_Xₐᵀξ(Xₐᵀ, ξ; t=T, 𝐓ₐs=𝐓ₐs, P_ξ₀ = 0)
    return ℙ_Xₐᵀ❘ξ(Xₐᵀ, ξ; t=t, 𝐓ₐs=𝐓ₐs, P_ξ₀ = P_ξ₀) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)
end

function ℙ_Xₐᵀ(Xₐᵀ; t=T, 𝐓ₐs=𝐓ₐs, P_ξ₀ = 0)
    return ∑(ℙ_Xₐᵀξ(Xₐᵀ, ξ; t=t, 𝐓ₐs=𝐓ₐs, P_ξ₀ = P_ξ₀) for ξ ∈ Ξ)
end

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
function 𝕀ₐᵀ(t,𝐓ₐs, P_ξ₀)
    I = ∑(ℙ_Xₐᵀξ(X,ξ;t,𝐓ₐs,P_ξ₀)* log₂(
                    ℙ_Xₐᵀξ(X,ξ;t,𝐓ₐs,P_ξ₀) / (ℙ_Xₐᵀ(X;t,𝐓ₐs,P_ξ₀) * ℙ_ξ(ξ;P_ξ₀))
                ) 
        for X ∈ 𝐗, ξ ∈ Ξ)
    return I 
end

function 𝕀ₚᵀ(t,𝐓ₚs, P_ξ₀, nₜₕ)
    I = ∑(ℙ_Xₚᵀξ(X,ξ;t,𝐓ₚs,P_ξ₀,nₜₕ)* log₂(
                    ℙ_Xₚᵀξ(X,ξ;t,𝐓ₚs,P_ξ₀,nₜₕ) / (ℙ_Xₚᵀ(X;t,𝐓ₚs,P_ξ₀,nₜₕ) * ℙ_ξ(ξ;P_ξ₀))
                ) 
        for X ∈ 𝐗, ξ ∈ Ξ)
    return I 
end
# channel capacity
using Optim
function CCₐᵀ(;t=T, 𝐓ₐs=𝐓ₐs)
    # Initial guess for p₀
    initial_guess = [0.5]
    lower_bounds = [0.0]  # Lower bound for p₀
    upper_bounds = [1.0]  # Upper bound for p₀
    loss = p -> -𝕀ₐᵀ(t,𝐓ₐs,p[1])
    # Perform the optimization
    opt_result = optimize(loss, lower_bounds, upper_bounds ,initial_guess, Fminbox(LBFGS()); inplace = false)

    # Extract optimal input distribution and channel capacity
    optimal_p₀ = opt_result.minimizer[1]
    channel_capacity = -opt_result.minimum  # Negate to get the maximum

    return (channel_capacity, optimal_p₀)
end

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
𝐂ₐᵀ = [CCₐᵀ(t=t, 𝐓ₐs=𝐓ₐs)[1] for t in times_refined]
𝐂ₚᵀs = [[CCₚᵀ(t=t, 𝐓ₚs=𝐓ₚs, nₜₕ=n)[1] for t in times_refined] for n in 1:100]

using PGFPlotsX


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


plot0 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="cell contact time \$T\$",
        ylabel="channel capacity \$C\$",
        legend_pos="north west",
        xmode="log",
        legend_cell_align="left",
        legend_style={fill="none", draw="none"},
    },
    Plot({mark="*",color="red", mark_repeat=2}, Table(times_refined, 𝐂ₐᵀ)),    
    Plot({mark="square",color="blue"}, Table(times, 𝐂ₚ)),
    LegendEntry("\$C\\big(\\xi; X_{\\rm a}\\big)\$"),
    LegendEntry("\$C\\big(\\xi; {\\rm P}(T)\\big)\$"),
)

pgfsave("FIG/multiround_simulation_fpt_analysis_$(label).pgf", plot0)
pgfsave("FIG/multiround_simulation_fpt_analysis_$(label).pdf", plot0)
pgfsave("FIG/multiround_simulation_fpt_analysis_$(label).svg", plot0)
# call pdftops to convert pdf to eps
run(`pdftops -eps FIG/multiround_simulation_fpt_analysis_$(label).pdf`)
@info "Figure saved to FIG/multiround_simulation_fpt_analysis_$(label).pdf"

