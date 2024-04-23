using Distributed
using BSON
using StatsBase
using ArgParse
addprocs(10)
# addprocs([("ryzen.lan",8)])
@everywhere begin 
    using ProgressMeter
end
# @everywhere profile = "PROFILES/large_kp.jl"

# Parse command line arguments
s = ArgParseSettings(description = "provide profile files")
@add_arg_table! s begin
    "--profile"
        required=false
        default = "PROFILES/base.jl"
end

args = parse_args(s)
profile = args["profile"]
label = split(profile, "/")[end] |> x -> split(x, ".jl")[1]

tasks = [(substrate, id, profile) for substrate in ["k", "q"] for id in 1:100]
@everywhere function execute_task(task)
    substrate, id, profile = task
    cmd = `julia multiround_simulation_single_thread.jl $substrate $id $profile`
    run(cmd)
end

progress_pmap(
    x -> execute_task(x),
    tasks;
    progress=Progress(
        length(tasks);
        desc="Simulating",
        color=:white,
        showspeed=true,
        barglyphs=BarGlyphs("[=> ]")
    )
)

rmprocs(workers())
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

function ℙ_ξ(ξ; P_ξ₀ = 0)
    if ξ == 0
        return P_ξ₀
    else
        return 1 - P_ξ₀
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
function log₂(x)
    if x == 0 || isnan(x)
        return 0
    else
        return log2(x)
    end
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
using Optim
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


### Product-based discrimination and Channel Capacity  ####
### Can only calculate for a specific time point t. 
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


### Visualization ####
using PGFPlotsX

# survival probability
plot1 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="cell contact time T",
        ylabel="survival probability",
        legend_pos="north west",
        xmode="log"
    },
    PlotInc({mark="o",color="green!50!black"}, Table(times[2:end], P_τₐs_q[2:end])),
    PlotInc({mark="square",color="yellow!50!black"}, Table(times[2:end], P_τₐs_k[2:end])),
    LegendEntry("self ligand"),
    LegendEntry("mutant ligand"),
)
pgfsave("DATA/survival_probability_simulation_$(label).pgf", plot1)
pgfsave("DATA/survival_probability_simulation_$(label).pdf", plot1)

# mutual information and channel capacity are nearly identical with p = 0.5
𝕀ₐs = [𝕀ₐ(t,τₐs,0.5) for t in times]
CCₐs = [CCₐ(t=t, τₐs=τₐs) for t in times]
plot2 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="cell contact time T",
        ylabel="mutual information and channel capacity",
        legend_pos="north west",
        xmode="log"
    },
    Plot({mark="o",color="blue",only_marks}, Table(times[2:end], 𝕀ₐs[2:end])),
    Plot({mark="*",color="blue",dashed}, Table(times[2:end], [CCₐs[i][1] for i in 2:length(times)])),
    Plot({mark="square",color="red",only_marks}, Table(times[2:end], 𝐈ₚ[2:end])),
    Plot({mark="square*",color="red",dashed}, Table(times[2:end], 𝐂ₚ[2:end])),
    LegendEntry("𝕀ₐ"),
    LegendEntry("CCₐ"),
    LegendEntry("𝕀ₚ"),
    LegendEntry("CCₚ"),
)
pgfsave("DATA/mutual_information_and_channel_capacity_simulation_$(label).pgf", plot2)
pgfsave("DATA/mutual_information_and_channel_capacity_simulation_$(label).pdf", plot2)

using DataFrames
df = DataFrame(
    t = times[2:end],
    mutual_information_a = 𝕀ₐs[2:end],
    CC_a = [CCₐs[i][1] for i in 2:length(times)],
    mutual_information_p = 𝐈ₚ[2:end],
    CC_p = 𝐂ₚ[2:end]
)
using CSV
CSV.write("DATA/numerical_mutual_information_and_channel_capacity_$(label).csv", df)
