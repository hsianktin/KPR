@info "This script assumes simulation results are already available in DATA/"
@info """
The purpose of this script is to obtain 
the channel capacity of the product/substrate channels.
- CCₚ is the channel capacity of the product number distribution nₜ | input ξ channel
- CCₚᵀ is the channel capacity of the FPT of product reaching nₜₕ < T | input ξ channel
"""
using BSON
using StatsBase

using ProgressMeter
profile = "PROFILES/large_kp.jl"
tasks = [(substrate, id) for substrate in ["k", "q"] for id in 1:99]
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

Ns_k = [[ϕ[1] for ϕ in Φ] for Φ in Φs_k] # N_repeat samples, each sample is a time trajectory
Nt_k = [[Ns_k[i][t] for i in eachindex(Ns_k)] for t in eachindex(Ts_k[1])]
𝔼nₜ❘k = [mean(Nt) for Nt in Nt_k]
𝕍nₜ❘k = [var(Nt) for Nt in Nt_k]

Ns_q = [[ϕ[1] for ϕ in Φ] for Φ in Φs_q] # N_repeat samples, each sample is a time trajectory
Nt_q = [[Ns_q[i][t] for i in eachindex(Ns_q)] for t in eachindex(Ts_q[1])]
𝔼nₜ❘q = [mean(Nt) for Nt in Nt_q]
𝕍nₜ❘q = [var(Nt) for Nt in Nt_q]

time = Ts_k[1]

# η_FLD = |𝔼nₜ∣q - 𝔼nₜ∣k|/(√𝕍nₜ∣q + √𝕍nₜ∣k)
η_FLD = [abs(𝔼nₜ❘q[i] - 𝔼nₜ❘k[i]) / (sqrt(𝕍nₜ❘q[i]) + sqrt(𝕍nₜ❘k[i]))     for i in eachindex(time)]
# replace NaN with 0
η_FLD = [ifelse(isnan(x), 0, x) for x in η_FLD]


𝐓ₚs_k =  [result[4] for result in results_k]
𝐓ₚs_q =  [result[4] for result in results_q]
𝐓ₚs = [𝐓ₚs_k, 𝐓ₚs_q]

##################################################################
############## Additional Ploting for Illustration ###############
##################################################################
using CSV, DataFrames
using PGFPlotsX
tikzpicture = @pgf TikzPicture({baseline})
axis = @pgf Axis({
    width="3in",
    height="3in",
    xlabel="contact time \$T\$",
    ylabel="number of products \${\\rm P}_T\$",
    xmode="log",
    ymode="log",
    xmin=1,
    ymin=1,
    # shift y label to right
    ylabel_style={yshift="-0.1cm"},
    "after end axis/.code" = {
        raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(a)};"
    },
    legend_cell_align="left",
    # legend_pos="north west",
    legend_style={
        at="{(0.15,0.95)}",
        anchor="north west",
        fill="none",
        draw="none",
    }
})
# plot()
counter=0
for Nk in Ns_k
    global counter
    if counter < 20
        # filter zeros
        times = time[Nk .> 0]
        Nk = Nk[Nk .> 0]
        # plot!(times, Nk, color="blue", alpha=0.1)
        
        if counter == 1
            push!(axis, @pgf Plot({mark="square",color="blue!20!white"}, 
                Table(times, Nk))
            )
        else
            push!(axis, @pgf Plot({mark="square",color="blue!20!white",forget_plot=true}, 
                Table(times, Nk))
            )
        end
        
        # save to csv as df 
        df = DataFrame(time=times, N=Nk)
        CSV.write("DATA/multiround_simulation_product_analysis_Nk_$(counter).csv", df)
    end
    counter += 1
end
axis


counter=0
for Nq in Ns_q
    global counter
    if counter < 20
        # filter zeros
        times = time[Nq .> 0]
        Nq = Nq[Nq .> 0]
        if counter == 1
            push!(axis, @pgf Plot({mark="o",color="red!20!white"}, 
                Table(times, Nq))
            )
        else
            push!(axis, @pgf Plot({mark="o",color="red!20!white",forget_plot=true}, 
                Table(times, Nq))
            )
        end
        df = DataFrame(time=times, N=Nq)
        CSV.write("DATA/multiround_simulation_product_analysis_Nq_$(counter).csv", df)
    end
    counter += 1
end
axis
# xlabel!("time")
# ylabel!("product number")




# implements the optimal threshold
function Pₜₕꜛ(T, pars_k, pars_q)
    K = pars_k.k₁/(pars_k.k₁ + pars_k.k₋₁)
    Q = pars_q.k₁/(pars_q.k₁ + pars_q.k₋₁)
    τ = pars_k.τ
    return pars_k.kₚ * T * √(K * Q * exp(-pars_k.k₋₁ *  τ - pars_q.k₋₁ * τ))
end
# plot the threshold
# plot!(time, [Pₜₕꜛ(t, k_pars, q_pars) for t in time], color="black", linestyle=:dash, label="threshold")
# save to csv as df
df = DataFrame(time=time, threshold=[Pₜₕꜛ(t, k_pars, q_pars) for t in time])
# filter df such that threshold > 1
df = df[df.threshold .> 1, :]
push!(axis, @pgf Plot({mark="none", color="black", dashed},
    Table(df.time, df.threshold))
)
# also plot a constant line where the threshold is equal to the end of df.threshold
push!(axis, @pgf Plot({mark="none", color="purple", dashed,thick}, Coordinates([(1,200), (1e6,200)]))
)

axis
CSV.write("DATA/multiround_simulation_product_analysis_threshold.csv", df)

push!(axis, LegendEntry("correct substrate"), 
)
push!(axis, LegendEntry("incorrect substrate"), 
)
push!(axis, LegendEntry("dynamic threshold"), 
)
push!(axis, LegendEntry("static threshold"), 
)

push!(tikzpicture, axis)
pgfsave("FIG/dynamical_threshold.pgf", tikzpicture)
pgfsave("FIG/dynamical_threshold.pdf", tikzpicture)
@info  "FIG/dynamical_threshold.pdf saved"



#### Plot of Xₐ
function stopping_time_to_trajectory(τₐ, time)
    return [t < τₐ ? 0 : 1 for t in time]
end

function stopping_time_to_trajectories(τₐs, time)
    return [stopping_time_to_trajectory(τₐ, time) for τₐ in τₐs]
end

tikzpicture = @pgf TikzPicture({baseline})
axis = @pgf Axis({
    width="3in",
    height="3in",
    xlabel="time since initial contact \$t\$",
    ylabel="activation status \$X_{\\rm a}(t)\$",
    xmode="log",
    # ymode="log",
    xmin=1,
    # ymin=1,
    ymin=-0.2,
    ymax=1.6,
    # yticks only at 0, 1
    ytick="{0,1}",
    # shift y label to right
    ylabel_style={yshift="-0.1cm"},
    "after end axis/.code" = {
        raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(b)};"
    },
    legend_cell_align="left",
    legend_pos="north west",
    legend_style={
        at="{(0.15,0.95)}",
        anchor="north west",
        fill="none",
        draw="none",
    }
})

# k substrate
counter=0
for τₐ in τₐs_k
    global counter
    if counter < 20
        # filter zeros
        # times = time[τₐ .> 0]
        # τₐ = τₐ[τₐ .> 0]
        if counter == 1
            push!(axis, @pgf Plot({mark="square",color="blue!40!white"}, 
                Table(time, stopping_time_to_trajectory(τₐ, time))
            )
            )
        else
            push!(axis, @pgf Plot({mark="square",color="blue!40!white",forget_plot=true}, 
                Table(time, stopping_time_to_trajectory(τₐ, time))
            )
            )
        end
        # df = DataFrame(time=times, X=stopping_time_to_trajectory(τₐ, time))
        # CSV.write("DATA/multiround_simulation_product_analysis_τₐ_k_$(counter).csv", df)
    end
    counter += 1
end
# q substrate
counter=0
for τₐ in τₐs_q
    global counter
    if counter < 20
        # filter zeros
        # times = time[τₐ .> 0]
        # τₐ = τₐ[τₐ .> 0]
        if counter == 1
            push!(axis, @pgf Plot({mark="o",color="red!40!white"}, 
                Table(time, stopping_time_to_trajectory(τₐ, time))
            )
            )
        else
            push!(axis, @pgf Plot({mark="o",color="red!40!white",forget_plot=true}, 
                Table(time, stopping_time_to_trajectory(τₐ, time))
            )
            )
        end
        # df = DataFrame(time=times, X=stopping_time_to_trajectory(τₐ, time))
        # CSV.write("DATA/multiround_simulation_product_analysis_τₐ_q_$(counter).csv", df)
    end
    counter += 1
end

# draw a vertical line at t = 400 in axis
push!(axis, @pgf Plot({mark="none", color="black", dashed, thick}, Coordinates([(300,0), (300,1)]))
)

# legend
push!(axis, LegendEntry("correct substrate"), 
)
push!(axis, LegendEntry("incorrect substrate"), 
)

push!(axis, LegendEntry("cell contact time \$ T\$"))

axis
push!(tikzpicture, axis)
pgfsave("FIG/FPT_schematic.pgf", tikzpicture)
pgfsave("FIG/FPT_schematic.pdf", tikzpicture)
@info  "FIG/FPT_schematic.pdf saved"

# histogram of P at T 
NT_k = [Ns_k[i][24] for i in eachindex(Ns_k)]
NT_q = [Ns_q[i][24] for i in eachindex(Ns_q)]

# histogram of P at T using PGFPlotsX
tikzpicture = @pgf TikzPicture({baseline})
axis = @pgf Axis({
    width="3in",
    height="3in",
    xlabel="number of products \${\\rm P}_T\$",
    ylabel="frequency",
    xmin=0,
    ymin=0,
    ybar,
    # xmode="log",
    # no y ticks 
    ytick="{}",
    yticklabels="{}",
    # shift y label to right
    ylabel_style={yshift="-0.1cm"},
    # y ticks using scientific notation
    # yticklabel_style = "{/pgf/number format/sci}",  # Use scientific notation
    "after end axis/.code" = {
        raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(c)};"
    },
    legend_cell_align="left",
    legend_pos="north west",
    legend_style={
        at="{(0.15,0.95)}",
        anchor="north west",
        fill="none",
        draw="none",
    }
})
# k substrate
# determine the bin counts
# load packages
using StatsBase
hist = fit(Histogram, NT_k, nbins=30)
# convert to bin position and bin count
bin_positions = hist.edges[1][1:end-1]
bin_counts = hist.weights
# normalize the bin counts 
bin_counts = bin_counts ./ sum(bin_counts)
# weight the bin counts by bin width
bin_counts = bin_counts ./ diff(hist.edges[1])
# plot the histogram in PGFPlotsX and adjust width of the bars
push!(axis, @pgf Plot({ybar, bar_width="0.2cm", color="blue!80!white", fill="blue!80!white"}, 
    Table(bin_positions, bin_counts))
)

# q substrate
hist = fit(Histogram, NT_q, nbins=30)
bin_positions = hist.edges[1][1:end-1]
bin_counts = hist.weights
bin_counts = bin_counts ./ sum(bin_counts)
bin_counts = bin_counts ./ diff(hist.edges[1])
push!(axis, @pgf Plot({ybar, bar_width="0.05cm", color="red!80!white", fill="red!80!white"}, 
    Table(bin_positions, bin_counts))
)
push!(axis, LegendEntry("correct substrate"), 
)
push!(axis, LegendEntry("incorrect substrate"), 
)

push!(tikzpicture, axis)
pgfsave("FIG/histogram_of_P.pgf", tikzpicture)
pgfsave("FIG/histogram_of_P.pdf", tikzpicture)
@info  "FIG/histogram_of_P.pdf saved"


##################################################################
function Xₜ(𝐓ₚ, pars_k, pars_q)
    # reference time based on Pₜₕꜛ
    # find the time T such that Pₜₕꜛ(T, pars_k, pars_q) = 1, 2, 3,..., 100
    K = pars_k.k₁/(pars_k.k₁ + pars_k.k₋₁)
    Q = pars_q.k₁/(pars_q.k₁ + pars_q.k₋₁)
    τ = pars_k.τ
    𝐓ₜₕ = [
        Pₜₕ / (pars_k.kₚ * √(K * Q * exp(-pars_k.k₋₁ *  τ - pars_q.k₋₁ * τ))) 
        for Pₜₕ in 1:100
    ]
    response = 0
    for (Tₚ, Tₜₕ) in zip(𝐓ₚ[10:end], 𝐓ₜₕ[10:end])
        if Tₚ < Tₜₕ # the number of products exceeds the threshold at Tₜₕ
            response = 1
            break
        end
    end
    return response 
end
function Tths(pars_k, pars_q)
    # reference time based on Pₜₕꜛ
    # find the time T such that Pₜₕꜛ(T, pars_k, pars_q) = 1, 2, 3,..., 100
    K = pars_k.k₁/(pars_k.k₁ + pars_k.k₋₁)
    Q = pars_q.k₁/(pars_q.k₁ + pars_q.k₋₁)
    τ = pars_k.τ
    𝐓ₜₕ = [
        Pₜₕ / (pars_k.kₚ * √(K * Q * exp(-pars_k.k₋₁ *  τ - pars_q.k₋₁ * τ))) 
        for Pₜₕ in 1:100
    ]
    return 𝐓ₜₕ
end


∑ = sum  # Alias for sum function

function P_τₚᵀₜₕ(𝐓ₚs;  ξ = 0, pars_k = k_pars, pars_q=q_pars) # the probability of τₚᵀ < T
    if ξ == 0
        # incorrect substrate
        𝐓ₚs = 𝐓ₚs[2]
    else
        # correct substrate
        𝐓ₚs = 𝐓ₚs[1]
    end
    𝐗ₜ = [Xₜ(𝐓ₚ, pars_k, pars_q) for 𝐓ₚ in 𝐓ₚs]
    return ∑(𝐗ₜ)/length(𝐗ₜ)
end

function Xₜ_randT(𝐓ₚ, pars_k, pars_q,t)
        # reference time based on Pₜₕꜛ
    # find the time T such that Pₜₕꜛ(T, pars_k, pars_q) = 1, 2, 3,..., 100
    K = pars_k.k₁/(pars_k.k₁ + pars_k.k₋₁)
    Q = pars_q.k₁/(pars_q.k₁ + pars_q.k₋₁)
    τ = pars_k.τ
    𝐓ₜₕ = [
        Pₜₕ / (pars_k.kₚ * √(K * Q * exp(-pars_k.k₋₁ *  τ - pars_q.k₋₁ * τ))) 
        for Pₜₕ in 1:100
    ]
    response = 0
    for (Tₚ, Tₜₕ) in zip(𝐓ₚ[10:end], 𝐓ₜₕ[10:end])
        if Tₚ < Tₜₕ && Tₚ < t # the number of products exceeds the threshold at Tₜₕ
            response = 1
            break
        end
    end
    return response
end


function P_τₚᵀₜₕ_randT(𝐓ₚs;  ξ = 0, pars_k = k_pars, pars_q=q_pars, Ts_k=rand_Ts_k, Ts_q=rand_Ts_q) # the probability of τₚᵀ < T
    if ξ == 0
        # incorrect substrate
        𝐓ₚs = 𝐓ₚs[2]
        Ts = Ts_q
    else
        # correct substrate
        𝐓ₚs = 𝐓ₚs[1]
        Ts = Ts_k
    end
    𝐗ₜ = [Xₜ_randT(𝐓ₚ, pars_k, pars_q, T) for (𝐓ₚ, T) in zip(𝐓ₚs, Ts)]
    return ∑(𝐗ₜ)/length(𝐗ₜ)
end


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

using Plots
begin
i=11
scatter(𝐓ₚs_q[i], label="sample 1")
scatter!(Tths(k_pars, q_pars), label="threshold")
title!("$(Xₜ(𝐓ₚs_q[i], k_pars, q_pars))")
end
P_τₚᵀₜₕ_q =    P_τₚᵀₜₕ(𝐓ₚs; ξ = 0)
P_τₚᵀₜₕ_k =    P_τₚᵀₜₕ(𝐓ₚs; ξ = 1)

function ℙ_Xₚ❘ξ(X, ξ; 𝐓ₚs=𝐓ₚs , pars_k = k_pars, pars_q=q_pars)
    if X == 1
        return P_τₚᵀₜₕ(𝐓ₚs; ξ = ξ, pars_k = pars_k, pars_q = pars_q)
    else
        return 1 - P_τₚᵀₜₕ(𝐓ₚs; ξ = ξ, pars_k = pars_k, pars_q = pars_q)
    end
end
function ℙ_Xₚ❘ξ_randT(X, ξ; Ts_k = rand_Ts_k, Ts_q=rand_Ts_q, 𝐓ₚs=𝐓ₚs , pars_k = k_pars, pars_q=q_pars)
    if X == 1
        return P_τₚᵀₜₕ_randT(𝐓ₚs; ξ = ξ, pars_k = pars_k, pars_q = pars_q, Ts_k = Ts_k, Ts_q = Ts_q)
    else
        return 1 - P_τₚᵀₜₕ_randT(𝐓ₚs; ξ = ξ, pars_k = pars_k, pars_q = pars_q, Ts_k = Ts_k, Ts_q = Ts_q)
    end
end

T = k_pars.T

rand_Ts_k = [rand() * T for i in 𝐓ₚs_k]
rand_Ts_q = [rand() * T for i in 𝐓ₚs_q]

Ts_k = rand_Ts_k
Ts_q = rand_Ts_q
function ℙ_Xₚᵀ❘ξ_randT(Xₚᵀ, ξ; Ts_k = rand_Ts_k, Ts_q=rand_Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0, nₜₕ = 1)
    if ξ == 0
        # incorrect substrate
        𝐓ₚ = 𝐓ₚs[2]
        Ts = Ts_q
    else
        # correct substrate
        𝐓ₚ = 𝐓ₚs[1]
        Ts = Ts_k
    end
    result = [Tₚ[nₜₕ] < t for (Tₚ,t) in zip(𝐓ₚ, Ts)]
    if Xₚᵀ == 1 
        return ∑(result)/length(result)
    else
        return 1 - ∑(result)/length(result)
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

function ℙ_Xₚᵀξ(Xₚᵀ, ξ; t=T, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0, nₜₕ = 1) # joint probability 
    return ℙ_Xₚᵀ❘ξ(Xₚᵀ, ξ; t=t, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ = nₜₕ) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)
end

function ℙ_Xₚᵀξ_randT(Xₚᵀ, ξ; Ts_k = rand_Ts_k, Ts_q=rand_Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0.5, nₜₕ = 1) # joint probability 
    return ℙ_Xₚᵀ❘ξ_randT(Xₚᵀ, ξ; Ts_k = Ts_k, Ts_q = Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ = nₜₕ) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)
end

function ℙ_Xₚᵀ(Xₚᵀ; t=T, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0, nₜₕ = 1)
    return ∑(ℙ_Xₚᵀξ(Xₚᵀ, ξ; t=t, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ=nₜₕ) for ξ ∈ Ξ)
end

function ℙ_Xₚᵀ_randT(Xₚᵀ; Ts_k = rand_Ts_k, Ts_q=rand_Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0, nₜₕ = 1)
    return ∑(ℙ_Xₚᵀξ_randT(Xₚᵀ, ξ; Ts_k = Ts_k, Ts_q = Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ=nₜₕ) for ξ ∈ Ξ)
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

function 𝕀ₚᵀ_randT(;t=T, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0.5, nₜₕ = 1, Ts_k=rand_Ts_k, Ts_q=rand_Ts_q)
    I = ∑(ℙ_Xₚᵀξ_randT(Xₚᵀ, ξ; Ts_k = Ts_k, Ts_q = Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ=nₜₕ) * log₂(
        ℙ_Xₚᵀξ_randT(Xₚᵀ, ξ; Ts_k = Ts_k, Ts_q=Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ = 1)
                    / (ℙ_Xₚᵀ_randT(Xₚᵀ; Ts_k = Ts_k, Ts_q=Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = P_ξ₀, nₜₕ = 1)* ℙ_ξ(ξ;P_ξ₀))
                ) 
        for Xₚᵀ ∈ 𝐗, ξ ∈ Ξ)
    return I 
end

𝕀ₚᵀ(1e6, 𝐓ₚs, 0.5, 10)
ℙ_Xₚᵀξ_randT(1,1,; nₜₕ=1, )
ℙ_Xₚᵀξ_randT(1, 1; Ts_k = rand_Ts_k, Ts_q = rand_Ts_q, 𝐓ₚs=𝐓ₚs, P_ξ₀ = 0.5, nₜₕ=15)
function 𝕀ₚ(𝐓ₚs, P_ξ₀)
    I = ∑(ℙ_Xₚ❘ξ(X,ξ; 𝐓ₚs=𝐓ₚs) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)  * log₂(
        ℙ_Xₚ❘ξ(X,ξ; 𝐓ₚs=𝐓ₚs) / ∑(
            [
                ℙ_Xₚ❘ξ(X,ξ′; 𝐓ₚs=𝐓ₚs) * ℙ_ξ(ξ′; P_ξ₀ = P_ξ₀) 
                for ξ′ ∈ Ξ
            ]
        )
    ) for ξ ∈ Ξ, X ∈ 𝐗)
    return I
end

function 𝕀ₚ_randT(𝐓ₚs, P_ξ₀, Ts_k, Ts_q)
    I = ∑(ℙ_Xₚ❘ξ_randT(X,ξ; Ts_k = Ts_k, Ts_q = Ts_q, 𝐓ₚs=𝐓ₚs) * ℙ_ξ(ξ; P_ξ₀ = P_ξ₀)  * log₂(
        ℙ_Xₚ❘ξ_randT(X,ξ; Ts_k = Ts_k, Ts_q = Ts_q, 𝐓ₚs=𝐓ₚs) / ∑(
            [
                ℙ_Xₚ❘ξ_randT(X,ξ′; Ts_k = Ts_k, Ts_q = Ts_q, 𝐓ₚs=𝐓ₚs) * ℙ_ξ(ξ′; P_ξ₀ = P_ξ₀) 
                for ξ′ ∈ Ξ
            ]
        )
    ) for ξ ∈ Ξ, X ∈ 𝐗)
    return I
end

push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usetikzlibrary{patterns}")
mutual_information_of_dynamical_threshold_fixed_T = 𝕀ₚ(𝐓ₚs, 0.5)
mutual_information_of_dynamical_threshold_random_T = 𝕀ₚ_randT(𝐓ₚs, 0.5, rand_Ts_k, rand_Ts_q)
# now consider the possibility of randomly masking the total time T
Ts = [10.0^i for i in 0:0.2:6]
maximum_mutual_information_of_fixed_threshold_fixed_T = maximum([𝕀ₚᵀ(T, 𝐓ₚs, 0.5, nₜₕ) for nₜₕ in 1:100, T in Ts])
maximum_mutual_information_of_fixed_threshold_random_T = maximum([𝕀ₚᵀ_randT(𝐓ₚs=𝐓ₚs, Ts_k=rand_Ts_k, Ts_q = rand_Ts_q, nₜₕ=n,P_ξ₀=0.5) for n in 1:100])

# use a grouped bar plots to visualize the above 4 quantities in PGFPlotsX 
using PGFPlotsX
tikzpicture = @pgf TikzPicture({baseline})
axis = @pgf Axis({
    width="3in",
    height="3in",
    ybar,
    ylabel= "mutual information \$\\mathcal{I}\$",
    xtick="{0.93,1.27}",
    xticklabels={
        "static threshold",
        "dynamic threshold"
    },
    xlabel="scenario",
    # nodes_near_coords,
    # tilt the xtick labels
    xmin=0.75,
    xmax=1.45,
    ymin=0,
    ymax=1.5,
    # shift y label to right
    ylabel_style={yshift="-0.1cm"},
    # do not show numbers on top of bars
    bar_width="15pt",
    legend_cell_align="left",
    legend_pos="north east",
    legend_style={
        # at="{(0.5,1.0)}",
        # anchor="north",
        fill="none",
        draw="none",
        # legend_columns=2,
    },
    "after end axis/.code" = {
        raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(b)};"
    }
},
Plot({color="green", fill="green!60!black!20!white"},Coordinates(
    [
        1
    ],
    [
        maximum_mutual_information_of_fixed_threshold_fixed_T,
    ]
)),
Plot({color="green", fill="green!60!black!20!white",pattern={"north east lines"}},Coordinates(
    [
        1.0
    ],
    [
        maximum_mutual_information_of_fixed_threshold_random_T,
    ]
)),
Plot({color="orange", fill="orange!60!black!20!white"},Coordinates(
    [
        1.2
    ],
    [
        mutual_information_of_dynamical_threshold_fixed_T,
    ]
)),
Plot({color="orange", fill="orange!60!black!20!white",pattern={"north east lines"}},Coordinates(
    [
        1.2
    ],
    [
        mutual_information_of_dynamical_threshold_random_T,
    ]
)),
# legend for each bar
LegendEntry("fixed \$P_{\\rm th}\$, fixed \$T\$"),
LegendEntry("fixed \$P_{\\rm th}\$, random \$T\$"),
LegendEntry("dynamic \$P_{\\rm th}\$, fixed \$T\$"),
LegendEntry("dynamic \$P_{\\rm th}\$, random \$T\$"),

)
push!(tikzpicture, axis)
pgfsave("FIG/mutual_information_of_random_and_fixed_T.pgf", tikzpicture)
pgfsave("FIG/mutual_information_of_random_and_fixed_T.pdf", tikzpicture)
@info "FIG/mutual_information_of_random_and_fixed_T.pdf saved"


# channel capacity
using Optim
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

#### Plot of η_FLD vs time and CCₚ vs time
plot_comparison_1 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="cell contact time \$T\$",
        ylabel="channel capacity \$C\$",
        legend_pos="north west",
        xmode="log",
        legend_cell_align="left",
        legend_style={
            anchor="north west",
            at="(0.02,0.98)",
            draw="none",
            fill="none",
        },
        ymin = 0,
        ymax = 1,
        xmin = 1,
        xmax = 1e6,
    },
    PlotInc({mark="square",color="blue"}, Table(times, 𝐂ₚ)),
    LegendEntry("\$C\\big(\\xi; {\\rm P}(T)\\big)\$"),
)
plot_comparison_2 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="cell contact time \$T\$",
        ylabel="FLD \$\\eta_{\\rm FLD}(T)\$",
        legend_style={
            anchor="north west",
            at="(0.02,0.90)",
            draw="none",
            fill="none",
        },
        xmode="log",
        # ymode="log",
        ymin = 0,
        ymax = 1,
        xmin = 1,
        xmax = 1e6,
        legend_cell_align="left",
        axis_y_line = "right"
    },
    PlotInc({mark="o",color="black"}, Table(times, η_FLD)),
    LegendEntry("\$\\eta_{\\rm FLD}(T)\$"),
)

overlay_plot = TikzPicture(plot_comparison_1, plot_comparison_2)
pgfsave("FIG/multiround_simulation_product_FLD_mutual_information_$(label)_overlay.pgf", overlay_plot)
pgfsave("FIG/multiround_simulation_product_FLD_mutual_information_$(label)_overlay.pdf", overlay_plot)
pgfsave("FIG/multiround_simulation_product_FLD_mutual_information_$(label)_overlay.svg", overlay_plot)
@info "Figure saved to FIG/multiround_simulation_product_FLD_mutual_information_$(label)_overlay.pdf"
# combine two plots into one frame
plot_comparison = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="cell contact time \$T\$",
        legend_pos="north west",
        xmode="log",
        legend_style={
            anchor="north west",
            at="(0.02,0.9)",
            draw="none",
            fill="none",
        },
        legend_cell_align="left",
    },
    plot_comparison_1,
    plot_comparison_2,
)




plot1 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="cell contact time \$T\$",
        ylabel="channel capacity \$C\$",
        legend_pos="north west",
        xmode="log",
        legend_cell_align="left",
        legend_style="{
            draw=none,
            fill=none,
        }"
    },
    PlotInc({mark="o",color="black!10!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[1])),
    PlotInc({mark="o",color="black!20!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[2])),
    PlotInc({mark="o",color="black!30!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[5])),
    PlotInc({mark="o",color="black!50!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[10])),
    PlotInc({mark="o",color="black!70!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[50])),
    PlotInc({mark="o",color="black!100!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[100])),    
    PlotInc({mark="square",color="blue"}, Table(times, 𝐂ₚ)),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 1)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 2)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 5)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 10)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 50)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 100)\$"),
    LegendEntry("\$C\\big(\\xi; {\\rm P}(T)\\big)\$"),
)

pgfsave("FIG/multiround_simulation_product_fpt_analysis_$(label).pgf", plot1)
pgfsave("FIG/multiround_simulation_product_fpt_analysis_$(label).pdf", plot1)
pgfsave("FIG/multiround_simulation_product_fpt_analysis_$(label).svg", plot1)
@info "Figure saved to FIG/multiround_simulation_product_fpt_analysis_$(label).pdf"

plot1 = @pgf Axis(
    {
        width="4in",
        height="4in",
        xlabel="cell contact time \$T\$",
        ylabel="channel capacity \$C\$",
        legend_pos="north west",
        xmode="log",
        legend_cell_align="left",
        legend_style="{
            draw=none,
            fill=none,
        }"
    },
    PlotInc({mark="o",color="black!10!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[1])),
    PlotInc({mark="o",color="black!20!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[2])),
    PlotInc({mark="o",color="black!30!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[5])),
    PlotInc({mark="o",color="black!50!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[10])),
    PlotInc({mark="o",color="black!70!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[50])),
    PlotInc({mark="o",color="black!100!red",no_marks}, Table(times_refined, 𝐂ₚᵀs[100])),    
    PlotInc({mark="square",color="blue"}, Table(times, 𝐂ₚ)),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 1)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 2)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 5)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 10)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 50)\$"),
    LegendEntry("\$C(\\xi; X_{\\rm th} \\mid {\\rm P}_{\\rm th} = 100)\$"),
    LegendEntry("\$C\\big(\\xi; {\\rm P}(T)\\big)\$"),
)
pgfsave("FIG/fig_poster_4.pdf", plot1)