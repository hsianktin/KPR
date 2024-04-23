@info "This script assumes simulation results are already available in DATA/"
@info """
This script analyze the dependence of channel capacity on processing time.
"""
using BSON
using StatsBase

using ProgressMeter
group_label = "kp_1"
τs = collect(0.0:0.2:5.0) ∪ collect(0.0:0.1:1.0) |> sort
profiles = ["PROFILES/profile_$(group_label)_τ_$(τ).jl" for τ in τs]
tasks = [(substrate, id) for substrate in ["k", "q"] for id in 1:99]
labels = [split(profile, "/")[end] |> x -> split(x, ".jl")[1] for profile in profiles]

# collect results
Results_k = []
Results_q = []
for label in labels
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
push!(Results_k, results_k)
push!(Results_q, results_q)
end
include("multiround_simulation_utils.jl");
include("multiround_simulation_product_analysis_utils.jl");
# include(profile);
times = [10.0^i for i in 0:0.2:3]

using DataFrames, CSV
using PGFPlotsX
using LaTeXStrings
using Optim



function handle_results(results_k,results_q,profile)
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
    times_refined = [10.0^i for i in 0:0.05:3]
    𝐂ₚᵀs = [[CCₚᵀ(t=t, 𝐓ₚs=𝐓ₚs, nₜₕ=n)[1] for t in times_refined] for n in 1:100]
    for (n, 𝐂ₚᵀ) in enumerate(𝐂ₚᵀs)
        for (m, 𝐂ₚᵀ) in enumerate(𝐂ₚᵀ)
            push!(df, (k_pars.τ, times_refined[m], 𝐂ₚᵀ, n, "τₚ"))
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
        push!(df, (k_pars.τ, times[n], 𝐂ₚ, 0, "P"))
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
        push!(df, (k_pars.τ, times[n], CCₐs[n][1], 0,"τₐ"))
    end
    for (t, η) in zip(times, η_FLD)
        push!(df, (k_pars.τ, t, η, 0, "η_FLD"))
    end
end
using ProgressMeter
if !isfile("DATA/multiround_simulation_product_CC_summary_$(group_label).csv")
    df = DataFrame(
        τ = Float64[],
        T = Float64[],
        CC = Float64[],
        nₜₕ = Int[],
        type = String[],
    )
    @showprogress for (n, profile) in enumerate(profiles)
        handle_results(Results_k[n], Results_q[n], profile)
    end
    CSV.write("DATA/multiround_simulation_product_CC_summary_$(group_label).csv", df)
else
    df = CSV.read("DATA/multiround_simulation_product_CC_summary_$(group_label).csv", DataFrame)
end



Ts = df.T |> unique
Ts = sort(Ts)[1:1:end]
using PGFPlotsX
Ts = [1000.0]
@showprogress for T ∈ Ts
    # select df with T = T
    df_T = df[df.T .== T, :]
    df_P = df_T[df_T.type .== "P", :]
    df_τₐ = df_T[df_T.type .== "τₐ", :]
    df_τₚ = df_T[df_T.type .== "τₚ", :]
    df_η_FLD = df_T[df_T.type .== "η_FLD", :]

    df_τₚ_1 = df_τₚ[df_τₚ.nₜₕ .== 1, :]
    df_τₚ_5 = df_τₚ[df_τₚ.nₜₕ .== 5, :]
    df_τₚ_10 = df_τₚ[df_τₚ.nₜₕ .== 10, :]
    df_τₚ_20 = df_τₚ[df_τₚ.nₜₕ .== 20, :]
    df_τₚ_30 = df_τₚ[df_τₚ.nₜₕ .== 30, :]
    # if length(df_τₚ.T) == 0 || length(df_τₐ.T) == 0 || length(df_P.T) == 0
    #     continue
    # end
    # plot CC vs τ

    plot1 = @pgf Axis(
        {
            width="5in",
            height="3in",
            xlabel="processing time \$\\tau\$",
            ylabel="channel capacity \$C\$",
            # legend_pos="outer north east",
            legend_style={
                anchor="north east",
                at="(1.0,1.0)",
                draw="none",
                fill="none",
                # two columns
                legend_columns="1",
            },
            legend_cell_align="left",
            # title="T=$T",
            ymax=1.19,
            ymin=-0.1,
            xmin=-0.4,
            xmax=5.4,
        },
        # PlotInc({
        #     mark = "o",
        #     color = "red",
        # },
        # Table(
        #     df_τₐ.τ,
        #     df_τₐ.CC,
        # )
        # ),
        PlotInc({
            mark = "square",
            color = "blue",
        },
        Table(
            df_P.τ,
            df_P.CC,
        )
        ),
        Plot({
            mark = "triangle*",
            color = "black!00!orange",
        },
        Table(
            df_τₚ_1.τ,
            df_τₚ_1.CC,
        )
        ),
        Plot({
            mark = "triangle*",
            color = "black!20!orange",
        },
        Table(
            df_τₚ_5.τ,
            df_τₚ_5.CC,
        )
        ),
        Plot({
            mark = "triangle*",
            color = "black!40!orange",
        },
        Table(
            df_τₚ_10.τ,
            df_τₚ_10.CC,
        )
        ),
        Plot({
            mark = "triangle*",
            color = "black!60!orange",
        },
        Table(
            df_τₚ_20.τ,
            df_τₚ_20.CC,
        )
        ),
        Plot({
            mark = "triangle*",
            color = "black!80!orange",
        },
        Table(
            df_τₚ_30.τ,
            df_τₚ_30.CC,
        )
        ),

        # LegendEntry("\$C(\\xi; X_{\\rm a})\$"),
        LegendEntry("\$C(\\xi; {\\rm P})\$"),
        LegendEntry("\$C(\\xi; X_{\\rm th}  \\mid P_{\\rm th}=1) \$"),
        LegendEntry("\$C(\\xi; X_{\\rm th}  \\mid P_{\\rm th}=5) \$"),
        LegendEntry("\$C(\\xi; X_{\\rm th}  \\mid P_{\\rm th}=10) \$"),
        LegendEntry("\$C(\\xi; X_{\\rm th}  \\mid P_{\\rm th}=20) \$"),
        LegendEntry("\$C(\\xi; X_{\\rm th}  \\mid P_{\\rm th}=30) \$"),
    )
    pgfsave("FIG/multiround_simulation_product_CC_summary_$(group_label)_T_$(T).pdf", plot1)
    # use pdftops to convert pdf to eps by accessing external command
    run(`pdftops -eps FIG/multiround_simulation_product_CC_summary_$(group_label)_T_$(T).pdf FIG/multiround_simulation_product_CC_summary_$(group_label)_T_$(T).eps`)
    pgfsave("FIG/multiround_simulation_product_CC_summary_$(group_label)_T_$(T).pgf", plot1)
    pgfsave("FIG/multiround_simulation_product_CC_summary_$(group_label)_T_$(T).svg", plot1)
    @info "figure saved to FIG/multiround_simulation_product_CC_summary_$(group_label)_T_$(T).pdf"
    # pgfsave("FIG/multiround_simulation_product_CC_FLT_summary_$(group_label)_T_$(T).pdf", overlay_plot)
    # pgfsave("FIG/multiround_simulation_product_CC_FLT_summary_$(group_label)_T_$(T).pgf", overlay_plot)
    # pgfsave("FIG/multiround_simulation_product_CC_FLT_summary_$(group_label)_T_$(T).svg", overlay_plot)
end

@info """
Comparison between channel capacity and FLD
with respect to processing time τ at various T.
"""
begin
plot_comparison_1 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="processing time \$\\tau\$",
        ylabel="channel capacity \$C\$",
        # legend_pos="outer north east",
        legend_style={
            anchor="north west",
            at="(1.12,1.0)",
            draw="none",
            fill="none",
        },
        legend_cell_align="left",
        # title="T=$T",
        ymax=1.0,
        ymin=0.0,
        xmin=0.0,
        xmax=5.0,
    },
 
    # LegendEntry("P"),
)
plot_comparison_2 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="processing time \$\\tau\$",
        ylabel="FLD \$\\eta_{\\rm FLD}(T)\$",
        xmin=0.0,
        xmax=5.0,
        legend_style={
            anchor="north west",
            at="(1.12,0.70)",
            draw="none",
            fill="none",
        },
        # xmode="log",
        legend_cell_align="left",
        axis_y_line = "right"
    },
)


Ts = df.T |> unique
Ts = sort(Ts)[1:1:end]
Ts_selected = Ts[41:16:81]
for T in Ts_selected
    df_T = df[df.T .== T, :]
df_P = df_T[df_T.type .== "P", :]
df_η_FLD = df_T[df_T.type .== "η_FLD", :]

if length(df_P.T) == 0 || length(df_η_FLD.T) == 0 
    continue
end

push!(plot_comparison_1, @pgf Plot({
        mark = "square",
        color = "blue!$(100*(log(T)-log(minimum(Ts_selected)))/(log(maximum(Ts_selected))-log(minimum(Ts_selected))))!black",
    },
    Table(
        df_P.τ,
        df_P.CC,
    )
    ))
push!(plot_comparison_1, LegendEntry("\$C(\\xi; {\\rm P}\\mid T=10^{$(log10(T))})\$"))
push!(plot_comparison_2, @pgf Plot({
        mark = "*",
        color = "red!$(100*(log(T)-log(minimum(Ts_selected)))/(log(maximum(Ts_selected))-log(minimum(Ts_selected))))!black",
    },
    Table(
        df_η_FLD.τ,
        df_η_FLD.CC,
    )
    ))
push!(plot_comparison_2, LegendEntry("\$\\eta_{\\rm FLD}(T=10^{$(log10(T))})\$"))
end
plot_comparison_1
plot_comparison_2


overlay_plot = TikzPicture(plot_comparison_1, plot_comparison_2)
pgfsave("FIG/multiround_simulation_product_CC_summary_$(group_label)_comparison.pdf", overlay_plot)
pgfsave("FIG/multiround_simulation_product_CC_summary_$(group_label)_comparison.pgf", overlay_plot)
pgfsave("FIG/multiround_simulation_product_CC_summary_$(group_label)_comparison.svg", overlay_plot)
@info "comparison figure saved to FIG/multiround_simulation_product_CC_summary_$(group_label)_comparison.pdf"

end

@info """
Comparison between channel capacity of products 
and channel capacity of first activation time
with respect to processing time τ at various T.
"""
begin
plot_comparison_1 = @pgf Axis(
    {
        width="5in",
        height="3in",
        xlabel="processing time \$\\tau\$",
        ylabel="channel capacity \$C\$",
        # legend_pos="outer north east",
        legend_style={
            anchor="north west",
            at="(1.02,1.0)",
            draw="none",
            fill="none",
        },
        legend_cell_align="left",
        # title="T=$T",
        ymax=1.0,
        ymin=0.0,
        xmin=0.0,
        xmax=5.0,
    },
 
    # LegendEntry("P"),
)

Ts = df.T |> unique
Ts = sort(Ts)[1:1:end]
Ts_selected = Ts[41:12:81]
for T in Ts_selected
    df_T = df[df.T .== T, :]
df_P = df_T[df_T.type .== "P", :]
df_a = df_T[df_T.type .== "τₐ", :]

if length(df_P.T) == 0 || length(df_a.T) == 0 
    continue
end

push!(plot_comparison_1, @pgf Plot({
        mark = "square",
        color = "blue!$(100*(log(T)-log(minimum(Ts_selected)))/(log(maximum(Ts_selected))-log(minimum(Ts_selected))))!black",
    },
    Table(
        df_P.τ,
        df_P.CC,
    )
    ))
push!(plot_comparison_1, LegendEntry("\$C(\\xi;~ {\\rm P} \\mid T=10^{$(log10(T))})\$"))
end

for T in Ts_selected
    df_T = df[df.T .== T, :]
    df_P = df_T[df_T.type .== "P", :]
    df_a = df_T[df_T.type .== "τₐ", :]

    if length(df_P.T) == 0 || length(df_a.T) == 0 
        continue
    end
    push!(plot_comparison_1, @pgf Plot({
            mark = "*",
            color = "red!$(100*(log(T)-log(minimum(Ts_selected)))/(log(maximum(Ts_selected))-log(minimum(Ts_selected))))!black",
        },
        Table(
            df_a.τ,
            df_a.CC,
        )
        ))
    push!(plot_comparison_1, LegendEntry("\$C(\\xi; X_{\\rm a} \\mid T=10^{$(log10(T))})\$"))
end

plot_comparison_1
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_comparison.pdf", plot_comparison_1)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_comparison.pgf", plot_comparison_1)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_comparison.svg", plot_comparison_1)
@info "comparison figure saved to FIG/multiround_simulation_product_CC_P_A_$(group_label)_comparison.pdf"

end

begin
    plot_comparison_1 = @pgf Axis(
        {
            width="4in",
            height="4in",
            xlabel="processing time \$\\tau\$",
            ylabel="channel capacity \$C\$",
            # legend_pos="outer north east",
            legend_style={
                # anchor at top center
                anchor="north",
                at="(0.5,1.0)",
                draw="none",
                fill="none",
                # two columns
                legend_columns="2",
            },
            legend_cell_align="left",
            # title="T=$T",
            ymax=1.5,
            ymin=0.0,
            xmin=0.0,
            xmax=5.0,
        },
     
        # LegendEntry("P"),
    )
    
    Ts_selected = Ts[41:12:81]
    for T in Ts_selected
        df_T = df[df.T .== T, :]
    df_P = df_T[df_T.type .== "P", :]
    df_a = df_T[df_T.type .== "τₐ", :]
    
    if length(df_P.T) == 0 || length(df_a.T) == 0 
        continue
    end
    
    push!(plot_comparison_1, @pgf Plot({
            mark = "square",
            color = "blue!$(100*(log(T)-log(minimum(Ts_selected)))/(log(maximum(Ts_selected))-log(minimum(Ts_selected))))!black",
        },
        Table(
            df_P.τ,
            df_P.CC,
        )
        ))
    push!(plot_comparison_1, LegendEntry("\$C(\\xi;~ {\\rm P} \\mid T=10^{$(log10(T))})\$"))
    df_T = df[df.T .== T, :]
    df_P = df_T[df_T.type .== "P", :]
    df_a = df_T[df_T.type .== "τₐ", :]

    if length(df_P.T) == 0 || length(df_a.T) == 0 
        continue
    end
    push!(plot_comparison_1, @pgf Plot({
            mark = "*",
            color = "red!$(100*(log(T)-log(minimum(Ts_selected)))/(log(maximum(Ts_selected))-log(minimum(Ts_selected))))!black",
        },
        Table(
            df_a.τ,
            df_a.CC,
        )
        ))
    push!(plot_comparison_1, LegendEntry("\$C(\\xi; X_{\\rm a} \\mid T=10^{$(log10(T))})\$"))

    end
    
    # for T in Ts_selected
    # end
    
    plot_comparison_1
    pgfsave("FIG/fig_poster_3.pdf", plot_comparison_1)
    pgfsave("FIG/fig_poster_3.pgf", plot_comparison_1)
    pgfsave("FIG/fig_poster_3.svg", plot_comparison_1)
    @info "comparison figure saved to FIG/fig_poster_3.pdf"
    
    end
    

Df_P = df[df.type .== "P", :]
Df_a = df[df.type .== "τₐ", :]
# for each Df_P, Df_a; and for each T ∈ Df_P.T, 
# find τ that maximizes CC_P and CC_a, respectively
# and plot CC_P(τ) and CC_a(τ) vs T
Ts_raw = Df_P.T |> unique
Ts = Ts_raw[4:1:20]

function get_τ_max(df, T)
    df_T = df[df.T .== T, :]
    τ_max_P = median(df_T[df_T.CC .≈ maximum(df_T.CC), :].τ)
    return τ_max_P
end

function get_T_max(df, τ)
    df_τ = df[df.τ .== τ, :]
    T_max_P = median(df_τ[df_τ.CC .== maximum(df_τ.CC), :].T)
    T_max_P = (df_τ[df_τ.CC .≈ maximum(df_τ.CC), :].T[end])
    return T_max_P
end

τs = df.τ |> unique

τ_max_Ps = [get_τ_max(Df_P, T) for T in Ts]
T_max_Ps = [get_T_max(Df_P, τ) for τ in τs]
τ_max_As = [get_τ_max(Df_a, T) for T in Ts]
T_max_As = [get_T_max(Df_a, τ) for τ in τs]


plt= @pgf Axis({
    width="2.5in",
    height="2.5in",
    xlabel="cell contact time \$T\$",
    ylabel="optimal processing time \$\\tau^*\$",
    legend_pos="north west",
    legend_cell_align="left",
    legend_style={
        anchor="north west",
        draw="none",
        fill="none",
    },
    xmode="log",
},
PlotInc({
    mark = "o",
    color = "blue",
}, Table(
    Ts,
    τ_max_Ps,
)),
PlotInc({
    mark = "square",
    color = "red",
}, Table(
    Ts,
    τ_max_As,
)),
LegendEntry("\$X_{\\rm th}^*\$"),
LegendEntry("\$\\tau_{\\rm a}^*\$"),
)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_τ_max.pdf", plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_τ_max.pgf", plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_τ_max.svg", plt)
@info "τ_max figure saved to FIG/multiround_simulation_product_CC_P_A_$(group_label)_τ_max.pdf"

# plot T_max 
plt= @pgf Axis({
    width="2.5in",
    height="2.5in",
    xlabel="processing time \$\\tau\$",
    ylabel="optimal cell contact time \$T^*\$",
    legend_pos="north west",
    legend_cell_align="left",
    legend_style={
        anchor="north west",
        draw="none",
        fill="none",
    },
    ymode="log",
},
PlotInc({
    mark = "o",
    color = "blue",
}, Table(
    τs,
    T_max_Ps,
)),
PlotInc({
    mark = "square",
    color = "red",
}, Table(
    τs,
    T_max_As,
)),
LegendEntry("\$T^*_{\\rm P}\$"),
LegendEntry("\$T^*_{\\rm a}\$"),
)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_T_max.pdf", plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_T_max.pgf", plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_T_max.svg", plt)

using PGFPlotsX

plt = @pgf begin
    GroupPlot({
        group_style = {
            group_size = "2 by 1",
            xlabels_at = "edge bottom",
            ylabels_at = "edge left",
            horizontal_sep = "1.5cm",
            vertical_sep = "1.5cm"
        },
        width = "2.5in",
        height = "2.5in"
    }, 

    {
        xlabel = "cell contact time \$T\$",
        ylabel = "optimal processing time \$\\tau^*\$",
        # legend_pos = "north east",
        legend_cell_align = "left",
        legend_style = {
            anchor = "north west",
            at = "(0.05,0.85)",
            draw = "none",
            fill = "none",
        },
        xmode = "log",
        "enlarge y limits" = 0.1,
        "after end axis/.code" = {
            raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(a)};"
        }
    },
    PlotInc({
        mark = "o",
        color = "blue",
    }, Table(
        Ts,
        τ_max_Ps,
    )),
    PlotInc({
        mark = "square",
        color = "red",
    }, Table(
        Ts,
        τ_max_As,
    )),
    LegendEntry("\$X_{\\rm th}^*\$"),
    LegendEntry("\$\\tau_{\\rm a}^*\$"),

    {
        xlabel = "processing time \$\\tau\$",
        ylabel = "optimal cell contact time \$T^*\$",
        # legend_pos = "south east",
        legend_cell_align = "left",
        legend_style = {
            anchor = "north west",
            at = "(0.05,0.85)",
            draw = "none",
            fill = "none",
        },
        xtick = "{0,1,2,3,4,5}",
        ymin=1,
        ymax=1e7,
        ymode = "log",
        "after end axis/.code" = {
            raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(b)};"
        }
    },
    PlotInc({
        mark = "o",
        color = "blue",
    }, Table(
        τs,
        T_max_Ps,
    )),
    PlotInc({
        mark = "square",
        color = "red",
    }, Table(
        τs,
        T_max_As,
    )),
    LegendEntry("\$T^*_{\\rm P}\$"),
    LegendEntry("\$T^*_{\\rm a}\$")
    )
end

pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.pdf",plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.pgf",plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.svg",plt)
@info "Grouped plots with TikZ labels saved to FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.pdf"

function N(T, pars::Pars)
    @unpack k₁, k₋₁, τ = pars  # Unpack parameters from the struct
    term1 = 1 / k₁
    term2 = 1 / k₋₁ * (1 - (τ * k₋₁ + 1) * exp(-k₋₁ * τ)) / (1 - exp(-k₋₁ * τ))
    N_approx = T / (term1 + term2)
    return N_approx
end

function Tₙ(N, pars::Pars)
    @unpack k₁, k₋₁, τ = pars  # Unpack parameters from the struct
    term1 = 1 / k₁
    if τ == 0
        term2 = 0
    else
        term2 = 1 / k₋₁ * (1 - (τ * k₋₁ + 1) * exp(-k₋₁ * τ)) / (1 - exp(-k₋₁ * τ))
    end
    Tₙ_approx = N * (term1 + term2)
    return Tₙ_approx
end

function N₊(k_pars::Pars, q_pars::Pars)
    @unpack k₁, k₋₁, τ = k_pars  # Unpack parameters from the struct
    q₋₁ = q_pars.k₋₁
    if τ == 0
        return 1
    else
        return (q₋₁ - k₋₁) * τ * exp(k₋₁ * τ) / (1 - exp(-(q₋₁ - k₋₁) * τ))
    end
end

function T₊(k_pars::Pars, q_pars::Pars)
    return Tₙ(N₊(k_pars, q_pars), k_pars)
end


function τᵖₘₐₓ(k_pars, q_pars)
    @unpack k₁, k₋₁, τ = k_pars  # Unpack parameters from the struct
    q₋₁ = q_pars.k₋₁
    q₁ = q_pars.k₁
    K = k₁ / (k₁ + k₋₁)
    Q = q₁ / (q₁ + q₋₁)
    return 2 / (q₋₁ - k₋₁) * log(sqrt(Q / K) * (q₋₁ / k₋₁))
end

function τᵃₘₐₓ(T,k_pars, q_pars)
    @unpack k₁, k₋₁, τ = k_pars  # Unpack parameters from the struct
    q₋₁ = q_pars.k₋₁
    q₁ = q_pars.k₁
    Nₜ = T / k₁^(-1)
    return (log(Nₜ) + log(q₋₁/k₋₁))/q₋₁
end
K_pars = []
Q_pars = []
for profile in profiles
    include(profile)
    push!(K_pars, deepcopy(k_pars))
    push!(Q_pars, deepcopy(q_pars))
end
k_par = K_pars[1]
q_par = Q_pars[1]
τᵃₘₐₓs = [τᵃₘₐₓ(T, k_par, q_par) for T in Ts]
τᵖₘₐₓs = [τᵖₘₐₓ(k_par, q_par) for T in Ts]
T₊s = [T₊(K_pars[n], Q_pars[n]) for n in eachindex(profiles)]
plt = @pgf begin
    GroupPlot({
        group_style = {
            group_size = "2 by 1",
            xlabels_at = "edge bottom",
            ylabels_at = "edge left",
            horizontal_sep = "1.5cm",
            vertical_sep = "1.5cm"
        },
        width = "2.5in",
        height = "2.5in"
    }, 


    {
        xlabel = "processing time \$\\tau\$",
        ylabel = "optimal cell contact time \$T^{\\rm o}\$",
        # legend_pos = "south east",
        legend_cell_align = "left",
        legend_style = {
            anchor = "north west",
            at = "(0.15,0.75)",
            draw = "none",
            fill = "none",
        },
        xtick = "{0,1,2,3,4,5}",
        ymin=1,
        ymax=1e7,
        ymode = "log",
        "after end axis/.code" = {
            raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(a)};"
        }
    },
    PlotInc({
        mark = "o",
        color = "blue",
    }, Table(
        τs,
        T_max_Ps,
    )),
    PlotInc({
        mark = "square",
        color = "red",
        only_marks = true,
        mark_size="1.5pt",
    }, Table(
        τs,
        T_max_As,
    )),
    PlotInc({
        no_marks,
        color = "red",
        "dashed",
        thick,
    }, Table(
        τs,
        T₊s,
    )),
    LegendEntry("\$T_{\\rm P}^{\\rm o}\$"),
    LegendEntry("\$T_{\\rm a}^{\\rm o}\$"),
    LegendEntry("\$\\widehat T_{\\rm a}^{\\rm o}\$"),
    {
        xlabel = "cell contact time \$T\$",
        ylabel = "optimal processing time \$\\tau^{\\rm o}\$",
        # legend_pos = "north east",
        legend_cell_align = "left",
        legend_style = {
            anchor = "north west",
            at = "(0.15,0.75)",
            draw = "none",
            fill = "none",
        },
        xmode = "log",
        "enlarge y limits" = 0.1,
        "after end axis/.code" = {
            raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(b)};"
        }
    },
    PlotInc({
        mark = "o",
        color = "blue",
        only_marks,
    }, Table(
        Ts,
        τ_max_Ps,
    )),
    PlotInc(
        {
            "no marks",
            color = "blue",
            thick,
            "dashed",
        },
        Table(
            Ts,
            τᵖₘₐₓs,
        )
    ),
    PlotInc({
        mark = "square",
        color = "red",
        mark_size="1.5pt",
        only_marks,
    }, Table(
        Ts,
        τ_max_As,
    )),
    # PlotInc(
    #     {
    #         "no marks",
    #         color = "red",
    #         thick,
    #         "dashed",
    #     },
    #     Table(
    #         Ts,
    #         τᵃₘₐₓs,
    #     )
    # ),
    LegendEntry("\$\\tau^{\\rm o}_{\\rm P}\$"),
    LegendEntry("\$\\hat \\tau^{\\rm o}_{\\rm P}\$"),
    LegendEntry("\$\\tau_{\\rm a}^{\\rm o}\$"),
    # LegendEntry("\$\\hat \\tau_{\\rm a}^*\$"),
    )
end

pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.pdf",plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.pgf",plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.svg",plt)
pgfsave("FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.tex",plt)

@info "Grouped plots with TikZ labels saved to FIG/multiround_simulation_product_CC_P_A_$(group_label)_grouped_plots.pdf"

# estimate mutual information from δ
using SpecialFunctions
function σₘₐₓ(T,k_pars, q_pars)
    @unpack k₁, k₋₁, τ = k_pars  # Unpack parameters from the struct
    @unpack kₚ = k_pars
    q₋₁ = q_pars.k₋₁
    q₁ = q_pars.k₁
    K = k₁ / (k₁ + k₋₁)
    Q = q₁ / (q₁ + q₋₁)
    return √(kₚ * T) * (√K * exp(- k₋₁*τ/2) - √Q * exp(- q₋₁*τ/2))
end

function εₘᵢₙ(σ)
    return 0.5 * (1 - erf(σ / √2))
end

function S_Bernoulli(p)
    # entropy of Bernoulli distribution with parameter p
    if p == 0 || p == 1
        return 0
    else
        return -p * log(p) - (1 - p) * log(1 - p)
    end
end

function Ĉ(T, k_pars, q_pars)
    return 1 - S_Bernoulli(εₘᵢₙ(σₘₐₓ(T, k_pars, q_pars)))/log(2)
end


Ĉs = [Ĉ(T, K_pars[16], Q_pars[16]) for T in Ts]
df_selected = df[df.τ .== 3.0, :]
df_p_selected = df_selected[df_selected.type .== "P", :]

Ĉₜs = [Ĉ(1000.0, k_par, q_par) for (k_par, q_par) in zip(K_pars, Q_pars)]
τₜs = [k_par.τ for k_par in K_pars]

df_τ_selected = df[df.T .== 1000.0, :]
df_p_τ_selected = df_τ_selected[df_τ_selected.type .== "P", :]


plt = @pgf begin
    Axis({
        width = "2.5in",
        height = "2.5in",
        xlabel = "cell contact time \$T\$",
        ylabel = "estimated channel capacity \$\\hat C\$",
        legend_pos = "north west",
        legend_cell_align = "left",
        legend_style = {
            anchor = "north west",
            draw = "none",
            fill = "none",
        },
        xmode = "log",
        "enlarge y limits" = 0.1,

    },
    PlotInc({
        mark = "o",
        color = "blue",
        mark_size = "1.5pt",
        no_marks,
        thick,
        dashed,
        mark_repeat = 5,
    }, Table(
        Ts,
        Ĉs,
    )),
    LegendEntry("\$\\hat C\$"),
    # plot actual C
    PlotInc({
        mark = "square",
        only_marks,
        mark_size="1.5pt",
        color = "blue",
        # "dashed",
    }, Table(
        df_p_selected.T,
        df_p_selected.CC,
    )),
    LegendEntry("\$C\$"),
    )

end

plt = @pgf begin
    Axis({
        width = "2.5in",
        height = "2.5in",
        xlabel = "processing time \$\\tau\$",
        ylabel = "estimated channel capacity \$\\hat C\$",
        legend_pos = "north east",
        legend_cell_align = "left",
        legend_style = {
            anchor = "north east",
            draw = "none",
            fill = "none",
        },
        # xmode = "log",
        "enlarge y limits" = 0.1,
    },
    PlotInc({
        mark = "o",
        color = "blue",
        mark_size = "1.5pt",
        only_marks,
    },Table(
        df_p_τ_selected.τ,
        df_p_τ_selected.CC,
    )),
    LegendEntry("\$\\hat C\$"),
    # plot actual C
    PlotInc({
        # mark = "square",
        no_marks,
        color = "blue",
        "dashed",
    },  Table(
        τₜs,
        Ĉₜs,
    )),
    LegendEntry("\$C\$"),
    )
end

plt = @pgf begin
    GroupPlot({
        group_style = {
            group_size = "2 by 1",
            xlabels_at = "edge bottom",
            ylabels_at = "edge left",
            horizontal_sep = "0.2cm",
            vertical_sep = "1.5cm"
        },
        width = "2.5in",
        height = "2.5in",
        "after end axis/.code" = {
            raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(a)};"
        }
    }, 
    {
        width = "2.5in",
        height = "2.5in",
        xlabel = "cell contact time \$T\$",
        ylabel = "channel capacity \$C\$",
        legend_pos = "south east",
        legend_cell_align = "left",
        legend_style = {
            anchor = "south east",
            draw = "none",
            fill = "none",
        },
        xmode = "log",
        ymin=-0.2,
        ymax=1.2,
        # "enlarge y limits" = 0.1,
    },
    PlotInc({
        mark = "o",
        color = "blue",
        mark_size = "1.5pt",
        no_marks,
        thick,
        dashed,
        mark_repeat = 5,
    }, Table(
        Ts,
        Ĉs,
    )),
    LegendEntry("\$\\hat C(\\xi; {\\rm P}(T))\$"),
    # plot actual C
    PlotInc({
        mark = "square",
        only_marks,
        mark_size="1.5pt",
        color = "blue",
        # "dashed",
    }, Table(
        df_p_selected.T,
        df_p_selected.CC,
    )),
    LegendEntry("\$C(\\xi; {\\rm P}(T))\$"),
    
    {
        width = "2.5in",
        height = "2.5in",
        xlabel = "processing time \$\\tau\$",
        # ylabel = "estimated channel capacity \$\\hat C\$",
        legend_pos = "north east",
        legend_cell_align = "left",
        legend_style = {
            anchor = "north east",
            draw = "none",
            fill = "none",
        },
        ymin=-0.2,
        ymax=1.2,
        yticklabels = {},
        "after end axis/.code" = {
            raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont] at (rel axis cs:.025,.975) {(b)};"
        }
        # xmode = "log",
    },
    # plot actual C
    PlotInc({
        # mark = "square",
        no_marks,
        color = "blue",
        "dashed",
        thick,
    },  Table(
        τₜs,
        Ĉₜs,
    )),
    PlotInc({
        mark = "square",
        color = "blue",
        mark_size = "1.5pt",
        only_marks,
    },Table(
        df_p_τ_selected.τ,
        df_p_τ_selected.CC,
    )),

    # LegendEntry("\$\\hat C\$"),
    # LegendEntry("\$C\$"),

    )
end
pgfsave("FIG/multiround_simulation_product_CC_$(group_label)_estimates.pgf",plt)
pgfsave("FIG/multiround_simulation_product_CC_$(group_label)_estimates.pdf",plt)
pgfsave("FIG/multiround_simulation_product_CC_$(group_label)_estimates.svg",plt)
# call pdftops to convert pdf to eps
run(`pdftops -eps FIG/multiround_simulation_product_CC_$(group_label)_estimates.pdf FIG/multiround_simulation_product_CC_$(group_label)_estimates.eps`)
@info "Grouped plots with TikZ labels saved to FIG/multiround_simulation_product_CC_$(group_label)_estimates.pdf"