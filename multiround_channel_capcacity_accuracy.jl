### DNA accuracy vs. processing time τ
function DNA_inaccuracy(k1,q1,k_1,q_1,τ)
    return q1 * exp(τ * (k_1-q_1))/(k1 + q1 * exp(τ * (k_1-q_1)))
end

τs = [exp(i) for i in 0:0.01:2.5]
DNA_inaccuracies = [DNA_inaccuracy(0.1,0.1,1,1/α,τ) for τ in τs]
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{amssymb}")

tikzpic = @pgf TikzPicture({baseline})
plot1 = @pgf Axis(
    {width="3.4in", 
    height="3.4in", 
    xlabel=raw"processing time $\tau$", 
    ylabel=raw"error probability $\mathbb{P}(t_{\rm p} \geq t_{\rm p'})$",
    legend_pos="north east",
    legend_style={font="\\large", fill="none", draw="none"},
    legend_cell_align="left",
    # large font size and tick label style
    tick_label_style={font="\\large"},
    x_label_style={font="\\Large"},
    y_label_style={font="\\Large"},
    # xmin=-5, xmax=5,
    # xmode="log",
    ymode="log",
    # enlarge_y_limits=0.1,
    #ymin=0.5,ymax=1,
    "after end axis/.code" = {
        raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont, fill=white] at (rel axis cs:.02,.98) {(a)};"
    },
    },
    # Plot({color="orange!50!black",mark="square", mark_repeat=5}, Table((τs), 1 .- Accuracy_numerical_maximum)),
    Plot({color="blue",ultra_thick,}, Table((τs), DNA_inaccuracies)),
    # LegendEntry("Channel Capacity"),
    # LegendEntry(L"1-\mathcal{A}^*(\tau)"),
    # LegendEntry(L"1-\mathcal{A}_{\mathrm{asymp}}(\tau)"),
)
push!(tikzpic, plot1)
pgfsave("FIG/multiround_DNA_accuracy_vs_tau.pgf", tikzpic)
pgfsave("FIG/multiround_DNA_accuracy_vs_tau.pdf", tikzpic)
pgfsave("FIG/multiround_DNA_accuracy_vs_tau.svg", tikzpic)

### TCR signaling accuracy & channel vs. processing time τ

using Optim

# Function to model conditional probability P(X|ξ)
# single binding

function ℙXₗξN(X, ξ; τ = 1, α = 0.5, N = 20)
    p_x1_given_xi1 = 1 - (1 - exp(- τ))^N
    p_x0_given_xi0 = (1 - exp(- τ/α))^N
    if ξ == 1
        if X == 1
            return p_x1_given_xi1
        else
            return 1 - p_x1_given_xi1
        end
    else
        if X == 0
            return p_x0_given_xi0
        else
            return 1 - p_x0_given_xi0
        end
    end
end
∑ = sum  # Alias for sum function
"""
    Accuracy(ℙXₗξN)
    Compute the accuracy of a binary symmetric channel given the conditional probability P(X|ξ) for i = 0,1.
    Accuarcy = ℙXₗξ(X=1,ξ=1) + ℙXₗξ(X=0,ξ=0)
"""
function accuracy(ℙXₗξ; τ = 1, α = 0.5, N = 20, ℙξ=0.5)
    ℙX1ξ1 = ℙXₗξ(1, 1; τ = τ, α = α, N = N) * ℙξ
    ℙX0ξ0 = ℙXₗξ(0, 0; τ = τ, α = α, N = N) * (1 - ℙξ)
    return ℙX1ξ1 + ℙX0ξ0
end



"""
    channel_capacity(ℙXₗξN)
    Compute the channel capacity of a binary symmetric channel given the conditional probability P(X|ξ) for i = 0,1.
    ℙXₗξ is a function that takes two arguments X and ξ and returns the conditional probability P(X|ξ).
"""
function channel_capacity(ℙXₗξ)
    # Objective function to maximize mutual information I(ξ; X)
    function mutual_information(p)
        ℙξ₀ = p[1]
        ℙξ₁ = 1 - ℙξ₀  # Since ℙξ₀ + ℙξ₁ = 1
        function ℙξ(ξ)
            if ξ == 0
                return ℙξ₀
            else
                return ℙξ₁
            end
        end
        # Calculate P(X=0) and P(X=1) using total probability
        ℙX0 = ℙXₗξ(0, 0) * ℙξ(0) + ℙXₗξ(0, 1) * ℙξ(1)
        ℙX1 = ℙXₗξ(1, 0) * ℙξ(0) + ℙXₗξ(1, 1) * ℙξ(1)
        function ℙX(X)
            if X == 0
                return ℙX0
            else
                return ℙX1
            end
        end
        function log₂(x)
            if x == 0
                return 0
            else
                return log2(x)
            end
        end
        function ℙXξ(X,ξ)
            return ℙXₗξ(X,ξ) * ℙξ(ξ)
        end
        I = ∑(ℙXξ(X,ξ) * log₂(ℙXξ(X,ξ) / (ℙX(X) * ℙξ(ξ))) 
                    for X ∈ [0,1], ξ ∈ [0,1])

        return -I  # Negate because Optim.jl minimizes
    end

    # Initial guess for p₀
    initial_guess = [0.5]
    lower_bounds = [0.0]  # Lower bound for p₀
    upper_bounds = [1.0]  # Upper bound for p₀
    
    # Perform the optimization
    opt_result = optimize(mutual_information, lower_bounds, upper_bounds ,initial_guess, Fminbox(LBFGS()); inplace = false)

    # Extract optimal input distribution and channel capacity
    optimal_p₀ = opt_result.minimizer[1]
    channel_capacity = -opt_result.minimum  # Negate to get the maximum

    return (channel_capacity, optimal_p₀)
end

T = [exp(i) for i in -5:0.05:5]
Ns = collect(1:100)
α = 0.5

CC = [channel_capacity((X,ξ) -> ℙXₗξN(X,ξ; τ = t, α = α, N=N))[1] for t in T, N in Ns]
# replace NaN by 0
CC = [if isnan(c) 0 else c end for c in CC]
sensitivity(τ,α,N) = 1 - (1 - exp(- τ))^N
specificity(τ,α,N) = (1 - exp(- τ/α))^(N)

sensitivities = [sensitivity(t, α, N) for t in T, N in Ns]
specificities = [specificity(t, α, N) for t in T, N in Ns]
accuracies = [accuracy(ℙXₗξN;τ=t, N=N) for t in T, N in Ns]

using PGFPlotsX
using LaTeXStrings
# pick N = 1, 5, 10, 50, 100
# plot ROC curve
# use the incremental colorscheme:
# \pgfplotscreateplotcyclelist{mark list 5}{%
# red!80!black,every mark/.append style={color=red!80!black},mark=o,mark size=1.8\\%
# color=brown!80!black,every mark/.append style={color=brown!80!black},mark=triangle*\\%
# color=black,every mark/.append style={color=black},mark=square,mark size=1.8\\%
# color=green!50!black,every mark/.append style={color=green!50!black},mark=diamond\\%
# blue,every mark/.append style={color=blue},mark=pentagon*\\%
# purple, every mark/.append style={color=purple},mark={10-pointed star}\\%
# }
# legend style = {
#     font = \footnotesize,
#     fill = none,
#     draw = none,
# },
# legend cell align={left},
plot = @pgf Axis({width="3.4in", 
            height="3.4in", 
            xlabel="1 - Specificity", 
            ylabel="Sensitivity",
            legend_pos="south east",
            legend_style={font="\\footnotesize", fill="none", draw="none"},
            legend_cell_align="left",
            xmin=0, xmax=1,
            ymin=0,ymax=1,},
    PlotInc({color="red!80!black", mark="none"}, Table(1 .- specificities[:,1], sensitivities[:,1])),
    PlotInc({color="brown!80!black", mark="none"}, Table(1 .- specificities[:,5], sensitivities[:,5])),
    PlotInc({color="black", mark="none"}, Table(1 .- specificities[:,10], sensitivities[:,10])),
    PlotInc({color="green!50!black", mark="none"}, Table(1 .- specificities[:,50], sensitivities[:,50])),
    PlotInc({color="blue", mark="none"}, Table(1 .- specificities[:,100], sensitivities[:,100])),
            LegendEntry("N = 1"),
    LegendEntry("N = 5"),
    LegendEntry("N = 10"),
    LegendEntry("N = 50"),
    LegendEntry("N = 100"),
)
pgfsave("FIG/multiround_ROC_curve.pgf", plot)
pgfsave("FIG/multiround_ROC_curve.pdf", plot)
pgfsave("FIG/multiround_ROC_curve.svg", plot)


# then do a heatmap plot of channel capacity, with x-axis = τ, y-axis = N
# use  /pgfplots/colormap/hsv2/.cd as the colormap
# heatmap in pgfplots is implemented in the form of a surf plot


T = [exp(i) for i in -5:0.05:5]
Ns = exp.([i for i in 0:0.1:10])
α = 0.5

CC = [channel_capacity((X,ξ) -> ℙXₗξN(X,ξ; τ = t, α = α, N=N))[1] for t in T, N in Ns]
Accuracy = [accuracy(ℙXₗξN; τ=t, N=N) for t in T, N in Ns]
# replace NaN by 0
CC = [if isnan(c) 0 else c end for c in CC]

using DataFrames
using CSV 
df = DataFrame(
    logT = [],
    logN = [],
    CC = [],
    Accuracy = []
)
for i in 1:length(T)
    for j in 1:length(Ns)
        push!(df, (log(T[i]), log(Ns[j]), CC[i,j], Accuracy[i,j]))
    end
end
CSV.write("DATA/heatmap_data.csv", df)

τs = [exp(i) for i in 0:0.01:2.5]
# assume k₋₁ = 1, q₋₁= 1/ α
N⁺(α,τ) = (1/α - 1) * τ * exp(τ)/ (1 - exp(-(1/α-1)τ))
accuracy(τ,α) = accuracy(ℙXₗξN; τ=τ, α=α, N=N⁺(α,τ))
α = 0.5
Accuracy⁺ = [accuracy(τ,α) for τ in τs]
N⁺s = [N⁺(α,τ) for τ in τs]
accuracy_asymptotic(τ,α) = 1 - (1/α -1) * τ * exp(-(1/α -1)*τ)/2
Accuracy_asymptotic = [accuracy_asymptotic(τ,α) for τ in τs]
Ns = exp.([i*log(10) for i in 0:0.1:40])
Accuracies = [accuracy(ℙXₗξN; τ=τ, N=N, α=α) for τ in τs, N in Ns]
# replace NaN by 0
Accuracies = [if isnan(c) 0 else c end for c in Accuracies]
Accuracy_numerical_maximum = [maximum([accuracy(ℙXₗξN; τ=τ, N=N, α=α) for N in Ns]) for τ in τs]

tikzpic = @pgf TikzPicture({baseline})
plot2 = @pgf Axis(
    {width="3.4in", 
    height="3.4in", 
    xlabel=raw"processing time $\tau$", 
    ylabel=raw"inaccuracy $1 - \mathcal{A}$",
    legend_pos="north east",
    legend_style={font="\\large", fill="none", draw="none"},
    legend_cell_align="left",
    # large font size and tick label style
    tick_label_style={font="\\large"},
    x_label_style={font="\\Large"},
    y_label_style={font="\\Large"},
    # xmin=-5, xmax=5,
    # xmode="log",
    ymode="log",
    # enlarge_x_limits=false,
    #ymin=0.5,ymax=1,
    "after end axis/.code" = {
        raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont, fill=white] at (rel axis cs:.02,.98) {(c)};"
    },

    },
    Plot({color="red",mark="square", mark_repeat=5,only_marks}, Table((τs), 1 .- Accuracy_numerical_maximum)),
    Plot({color="red",dashed,ultra_thick}, Table((τs), 1 .-Accuracy_asymptotic)),
    # LegendEntry("Channel Capacity"),
    LegendEntry(L"1-\mathcal{A}^*(\tau)"),
    LegendEntry(L"1-\mathcal{A}^*_{\mathrm{asymp}}(\tau)"),
)
push!(tikzpic, plot2)

pgfsave("FIG/multiround_accuracy_vs_tau.pgf", tikzpic)
pgfsave("FIG/multiround_accuracy_vs_tau.pdf", tikzpic)
pgfsave("FIG/multiround_accuracy_vs_tau.svg", tikzpic)

### total duration vs. processing time τ


function expected_t(k1, k_1, q1, q_1, τ)
    numerator = (k1/k_1) * (1 - exp(-k_1 * τ)) + (q1/q_1) * (1 - exp(-q_1 * τ)) + 1
    denominator = k1 * exp(-k_1 * τ) + q1 * exp(-q_1 * τ)

    exact = numerator / denominator

    # Approximation check
    # if k1 * exp(-k_1 * τ) >> q1 * exp(-q_1 * τ)
        approximate = exp(k_1 * τ) / k1 * (k1/k_1 + q1/q_1 + 1)
    # else
        # approximate = NaN # Not a valid approximation
    # end

    return exact, approximate
end

𝔼ts = [expected_t(0.1,1,0.1,1/α,τ)[1] for τ in τs]
𝔼̂ts = [expected_t(0.1,1,0.1,1/α,τ)[2] for τ in τs]


function optimal_T(α,τ)
    return N⁺(α,τ) / 0.1 # k1 = q1 = 0.1
end

Tꜛs = [optimal_T(α,τ) for τ in τs]

tikzpic = @pgf TikzPicture({baseline})
# plot3 = @pgf Axis(
#     {width="3.4in", 
#     height="3.4in", 
#     xlabel=raw"processing time $\tau$", 
#     ylabel=raw"total duration",
#     legend_pos="north east",
#     legend_style={font="\\large", fill="none", draw="none",at="(0.05,0.9)", anchor="north west"},
#     legend_cell_align="left",
#     # large font size and tick label style
#     tick_label_style={font="\\large"},
#     x_label_style={font="\\Large"},
#     y_label_style={font="\\Large"},
#     # xmin=-5, xmax=5,
#     # xmode="log",
#     ymode="log",
#     # enlarge_x_limits=false,
#     #ymin=0.5,ymax=1,
#     "after end axis/.code" = {
#         raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont, fill=white] at (rel axis cs:.02,.98) {(d)};"
#     },

#     },
#     Plot({color="red",mark="square", mark_repeat=5, only_marks}, Table((τs), 𝔼ts)),
#     Plot({color="red",mark="square", mark_repeat=5, no_marks}, Table((τs), 𝔼̂ts)),
#     Plot({color="blue",dashed,ultra_thick}, Table((τs), Tꜛs)),
#     # LegendEntry("Channel Capacity"),
#     LegendEntry(L"\mathbb{E}[t]"),
#     LegendEntry(raw"$\mathbb{E}[t]$ approximate" ),
#     LegendEntry(L"T^*"),
# )

plot3 = @pgf Axis(
    {width="3.4in", 
    height="3.4in", 
    xlabel=raw"inaccuracy $1 - \mathcal{A}$", 
    ylabel=raw"total duration",
    legend_pos="north east",
    legend_style={font="\\large", fill="none", draw="none",at="(0.95,0.95)", anchor="north east"},
    legend_cell_align="left",
    # large font size and tick label style
    tick_label_style={font="\\large"},
    x_label_style={font="\\Large"},
    y_label_style={font="\\Large"},
    # xmin=-5, xmax=5,
    xmode="log",
    ymode="log",
    # enlarge_x_limits=false,
    #ymin=0.5,ymax=1,
    "after end axis/.code" = {
        raw"\node[anchor=north west, font=\fontsize{12}{14}\selectfont, fill=white] at (rel axis cs:.02,.98) {(d)};"
    },

    },
    # Plot({color="red",mark="square", mark_repeat=5, no_marks}, Table((τs), 𝔼̂ts)),
    Plot({color="blue", mark="o", mark_repeat=5}, Table((DNA_inaccuracies), Tꜛs)),
    LegendEntry("DNA replication"),
    # LegendEntry("Channel Capacity"),
    Plot({color="red",mark="square", mark_repeat=5}, Table(1 .-(Accuracy⁺), 𝔼ts)),
    LegendEntry("TCR signaling"),
)

push!(tikzpic, plot3)
pgfsave("FIG/multiround_total_duration_vs_accuracy.pgf", tikzpic)
pgfsave("FIG/multiround_total_duration_vs_accuracy.pdf", tikzpic)
pgfsave("FIG/multiround_total_duration_vs_accuracy.svg", tikzpic)