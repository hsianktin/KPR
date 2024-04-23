using PGFPlotsX
using SpecialFunctions
q₁ = k₁ = 0.1
k₋₁ = 1
τ = 5
k₋₁⁺ = 1
kₚ = 0.01
T = 1e6

q₋₁s = [i for i in 0:0.01:3]

A(q₋₁) = k₁/(k₁ + k₁ * exp(-(q₋₁ - k₋₁) * τ) )

As = A.(q₋₁s)

tikzpicture = @pgf TikzPicture({baseline})
axis = @pgf Axis({
    width = "3in",
    height = "3in",
    xlabel = "unbinding rate of incorrect substrate \$q_{-1}\$ (s\$^{-1}\$)",
    ylabel = "accuracy",
    xmin = 0,
    xmax = 3,
    ymin = 0,
    ymax = 1,
    legend_style = "{fill=none, draw=none, at={(1,0)}, anchor=south east}",
    legend_cell_align = "left",
},
Plot({mark = "none", color = "blue"}, Table([q₋₁s, As])),
)


N_opt(q₋₁) = abs(q₋₁ - k₋₁) * τ * exp(minimum([k₋₁,q₋₁]) * τ) / (1 - exp(-abs(q₋₁ - k₋₁) * τ))
A_FPT(q₋₁) = 1/2 * (1 + (1 - exp(- q₋₁ * τ))^N_opt(q₋₁) - (1 - exp(- k₋₁ * τ))^N_opt(q₋₁))
A_p(q₋₁) = 1/2 * (1 + erf(sqrt(1/2*kₚ*T*k₁/(k₁+ k₋₁)*exp(-k₋₁ * τ))-sqrt(1/2*kₚ*T*q₁/(q₁+ q₋₁)*exp(-q₋₁ * τ))))

A_FPTs = A_FPT.(q₋₁s)
A_ps = A_p.(q₋₁s)

push!(axis, @pgf Plot(
    {mark = "none", color = "red"},
    Table([q₋₁s, A_FPTs])
))

push!(axis, @pgf Plot(
    {mark = "none", color = "green"},
    Table([q₋₁s, A_ps])
))

push!(axis, @pgf LegendEntry(
    "DNA"
))
push!(axis, @pgf LegendEntry(
    "FPT in TCR"
))
push!(axis, @pgf LegendEntry(
    "product in TCR"
))
axis
push!(tikzpicture, axis)

pgfsave("FIG/response_plot_dna_accuracy.tex", tikzpicture)

