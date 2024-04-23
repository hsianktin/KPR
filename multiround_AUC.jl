using QuadGK
# function AUC(N::Int; α=0.5)
#     total_sum = 0
#     for j in 0:N
#         for k in 0:N-1
#             total_sum += binomial(N, j) * binomial(N - 1, k) * (-1)^(j+k) / (α * j + k+1)
#         end
#     end
#     return 1 - N * total_sum
# end

function AUC_int(N; α=0.5)
    # (1 - (1 - y**alpha)**N) * N / alpha * (1 - y)**(N/alpha - 1)
    integrand(y) = (1 - (1 - y^α)^N) * N * (1 - y)^(N - 1)
    return quadgk(integrand, 0, 1)[1]
end

Ns = 1:100
AUCs = [AUC_int(N) for N in Ns]
using PGFPlotsX

plot = @pgf Axis({
        xlabel = "N",
        ylabel = "AUC",
        width = "3.4in",
        height = "3.4in",
        legend_pos = "south east",
        xmin=0, xmax=100,
        ymin=0.5,ymax=1,
        label_style = {font = "\\Large"},
    },
    PlotInc(
        {
            color = "red!80!black",
            mark = "none",},
        Table(
             collect(Ns),
             AUCs,
        ),
    ),
)
pgfsave("./FIG/multiround_AUC.pgf", plot)