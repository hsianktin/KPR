using Integrals,Optim
∑ = sum

function mutual_information_continuous(ρ₀,ρ₁,p=0.5)
    ρ(x) = p*ρ₀(x) + (1-p)*ρ₁(x)
    function h₀(x)
        if x ≤ 0
            return 0
        else
            result = p* ρ₀(x) * log((p*ρ₀(x))/(p*ρ(x)))
            if isnan(result)
                return 0
            else
                return result
            end
        end
    end
    function h₁(x)
        if x ≤ 0
            return 0
        else
            result = (1-p)* ρ₁(x) * log(((1-p)*ρ₁(x))/((1-p)*ρ(x)))
            if isnan(result)
                return 0
            else
                return result
            end
        end
    end
    prob₀ = IntegralProblem((u,ϕ) -> h₀(u), -Inf, Inf)
    sol₀ = solve(prob₀, HCubatureJL())
    H_ρ₀ = sol₀[1]
    prob₁ = IntegralProblem((u,p) -> h₁(u), -Inf, Inf)
    sol₁ = solve(prob₁, HCubatureJL())
    H_ρ₁ = sol₁[1]
    return (H_ρ₀ + H_ρ₁)/log(2)
end

using Plots
xs = collect(-10:0.1:1000)
# plot(xs,h₀.(xs),label="ρ₀")
k₋₁ = 1
q₋₁ = 2
function ρ₀(x)
    if x ≤ 0
        return 0
    else
        return k₋₁ * exp(- k₋₁ * x) # exponential distribution
    end
end

function ρ₁(x)
    if x ≤ 0
        return 0
    else
        return q₋₁ * exp(- q₋₁ * x) # exponential distribution
    end
end

mutual_information_continuous(ρ₀,ρ₁,0.5)

function mutual_information_binary(P₀,P₁,p=0.5)
    P(x) = p*P₀(x) + (1-p)*P₁(x)
    function h₀(x)
        if x ≤ 0
            return 0
        else
            result = p* P₀(x) * log((p*P₀(x))/(p*P(x)))
            if isnan(result)
                return 0
            else
                return result
            end
        end
    end
    function h₁(x)
        if x ≤ 0
            return 0
        else
            result = (1-p)* P₁(x) * log(((1-p)*P₁(x))/((1-p)*P(x)))
            if isnan(result)
                return 0
            else
                return result
            end
        end
    end
    Ξ = [0,1]
    prob₀ = ∑(h₀.(Ξ))
    prob₁ = ∑(h₁.(Ξ))
    return (prob₀ + prob₁)/log(2)
end

function P₀(x, Th = 2)
    prob = IntegralProblem((u,ϕ) -> ρ₀(u), (0, Th), dim=1)
    sol = solve(prob, HCubatureJL(), reltol=1e-10)
    function filter(x)
        if x ≤ 0
            return 1e-10
        elseif x ≥ 1
            return 1- 1e-10
        else
            return x
        end
    end
    if isnan(sol[1])
        @error "sol[1] is NaN"
    end
    if x == 0
        return sol[1] |> filter
    else
        return 1- sol[1] |> filter
    end  
end

function P₁(x, Th = 2)
    prob = IntegralProblem((u,p) -> ρ₁(u), (Th, Inf), dim=1)
    sol = solve(prob, HCubatureJL(), reltol=1e-10)
    function filter(x)
        if x ≤ 0
            return 1e-10
        elseif x ≥ 1
            return 1- 1e-10
        else
            return x
        end
    end
    if x == 0
        return 1 - sol[1] |> filter
    else
        return sol[1] |> filter
    end    
end
# there is also another stochastic way for binary output 
function P₀′(x, k_f = 1)
    if x == 0
        return q₋₁/(q₋₁ + k_f)
    else
        return 1 - q₋₁/(q₋₁ + k_f)
    end
end

function P₁′(x, k_f = 1)
    if x == 0
        return k₋₁/(k₋₁ + k_f)
    else
        return 1 - k₋₁/(k₋₁ + k_f)
    end
end

Ths = collect(0.1:0.1:100)
MIbs = [mutual_information_binary(x->P₀(x,Th),
    x->P₁(x,Th),
0.5) for Th in Ths]
MIbs′ = [mutual_information_binary(x->P₀′(x,1/Th),
    x->P₁′(x,1/Th),
0.5) for Th in Ths]
 [mutual_information_binary(x->P₀′(x,1/Th),
    x->P₀′(x,1/Th),
0.5) for Th in Ths]

using PGFPlotsX
# add preamble using dsfont package
push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\usepackage{dsfont,amssymb}")
plot1 = @pgf Axis({
    xlabel = "mean threshold time \$\\mathbb{E}[\\tau]\$",
    ylabel = "mutual information",
    legend_pos = "north west",
    width = "3.4in",
    height = "3.4in",
    xmin = 0,
    xmax = 100,
    ymin = 0,
    ymax = 0.05,
    xmode = "log",
    ymode = "linear",
    legend_entries = {"\$\\mathcal{I}(\\xi;\\mathds{1}_{a> \\tau}\$)",
        "\$\\mathcal{I}(\\xi;\\mathds{1}_{a>\\tau_f})\$"},
    legend_columns = 2,
    legend_style = {draw = "none", fill = "none", at = "{(0.02,0.98)}", anchor = "north west"},
    },
    PlotInc({mark="none",color="blue"}, Table(Ths, MIbs)),
    PlotInc({mark="none",color="red"}, Table(Ths, MIbs′)),
)
pgfsave("FIG/mutual_information_binary.pgf", plot1)
pgfsave("FIG/mutual_information_binary.pdf", plot1)
pgfsave("FIG/mutual_information_binary.svg", plot1)

function supMIB(P₀,P₁)
        
    # use Optim to find optimal Th
    function loss(Th)
        return -mutual_information_binary(x->P₀(x,Th[1]),
            x->P₁(x,Th[1]),
        0.5)
    end
    initial_guess = [0.03]
    lower_bounds = [0.0]  # Lower bound for p₀
    upper_bounds = [100.0]  # Upper bound for p₀
    opt_result = optimize(loss, lower_bounds, upper_bounds ,initial_guess, Fminbox(LBFGS()); inplace = false)

    opt_result.minimizer[1]
    return -opt_result.minimum
end
supMIB(P₀,P₁)

supMIB(P₀′,P₁′)
P₀′(1,1)
P₁′(1,1)