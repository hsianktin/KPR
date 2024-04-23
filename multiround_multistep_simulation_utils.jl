using Random
using Statistics
using Parameters

# setup Parameters
@with_kw struct Pars
    k₁::Float64 = 0.1
    k₋₁::Float64 = 1
    τ::Float64 = 3 # mean delay time 
    N_f::Int = 10 # number of proofreading steps
    N_s::Int = 10 # threshold for the stabilization step 
    k₋₁⁺::Float64 = 1
    kₚ::Float64 = 0.01
    T::Float64 = 100
end

# ϕ = [n, X], n = number of signaling events, X = number of phosphorylated sites. X=0 means the kinase is unbound. X > 0 means the kinase is bound.
# setup transitions
function binding!(ϕ)
    if ϕ[2] == 0
        ϕ[2] = 1
    else
        @error "binding! should only be called when X = 0"
    end
end

function unbinding!(ϕ)
    if ϕ[2] > 0
        ϕ[2] = 0
    else
        @error "unbinding! should only be called when X > 0"
    end
end

function activation!(ϕ)
    ϕ[2] += 1
end

function signaling!(ϕ)
    ϕ[1] += 1
end

### Because deterministic delay is involved in this system,
### we introduce a queue to store the events that will happen in the future.
events = []


# setup stepping function
function gillespie!(ϕ, t; pars=Pars())
    # @show ϕ
    reactions = [binding!, unbinding!, activation!, signaling!]
    rates = zeros(4)
    X = ϕ[2] # current state
    if X == 0
        rates[1] = pars.k₁ # binding rate
    elseif X ≥ 1 && X ≤ pars.N_f
        rates[3] = pars.N_f / pars.τ # activation rate
    elseif X == pars.N_f + 1
        rates[4] = pars.kₚ # signaling rate
    end
    if X > pars.N_s
        rates[2] = pars.k₋₁⁺
    elseif X ≥ 1
        rates[2] = pars.k₋₁
    end
    # @show rates
    # find reactions with non-zero rates
    reactions = reactions[rates .> 0]
    rates = rates[rates .> 0]
    events = []
    # push random waiting time for each reaction to the queue
    for (i, r) in enumerate(reactions)
        push!(events, [randexp()/rates[i], r])
    end
    δt, i = findmin([e[1] for e in events])
    # update state
    events[i][2](ϕ)
    t += δt
    # @show ϕ
    return t
end

# setup simulation function
# function simulation(ϕ₀, pars=Pars(), times=[])
#     ϕ = deepcopy(ϕ₀)
#     t = 0.0
#     τₐ = +Inf # first activation time
#     𝐓ₚ = [+Inf for i in 1:100]
#     function update_FPT!(𝐓ₚ,ϕ,t)
#         nₚ = ϕ[1]
#         if 0 < nₚ ≤ length(𝐓ₚ)
#             if 𝐓ₚ[nₚ] == +Inf
#                 𝐓ₚ[nₚ] = t
#             end
#         end
#     end
#     Φ = [deepcopy(ϕ)]
#     T = [t]
#     while t < pars.T
#         t = gillespie!(ϕ, t; pars=pars) # inspect current state, update events, and react
#         if ϕ[2] == pars.N_f + 1 && τₐ == +Inf # record first activation time
#             τₐ = t
#         end
#         update_FPT!(𝐓ₚ,ϕ,t)
#         push!(Φ, deepcopy(ϕ))
#         push!(T, t)
#     end
#     # the final state t ≥ T, which we don't need
#     pop!(Φ)
#     pop!(T)
#     push!(T, pars.T)
#     push!(Φ, deepcopy(Φ[end]))
#     if length(times) > 0
#         # align the time points
#         Φ′ = [
#                 Φ[
#                     maximum(findall(T .<= t))
#                 ] # find the last state before t
#                 for t in times
#             ]
#         @assert times[end] <= pars.T
#         # manual GC
#         Φ = []
#         T = []
#         return Φ′, times, τₐ, 𝐓ₚ
#     else
#         return Φ, T, τₐ, 𝐓ₚ
#     end
# end

function simulation(ϕ₀, pars=Pars(), times=[]; dt=0.1)
    ϕ = deepcopy(ϕ₀)
    t = 0.0
    τₐ = +Inf # first activation time
    𝐓ₚ = [+Inf for i in 1:100]
    function update_FPT!(𝐓ₚ,ϕ,t)
        nₚ = ϕ[1]
        if 0 < nₚ ≤ length(𝐓ₚ)
            if 𝐓ₚ[nₚ] == +Inf
                𝐓ₚ[nₚ] = t
            end
        end
    end
    Φ = [deepcopy(ϕ)]
    T = [t]
    while t < pars.T
        t = gillespie!(ϕ, t; pars=pars) # inspect current state, update events, and react
        if ϕ[2] == pars.N_f + 1 && τₐ == +Inf # record first activation time
            τₐ = t
        end
        update_FPT!(𝐓ₚ,ϕ,t)
        if floor(Int,t / dt) > floor(Int,T[end] / dt)
            push!(Φ, deepcopy(ϕ))
            push!(T, t)
        end
        # push!(Φ, deepcopy(ϕ))
        # push!(T, t)
    end
    # the final state t ≥ T, which we don't need
    pop!(Φ)
    pop!(T)
    push!(T, pars.T)
    push!(Φ, deepcopy(Φ[end]))
    if length(times) > 0
        # align the time points
        Φ′ = [
                Φ[
                    maximum(findall(T .<= t))
                ] # find the last state before t
                for t in times
            ]
        # @assert times[end] <= pars.T
        # manual GC
        Φ = []
        T = []
        return Φ′, times, τₐ, 𝐓ₚ
    else
        return Φ, T, τₐ, 𝐓ₚ
    end
end

# setup simulation state variables (n,X)
ϕ₀ = [0,1]

# basic saliency test
# setup simulation parameters
# pars = Pars(
#     k₁ = 1,
#     k₋₁ = 0,
#     τ = 10,
#     N_f = 10,
#     N_s = 10,
#     k₋₁⁺ = 0,
#     kₚ = 1,
#     T = 1000
# )


# k_pars = Pars(
#     k₁ = 0.1,
#     k₋₁ = 1,
#     τ = 3,
#     N_f = 2000,
#     N_s = 2000,
#     k₋₁⁺ = 1,
#     kₚ = 0.01,
#     T = 1000
# )

# q_pars = Pars(
#     k₁ = 0.1,
#     k₋₁ = 2,
#     τ = 3,
#     N_f = 2000,
#     N_s = 2000,
#     k₋₁⁺ = 2,
#     kₚ = 0.01,
#     T = 1000
# )

# run simulation
# using Plots
# τₐs = []
# for _ in 1:1000
# Φ, T, τₐ, 𝐓ₚ = simulation(ϕ₀, pars)
# push!(τₐs, τₐ)
# end
# histogram(τₐs, bins=100, label="τₐ",title="N=10") 
# savefig("N=10.png")
# pars = Pars(
#     k₁ = 1,
#     k₋₁ = 0,
#     τ = 10,
#     N_f = 100,
#     N_s = 100,
#     k₋₁⁺ = 0,
#     kₚ = 10,
#     T = 100
# )
# τₐs = []
# for _ in 1:1000
# Φ, T, τₐ, 𝐓ₚ = simulation(ϕ₀, pars)
# push!(τₐs, τₐ)
# end
# histogram(τₐs, bins=100, label="τₐ",title="N=100") 
# savefig("N=100.png")

# # theoretical expectatation
# ta = []
# for _ in 1:1000
#     push!(ta, sum([randexp()/1  for _ in 1:10]))
# end
# histogram(ta, bins=100, label="τₐ",title="N=100")
# # sum([randexp()/10  for _ in 1:100])




# # # plot
# using Plots
# # plot(T, [ϕ[1] for ϕ in Φ], label="n")
# # # # ##### validation for first passage time #####
# # Φ̃ = [ϕ[2] for ϕ in Φ]
# # plot(T[2:end], Φ̃[2:end], xaxis=:log)
# # # # mark the first activation time τₐ on plot 
# # plot!([τₐ], [pars.N_f+1], label="τₐ",marksize=5, mark=:circle, markercolor=:red)
# # # # visual inspection of the first passage time confirms the FPT is indeed τₐ
# # τₐs = []
# # using ProgressMeter
# # Φ, T, τₐ, 𝐓ₚ = simulation(ϕ₀, pars)
# begin
# k_τₐs = []
# q_τₐs = []
# @showprogress for i in 1:1000
#     Φ, T, τₐ, 𝐓ₚ = simulation(ϕ₀, k_pars)
#     push!(k_τₐs, τₐ)
#     Φ, T, τₐ, 𝐓ₚ = simulation(ϕ₀, q_pars)
#     push!(q_τₐs, τₐ)
# end
# histogram(k_τₐs, bins=100, label="τₐ")
# histogram!(q_τₐs, bins=100, label="τₐ")
# end

# response_k = [τₐ == +Inf ? 0 : 1 for τₐ in k_τₐs]
# response_q = [τₐ == +Inf ? 0 : 1 for τₐ in q_τₐs]
# P_response_k = sum(response_k)/length(response_k)
# P_response_q = sum(response_q)/length(response_q)
# P11 = P_response_k
# P10 = P_response_q
# P01 = 1 - P_response_k
# P00 = 1 - P_response_q
# P_1 = 1/2
# P_0 = 1/2
# P1 = P11 * P_1 + P10 * P_0
# P0 = P01 * P_1 + P00 * P_0
# # compute the mutual information
# I = (P11 * P_1 * log(P11/(P1)) + P10 * log(P10/P1) + P01 * log(P01/P0) + P00 * log(P00/P0))
# # ##### validation for steady state distribution of states #####
# # # # find culumative dwell time in each state
# # Ts = [sum([diff(T)[i] for i in findall(Φ̃[1:end-1] .== i)]) for i in 0:pars.N_f+1]
# # 𝐏 = Ts/pars.T
# # 𝐏_normalized = [𝐏[1], sum(𝐏[2:end-1]), 𝐏[end]]
# # # 𝐏̂ = [T_0, T_1, T_2]/pars.T
# # # # theoretical distribution
# # weight = [1,
# #          (pars.k₁/pars.k₋₁) * (1 - exp(- pars.k₋₁ * pars.τ)), 
# #          (pars.k₁/pars.k₋₁⁺) * exp(- pars.k₋₁⁺ * pars.τ)]
# # 𝐏̃ = weight ./ sum(weight)
