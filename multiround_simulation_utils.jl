using Random
using Statistics
using Parameters

# setup Parameters
@with_kw struct Pars
    k₁::Float64 = 0.1
    k₋₁::Float64 = 1
    τ::Float64 = 3 # deterministic delay
    k₋₁⁺::Float64 = 1
    kₚ::Float64 = 0.01
    T::Float64 = 100
end


# setup transitions
function binding!(ϕ)
    if ϕ[2] == :E₊S
        ϕ[2] = :ES
    end
end

function unbinding!(ϕ)
    ϕ[2] = :E₊S
end

function activation!(ϕ)
    if ϕ[2] == :ES
        ϕ[2] = :ES⁺
    end
end

function signaling!(ϕ)
    if ϕ[2] == :ES⁺
        ϕ[1] += 1
    end
end

### Because deterministic delay is involved in this system,
### we introduce a queue to store the events that will happen in the future.
events = []


# setup stepping function
function gillespie!(ϕ, events, t; pars=Pars())
    reactions = [binding!, unbinding!, signaling!]
    rates = zeros(3)
    X = ϕ[2]
    if X == :E₊S
        rates[1] = pars.k₁
    elseif X == :ES
        rates[2] = pars.k₋₁
    elseif X == :ES⁺
        rates[2] = pars.k₋₁⁺
        rates[3] = pars.kₚ
    end
    # find reactions with non-zero rates
    reactions = reactions[rates .> 0]
    rates = rates[rates .> 0]
    # push random waiting time for each reaction to the queue
    for (i, r) in enumerate(reactions)
        push!(events, [randexp()/rates[i], r])
    end
    # check if activation! event is in the queue
    is_waiting_for_activation = false
    for e in events
        if e[2] == activation!
            is_waiting_for_activation = true
            break
        end
    end
    if X == :ES
        if !is_waiting_for_activation # this means ES is newly formed and no activation! event is in the queue
            push!(events, [pars.τ, activation!])
        end
    else
        # remove activation! event from the queue
        deleteat!(events, findall(e -> e[2] == activation!, events))
    end
    # find the next event
    δt, i = findmin([e[1] for e in events])
    # validate no other event happens at the same time
    for (j, e) in enumerate(events)
        if j != i && e[1] == δt
            error("two events happen at the same time")
        end
    end
    # react!
    # @show t, ϕ, events[i][2]
    events[i][2](ϕ)
    t += δt
    # remove the event from the queue
    deleteat!(events, i)
    # for the rest of the events, update their waiting time
    for e in events
        e[1] -= δt
    end
    # forget all events except activation! event due to memoryless property
    deleteat!(events, findall(e -> e[2] != activation!, events))
    return t
end

# setup simulation function
function simulation(ϕ₀, pars=Pars(), times=[])
    ϕ = deepcopy(ϕ₀)
    t = 0.0
    events = []
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
        t = gillespie!(ϕ, events, t; pars=pars) # inspect current state, update events, and react
        if ϕ[2] == :ES⁺ && τₐ == +Inf # record first activation time
            τₐ = t
        end
        update_FPT!(𝐓ₚ,ϕ,t)
        push!(Φ, deepcopy(ϕ))
        push!(T, t)
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
        @assert times[end] <= pars.T
        # manual GC
        Φ = []
        T = []
        return Φ′, times, τₐ, 𝐓ₚ
    else
        return Φ, T, τₐ, 𝐓ₚ
    end
end

# setup simulation state variables (n,X)
ϕ₀ = [0,:E₊S]

# setup simulation parameters
# pars = Pars(
#     k₁ = 0.1,
#     k₋₁ = 1,
#     τ = 3,
#     k₋₁⁺ = 1,
#     kₚ = 0.01,
#     T = 100000
# )

# run simulation
# Φ, T, τₐ = simulation(ϕ₀, pars)

# plot
# using Plots
# plot(T, [ϕ[1] for ϕ in Φ], label="n")
# ##### validation for first passage time #####
# Φ̃ = [ϕ[2] == :E₊S ? 0 : ϕ[2] == :ES ? 1 : 2 for ϕ in Φ]
# plot(T[2:end], Φ̃[2:end], xaxis=:log)
# # mark the first activation time τₐ on plot 
# plot!([τₐ], [2], label="τₐ",marksize=5, mark=:circle, markercolor=:red)
# visual inspection of the first passage time confirms the FPT is indeed τₐ

##### validation for steady state distribution of states #####
# # find culumative dwell time in each state
# T_0 = sum([diff(T)[i] for i in findall(Φ̃[1:end-1] .== 0)])
# T_1 = sum([diff(T)[i] for i in findall(Φ̃[1:end-1] .== 1)])
# T_2 = sum([diff(T)[i] for i in findall(Φ̃[1:end-1] .== 2)])
# 𝐏̂ = [T_0, T_1, T_2]/pars.T
# # theoretical distribution
# weight = [1,
#          (pars.k₁/pars.k₋₁) * (1 - exp(- pars.k₋₁ * pars.τ)), 
#          (pars.k₁/pars.k₋₁⁺) * exp(- pars.k₋₁⁺ * pars.τ)]
# 𝐏 = weight ./ sum(weight)
