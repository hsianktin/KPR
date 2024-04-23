begin 
    using ProgressMeter
    using BSON
end

substrate = ARGS[1]
id = ARGS[2]
# profile is either "PROFILES/base.jl" or ARGS[3], depending on
# whether we are provided with ARGS[3]
profile = length(ARGS) > 2 ? ARGS[3] : "PROFILES/base.jl"
# get base name of profile
label = split(profile, "/")[end] |> x -> split(x, ".jl")[1]
# remove .jl extension
begin
    include("multiround_simulation_utils.jl");

    include(profile);
end

function task_k()
    # run simulation
    Φ, T, τₐ, 𝐓ₚ = simulation(ϕ₀, k_pars, times)
    return Φ, T, τₐ, 𝐓ₚ
end

function task_q()
    # run simulation
    Φ, T, τₐ, 𝐓ₚ = simulation(ϕ₀, q_pars, times)
    return Φ, T, τₐ, 𝐓ₚ
end

N_samples = 100
if substrate == "k"
    task = task_k
else
    task = task_q
end
if !isfile("DATA/multiround_simulation_results_$(substrate)_$(id)_$(label).bson")
    results = map(
        x -> task(),
        1:N_samples;
    )

    BSON.@save "DATA/multiround_simulation_results_$(substrate)_$(id)_$(label).bson" results
else 
    @info "File already exists, skipping"
end
