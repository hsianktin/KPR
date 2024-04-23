using Statistics

function get_η_FLD(results_k, results_q)
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

    # η_FLD = |𝔼nₜ∣q - 𝔼nₜ∣k|/(√𝕍nₜ∣q + √𝕍nₜ∣k)
    η_FLD = [abs(𝔼nₜ❘q[i] - 𝔼nₜ❘k[i]) / (sqrt(𝕍nₜ❘q[i]) + sqrt(𝕍nₜ❘k[i])) for i in eachindex(𝔼nₜ❘q)]
    # replace NaN with 0
    η_FLD = [ifelse(isnan(x), 0, x) for x in η_FLD]
    return η_FLD
end