#!/usr/bin/env julia

using Distributions, DataFrames

# Fit a log-normal mixture from isolator transcript quantifications.

# Fit a normal mixture consisting of k components over observations xs.
function fit_norm_mix(xs::AbstractVector{Float64}, k::Integer)

    n = length(xs)
    ps = fill(1.0/k, k)
    zs = Array(Float64, (k, n))
    mus = Array(Float64, k)
    sds = Array(Float64, k)
    work = Array(Float64, n)

    xs_mu = mean(xs)
    xs_sd = std(xs)

    fill!(sds, xs_sd)

    # choose evenly spaced starting components spanning three sds to either side
    # of the mean
    if k == 1
        mus[1] = xs_mu
    else
        for i in 1:k
            mus[i] = xs_mu + (-1 + 2 * ((i-1)/(k-1))) * xs_sd
        end
    end

    eps = 1e-3
    ll0 = norm_mix_log_likelihood(xs, work, ps, mus, sds)
    it = 1
    while true
        fit_norm_mix_e_step(xs, ps, zs, mus, sds)
        fit_norm_mix_m_step(xs, ps, zs, mus, sds)

        ll = norm_mix_log_likelihood(xs, work, ps, mus, sds)
        #println("Iteration: ", it)
        #println("Log-likelihood: ", ll)
        #print("Mu: ")
        #print(mus)
        #println("Sd: ")
        #print(sds)
        #println("Ps: ")
        #print(ps)
        #println()

        if abs(ll - ll0) < eps
            break
        end

        it += 1
        ll0 = ll
    end

    return (ps, mus, sds)
end


function norm_mix_log_likelihood(xs::AbstractVector{Float64},
                                 work::AbstractVector{Float64},
                                 ps::AbstractVector{Float64},
                                 mus::AbstractVector{Float64},
                                 sds::AbstractVector{Float64})
    n = length(xs)
    k = length(ps)
    fill!(work, 0.0)
    for i in 1:k
        dist = Normal(mus[i], sds[i])
        for j in 1:n
            work[j] += ps[i] * pdf(dist, xs[j])
        end
    end

    ll = 0.0
    for j in 1:n
        ll += log(work[j])
    end

    return ll
end


function fit_norm_mix_e_step(xs::AbstractVector{Float64},
                             ps::AbstractVector{Float64},
                             zs::AbstractMatrix{Float64},
                             mus::AbstractVector{Float64},
                             sds::AbstractVector{Float64})
    k = length(ps)
    n = length(xs)

    for i in 1:k
        dist = Normal(mus[i], sds[i])
        for j in 1:n
            zs[i, j] = ps[i] * pdf(dist, xs[j])
        end
    end

    for j in 1:n
        zsum = 0.0
        for i in 1:k
            zsum += zs[i, j]
        end
        for i in 1:k
            zs[i, j] /= zsum
        end
    end
end

function fit_norm_mix_m_step(xs::AbstractVector{Float64},
                             ps::AbstractVector{Float64},
                             zs::AbstractMatrix{Float64},
                             mus::AbstractVector{Float64},
                             sds::AbstractVector{Float64})
    k = length(ps)
    n = length(xs)

    for i in 1:k
        zsum = 0.0
        mus[i] = 0.0
        for j in 1:n
            zsum += zs[i, j]
            mus[i] += zs[i, j] * xs[j]
        end
        mus[i] /= zsum

        sds[i] = 0.0
        for j in 1:n
            sds[i] += zs[i, j] * (xs[j] - mus[i])^2
        end
        sds[i] = sqrt(sds[i] / zsum)
        ps[i] = zsum / n
    end
end



# Fitting to cufflinks estimates. The big problem here is that cufflinks will output
# a ton of zeros which we can't fit log-normals to. If we discard the zeros we will
# get parameters that exaggerate the avg abundance. We can replace the zeros with
# some small value


function readfpkm(filename)
    df = readtable(filename, header=true, separator='\t')
    xs = df[:FPKM][sortperm(df[:tracking_id])]
    # Cufflinks outputs tons of zeros which won't work for fitting
    # parameters on a log scale. We replace those zeros with whatever
    # the smallest non-zero reported value is
    minval = Inf
    for x in xs
        if x != 0 && x > 1e-6 && x < minval
            minval = x
        end
    end
    minval *= 0.1

    xs[xs .< minval] = minval

    return log(1e6 * xs ./ sum(xs))
end


# upper quantile norm
function uqn(xs)
    return xs ./ quantile(xs, 0.8)
end


a1 = readfpkm("training/4-h7-day20.isoforms.fpkm_tracking")
a2 = readfpkm("training/5-h7-day20.isoforms.fpkm_tracking")
a3 = readfpkm("training/6-h7-day20.isoforms.fpkm_tracking")

b1 = readfpkm("training/1-h7-1yr.isoforms.fpkm_tracking")
b2 = readfpkm("training/2-h7-1yr.isoforms.fpkm_tracking")
b3 = readfpkm("training/3-h7-1yr.isoforms.fpkm_tracking")

a_tpm = (a1 + a2 + a3) ./ 3
b_tpm = (b1 + b2 + b3) ./ 3
experiment_tpm = (a_tpm + b_tpm) ./ 2

println("Experiment Params:")
ps, mus, sds = fit_norm_mix(experiment_tpm, 2)
println("p: ", join(ps, ", "))
println("mu: ", join(mus, ", "))
println("sd: ", join(sds, ", "))

# We can't trust the fold-change estimates from cufflinks for
# low abundance transcripts.
idxs = experiment_tpm .>= quantile(experiment_tpm, 0.90)

println("Condition Params:")
#condition_tpm = vcat(a1[idxs] - experiment_tpm[idxs],
                     #a2[idxs] - experiment_tpm[idxs],
                     #a3[idxs] - experiment_tpm[idxs],
                     #b1[idxs] - experiment_tpm[idxs],
                     #b2[idxs] - experiment_tpm[idxs],
                     #b3[idxs] - experiment_tpm[idxs])
condition_tpm = vcat(a_tpm[idxs] - experiment_tpm[idxs],
                     b_tpm[idxs] - experiment_tpm[idxs])

ps, mus, sds = fit_norm_mix(condition_tpm, 2)
println("p: ", join(ps, ", "))
println("mu: ", join(mus, ", "))
println("sd: ", join(sds, ", "))

println("Replicate Params")
replicate_tpm = vcat(a1[idxs] - a_tpm[idxs],
                     a2[idxs] - a_tpm[idxs],
                     a3[idxs] - a_tpm[idxs],
                     b1[idxs] - b_tpm[idxs],
                     b2[idxs] - b_tpm[idxs],
                     b3[idxs] - b_tpm[idxs])

ps, mus, sds = fit_norm_mix(replicate_tpm, 2)
println("p: ", join(ps, ", "))
println("mu: ", join(mus, ", "))
println("sd: ", join(sds, ", "))


