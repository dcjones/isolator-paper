#!/usr/bin/env julia

using Distributions, DataFrames

function main()
    if length(ARGS) < 4
        println("Usage: seqsim-perturb.jl expression.tsv p1 mu1 sd1 [[lp2 mu2 sd2]...]")
        return
    end

    base_expr = readtable(ARGS[1], separator='\t', header=true)
    n = size(base_expr, 1)

    ps = Array(Float64, 0)
    mus = Array(Float64, 0)
    sds = Array(Float64, 0)
    for i in 2:3:length(ARGS)
        push!(ps,  parsefloat(ARGS[i]))
        push!(mus, parsefloat(ARGS[i+1]))
        push!(sds, parsefloat(ARGS[i+2]))
    end
    if !(length(ps) == length(mus) == length(sds))
        error("A matching number of mu and sigma parameters must be given.")
    end
    ps ./= sum(ps)

    xs = Array(Float64, n)
    component_dist = Categorical(ps)
    dists = [LogNormal(mu, sd) for (mu, sd) in zip(mus, sds)]
    for j in 1:n
        d = rand(dists[rand(component_dist)])
        xs[j] = base_expr[j, :tpm] * d
    end
    xs ./= sum(xs)

    println("transcript_id\ttpm")
    for (transcript_id, x) in zip(base_expr[:transcript_id], xs)
        println(transcript_id, "\t", 1e6 * x)
    end
end


main()


