#!/usr/bin/env julia

using Distributions

# This is a simple tool to aid in generating expression values which are then
# into flux simulator to generate reads.

function read_gtf_transcript_ids(input::IO)
    pat = r"transcript_id\s*\"?([^\"\t\n]+)\"?"
    data = readall(input)
    transcript_ids = Set()
    for m in eachmatch(pat, data)
        push!(transcript_ids, m.captures[1])
    end
    return transcript_ids
end


function main()
    if length(ARGS) < 4
        println("Usage: seqsim-generate.jl genes.gtf p1 mu1 sd1 [[p2 mu2 sd2]...]")
        return
    end

    transcript_ids = read_gtf_transcript_ids(open(ARGS[1]))
    n = length(transcript_ids)
    println(STDERR, "simulating expression for ", n, " transcripts")

    ps  = Array(Float64, 0)
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

    println("transcript_id\ttpm")
    component_dist = Categorical(ps)
    dists = [LogNormal(mu, sd) for (mu, sd) in zip(mus, sds)]
    xs = Array(Float64, n)
    for (j, transcript_id) in enumerate(transcript_ids)
        xs[j] = rand(dists[rand(component_dist)])
    end
    xs ./= sum(xs)

    for (transcript_id, x) in zip(transcript_ids, xs)
        println(transcript_id, "\t", 1e6 * x)
    end
end


main()


