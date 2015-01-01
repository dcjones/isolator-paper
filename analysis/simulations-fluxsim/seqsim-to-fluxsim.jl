#!/usr/bin/env julia

using DataFrames

function main()
    if length(ARGS) != 3
        println("Usage: seqsim-to-fluxsim.jl expression.tsv output.pro N")
        return
    end

    metadata = Dict()
    for line in eachline(open(ARGS[2]))
        row = split(strip(line), '\t')
        metadata[row[2]] = (row[1], row[3], row[4])
    end

    ex = readtable(ARGS[1], header=true, separator='\t')
    transcript_count = parseint(ARGS[3])

    # rewrite the .pro file
    output = open(ARGS[2], "w")

    for i in 1:size(ex, 1)
        transcript_id = ex[:transcript_id][i]

        if !haskey(metadata, transcript_id)
            continue
        end

        locus, biotype, tlen = metadata[transcript_id]
        print(output, locus, "\t", transcript_id, "\t", biotype, "\t", tlen, "\t")

        x = ex[:tpm][i] / 1e6

        # absolute expression, which is just relative expression times the total
        # number of transcripts in the cell
        xabs = iround(transcript_count * x)

        if xabs <= 0
            println(output, "0.0\t0\t0.0\t0")
        else
            println(output, x, "\t", xabs)
        end
    end
end


main()

