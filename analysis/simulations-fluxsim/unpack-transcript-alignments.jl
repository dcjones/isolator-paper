#!/usr/bin/env julia

# A simple program to turn transcript alignments in the SAM format into
# genome alignments.

immutable Interval
    start::Int
    stop::Int
end


function Base.isless(a::Interval, b::Interval)
    return (a.start < b.start) || (a.start == b.start && a.stop < b.stop)
end

type Transcript
    seqname::String
    strand::Char
    exons::Vector{Interval}

    function Transcript(seqname::String, strand::Char)
        return new(seqname, strand, Interval[])
    end
end


# Crude parsing of GTF files
function parse_gtf(input::IO)

    transcript_id_pat = r"transcript_id\s+\"?([^\"]+)"
    transcripts = (String => Transcript)[]

    #println(STDERR, "parsing GTF...")
    for (i, line) in enumerate(eachline(input))
        #if i % 100000 == 0
            #println(STDERR, "    ", i, " lines")
        #end

        if line[1] == '#'
            continue
        end

        row = split(strip(line), '\t')
        if row[3] != "exon"
            continue
        end

        m = match(transcript_id_pat, row[9])
        if m === nothing
            continue
        end

        transcript_id = m.captures[1]
        if !haskey(transcripts, transcript_id)
            seqname = row[1]
            strand = row[7][1]
            t = Transcript(seqname, strand)
            transcripts[transcript_id] = t
        else
            t = transcripts[transcript_id]
        end

        exonstart = parseint(row[4])
        exonend = parseint(row[5])
        push!(t.exons, Interval(exonstart, exonend))
    end

    for t in values(transcripts)
        sort!(t.exons)
    end

    return transcripts
end


function print_nonneg_decimal_int(io::IO, x::Int)
    n = ndigits(x)
    u = 10^(n - 1)
    for i in 1:n
        d, x = divrem(x, u)
        u = div(u, 10)
        write(io, '0' + d)
    end
end


const cigar_buf = IOBuffer()

function adjust_position_cigar(t::Transcript, readlen::Int, tpos::Int)
    off = 0
    gpos = 0
    for (i, exon) in enumerate(t.exons)
        off += exon.stop - exon.start + 1
        if tpos < off
            if gpos == 0
                gpos = exon.stop - (off - tpos)
            end

            oplen = min(readlen, off - tpos)
            print_nonneg_decimal_int(cigar_buf, oplen)
            print(cigar_buf, 'M')
            readlen -= oplen
            tpos += oplen

            if tpos > 0 && readlen > 0
                # we are at the last exon, no splicing allowed
                @assert i < length(t.exons)
                print_nonneg_decimal_int(cigar_buf, t.exons[i + 1].start - exon.stop - 1)
                print(cigar_buf, 'N')
            end
        end

        if readlen == 0
            break
        end
    end

    return (gpos + 1, takebuf_string(cigar_buf))
end


function adjust_position(t::Transcript, tpos::Int)
    off = 0
    gpos = 0
    for (i, exon) in enumerate(t.exons)
        off += exon.stop - exon.start + 1
        if tpos < off
            gpos = exon.stop - (off - tpos)
            break
        end
    end

    return gpos + 1
end


function process_transcript_alignments(input::IO, output::IO, transcripts)
    #println(STDERR, "processing sam...")

    tab_positions = Array(Int, 10)

    for (i, line) in enumerate(eachline(input))
        #if i % 100000 == 0
            #println(STDERR, "    ", i, " lines")
        #end

        if line[1] == '@'
            continue
        end

        j = 1
        i = 1
        while j <= 10 && i < length(line)
            p = search(line, '\t', i)
            if p == 0
                error("Premature end of line")
            end
            tab_positions[j] = p
            i = p + 1
            j += 1
        end

        if line[tab_positions[2]+1] == '*'
            continue
        end

        l = length(line[tab_positions[9]+1:tab_positions[10]-1]) # read length
        t = transcripts[line[tab_positions[2]+1:tab_positions[3]-1]] # transcript
        tpos = parseint(line[tab_positions[3]+1:tab_positions[4]-1]) - 1 # transcript-coordinate position

        flags = parseint(line[tab_positions[1]+1:tab_positions[2]-1])
        strand = flags & 16 == 0 ? '+' : '-'
        if strand == '+'
            strand = t.strand
        elseif t.strand == '+'
            strand = '-'
        elseif t.strand == '-'
            strand = '+'
        end

        if strand == '+'
            flags = (flags | 16) $ 16
        else
            flags |= 16
        end

        gpos, cigar = adjust_position_cigar(t, l, tpos)
        mgpos       = adjust_position(t,
            parseint(line[tab_positions[8]+1:tab_positions[9]-1]))

        print(output, line[1:tab_positions[1]])
        print_nonneg_decimal_int(output, int(flags))
        print(output,
              '\t', t.seqname, '\t')
        print_nonneg_decimal_int(output, int(gpos))
        print(output,
              line[tab_positions[4]:tab_positions[5]],
              cigar,
              line[tab_positions[6]:tab_positions[7]])
        print_nonneg_decimal_int(output, int(mgpos))
        print(output, line[tab_positions[8]:end])
    end
end


if length(ARGS) < 3
    println(STDERR, "Usage: unpack-transcript-alignments.jl genes.gtf in.sam out.sam")
    exit(1)
end

transcripts = parse_gtf(open(ARGS[1]))
process_transcript_alignments(open(ARGS[2]), open(ARGS[3], "w"), transcripts)


