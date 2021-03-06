using GZip

import Base: eof, close, position
export BamReader, close, value, eof, advance!, eachposition

type BamReader
    bamStream
    done::Bool
    tid_to_refname::Dict{Int32,ASCIIString}
end

type Cigar
    length::UInt32
    operation::UInt8
end

type Read
    data::Array{Int8}
end

function BamReader(bamFileName::UTF8String)
    BamReader(ascii(bamFileName))
end

function BamReader(bamFileName::ASCIIString)
    f = GZip.open(bamFileName)

    # make sure this is a BAM file
    code = read(f, UInt8, 4)
    @assert code == b"BAM\1"

    # get through the header data
    l_text = read(f, Int32)
    skip(f, l_text)

    # make sure the contigs match our reference
    n_ref = read(f, Int32)

    tid_to_refname = Dict{Int32,ASCIIString}()

    for j in 1:n_ref
        l_name = read(f, Int32)
        refName = convert(ASCIIString, read(f, UInt8, l_name)[1:end-1]) # ignore the null terminator
        l_ref = read(f, Int32)
        tid_to_refname[j] = refName
    end

    BamReader(f, false, tid_to_refname)
end

close(read::Read) = GZip.close(read.bamStream)
refid(read::Read) = reinterpret(Int32,read.data[1:4])[1] + 1
position(read::Read) = reinterpret(Int32,read.data[5:8])[1] + 1
l_read_name(read::Read) = read.data[9]
mq(read::Read) = read.data[10]
bin(read::Read) = reinterpret(Int16,read.data[11:12])[1]
n_cigar_op(read::Read) = reinterpret(Int16,read.data[13:14])[1]
flag(read::Read) = reinterpret(Int16,read.data[15:16])[1]
l_seq(read::Read) = reinterpret(Int32,read.data[17:20])[1]
next_refid(read::Read) = reinterpret(Int32,read.data[21:24])[1] + 1
next_position(read::Read) = reinterpret(Int32,read.data[25:28])[1] + 1
tlen(read::Read) = reinterpret(Int32,read.data[29:32])[1] + 1
read_name(read::Read) = ascii(reinterpret(UInt8, read.data[33:31+l_read_name(read)]))
is_paired(read::Read) = (flag(read::Read) & 1) > 0
is_proper_paired(read::Read) = (flag(read::Read) & 2) > 0
is_unmapped(read::Read) = (flag(read::Read) & 4) > 0
mate_unmapped(read::Read) = (flag(read::Read) & 8) > 0
is_reverse(read::Read) = (flag(read::Read) & 16) > 0
mate_reverse(read::Read) = (flag(read::Read) & 32) > 0
is_read1(read::Read) = (flag(read::Read) & 64) > 0
is_read2(read::Read) = (flag(read::Read) & 128) > 0
not_primary(read::Read) = (flag(read::Read) & 256) > 0
fail_qc(read::Read) = (flag(read::Read) & 512) > 0
is_duplicate(read::Read) = (flag(read::Read) & 1024) > 0
supplementary_alignment(read::Read) = (flag(read::Read) & 2048) > 0

function cigar(read::Read)
    start_pos = 33 + l_read_name(read)
    stop_pos = start_pos + n_cigar_op(read) * 4
    cigar_array = reinterpret(UInt32,read.data[start_pos:stop_pos])
    [(c & 15, c >> 4) for c in cigar_array]
end

function seq(read::Read, trans = b"=ACMGRSVTWYHKDBN")
    seq_length = l_seq(read)
    start_pos = 33 + l_read_name(read) + n_cigar_op(read) * 4
    stop_pos = start_pos + div(seq_length + 1, 2)
    seq = reinterpret(UInt8,read.data[start_pos:stop_pos])

    ret = Array{UInt8}(seq_length)
    for i = 1:seq_length
        s = seq[div(i + 1,2)]
        shift = (i % 2) * 4
        mask = 15 << shift
        c = ((s & mask) >> shift) + 1
        ret[i] = trans[c]
    end

    ascii(ret)
end

function qualities(read::Read)
    seq_length = l_seq(read)
    start_pos = 33 + l_read_name(read) + n_cigar_op(read) * 4 + div(seq_length + 1, 2)
    stop_pos = start_pos + seq_length - 1
    a = reinterpret(UInt8,read.data[start_pos:stop_pos])
    for i = 1:length(a)
        a[i] = a[i] + 33
    end

    ascii(a)
end

# reader functions
reference_name(reader::BamReader, tid::Int64) = reader.tid_to_refname[tid]
eof(reader::BamReader) = reader.position == -1

function nextread(r::BamReader)
    f = r.bamStream
    if peek(f) == -1 # eof does not work on the BAM files either in C++ or here (BAM vs. gzip issue?)
        r.done = true
        r.position = -1
        return
    end

    buf = Array(Int32, 1)
    gzread(f, pointer(buf), 4)
    block_size = buf[1]

    read = Array(Int8, block_size) #the raw bytes
    gzread(f, pointer(read), block_size)

    Read(read)
end

## here we want to update the reader
#eachposition(r::BamReader) = BamReaderIterator(r)
#immutable BamReaderIterator
#	reader::BamReader
#end
#Base.start(it::BamReaderIterator) = it.reader.position
#Base.done(it::BamReaderIterator, position) = position == -1
#function Base.next(it::BamReaderIterator, position)
#	pos = it.reader.position
#	advance!(it.reader)
#	pos,it.reader.position
#end
