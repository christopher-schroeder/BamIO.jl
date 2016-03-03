using GZip

import Base: eof, close, position
export BamReader, close, value, eof, advance!, eachposition

type BamReader
    bamStream
    done::Bool
    position::Int32
    read::Array{Int8}
end

type Cigar
    length::UInt32
    operation::UInt8
end

type Read
    data::Array{Int8}
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

    for j in 1:n_ref
        l_name = read(f, Int32)
        refName = convert(ASCIIString, read(f, UInt8, l_name)[1:end-1]) # ignore the null terminator
        l_ref = read(f, Int32)
    end

    BamReader(f, false, 0, Array{Int8}(0))
end
close(reader::Read) = GZip.close(reader.bamStream)
refid(reader::Read) = reinterpret(Int32,reader.data[1:4])[1] + 1
position(reader::Read) = reinterpret(Int32,reader.data[5:8])[1] + 1
l_read_name(reader::Read) = reader.data[9]
mq(reader::Read) = reader.data[10]
bin(reader::Read) = reinterpret(Int16,reader.data[11:12])[1]
n_cigar_op(reader::Read) = reinterpret(Int16,reader.data[13:14])[1]
flag(reader::Read) = reinterpret(Int16,reader.data[15:16])[1]
l_seq(reader::Read) = reinterpret(Int32,reader.data[17:20])[1]
next_refid(reader::Read) = reinterpret(Int32,reader.data[21:24])[1] + 1
next_position(reader::Read) = reinterpret(Int32,reader.data[25:28])[1] + 1
tlen(reader::Read) = reinterpret(Int32,reader.data[29:32])[1] + 1
read_name(reader::Read) = ascii(reinterpret(UInt8, reader.data[33:31+l_read_name(reader)]))

function cigar(reader::Read)
    start_pos = 33 + l_read_name(reader)
    stop_pos = start_pos + n_cigar_op(reader) * 4
    reinterpret(UInt32,reader.data[start_pos:stop_pos])
end

function seq(reader::Read, trans = b"=ACMGRSVTWYHKDBN")
    seq_length = l_seq(reader)
    start_pos = 33 + l_read_name(reader) + n_cigar_op(reader) * 4
    stop_pos = start_pos + div(seq_length + 1, 2)
    seq = reinterpret(UInt8,reader.data[start_pos:stop_pos])

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

function qualities(reader::Read)
    seq_length = l_seq(reader)
    start_pos = 33 + l_read_name(reader) + n_cigar_op(reader) * 4 + div(seq_length + 1, 2)
    stop_pos = start_pos + seq_length - 1
    a = reinterpret(UInt8,reader.data[start_pos:stop_pos])
    for i = 1:length(a)
        a[i] = a[i] + 33
    end

    ascii(a)
end

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
