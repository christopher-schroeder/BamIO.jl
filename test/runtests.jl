using BamIO
using FactCheck

reader = BamReader("data/small.bam")

read = nextread(reader)

facts("read1") do
    @fact read_name(read) --> "FC30W98HM_20090202:4:76:1751:1274"
    @fact flag(read) --> 16
    @fact refid(read) --> 1
    @fact reference_name(reader, refid(read)) --> "chr1"
    @fact position(read) --> 10056
    @fact mq(read) --> 1
    @fact cigar(read) --> UInt32[0x000001b0]
    @fact next_refid(read) --> 0
    @fact BamIO.bin(read) --> 4681
    @fact next_position(read) --> 0
    @fact seq(read) --> "AACCCTAACCCTAACCCTAACCCTAAC"
    @fact qualities(read) --> "[[OPPV[[TVV[[[TTP[[[[[[[[[["
    @fact l_read_name(read) --> 34
    @fact n_cigar_op(read) --> 1
    @fact l_seq(read) --> 27
end

read = nextread(reader)

facts("read2") do
    @fact read_name(read) --> "FC30W98HM_20090202:4:86:34:1700"
    @fact flag(read) --> 16
    @fact refid(read) --> 1
    @fact reference_name(reader, refid(read)) --> "chr1"
    @fact position(read) --> 10102
    @fact mq(read) --> 1
    @fact cigar(read) --> UInt32[0x000001b0]
    @fact next_refid(read) --> 0
    @fact BamIO.bin(read) --> 4681
    @fact next_position(read) --> 0
    @fact seq(read) --> "CTAACCCAACCCTAACCCTAACCCTAA"
    @fact qualities(read) --> "UY[[YYY[[U[O[[[[[[[[[[[[[[["
    @fact l_read_name(read) --> 32
    @fact n_cigar_op(read) --> 1
    @fact l_seq(read) --> 27
end
