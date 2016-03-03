module BamIO
export refid, position, l_read_name, mq, n_cigar_op, flag, l_seq, next_refid, next_position, tlen, read_name, cigar, seq, qualities, eof, nextread, reference_name
export is_paired, is_proper_paired, is_unmapped, mate_unmapped, is_reverse, mate_reverse, is_read1, is_read2, not_primary, fail_qc, is_duplicate, supplementary_alignment

include("bam.jl")

# package code goes here

end # module
