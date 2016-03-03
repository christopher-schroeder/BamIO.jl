# BamIO

A simple bam reader for julia

example

```julia
using BamIO

## BamReader

# test forward reads
reader = BamReader("small.bam")

for i = 1:10
    read = nextread(reader)
    print(read_name(read), "\t")
    print(flag(read), "\t")
    print(refid(read), "\t")
    print(position(read), "\t")
    print(mq(read), "\t")
    print(cigar(read), "\t")
    print(next_refid(read), "\t")
    print(next_position(read), "\t")
    print(seq(read), "\t")
    print(qualities(read), "\t")
    println()

    #print(BamIO.bin(read), "\t")
    #print(l_read_name(read), "\t")
    #print(n_cigar_op(read), "\t")
    #print(l_seq(read), "\t")
end
```
