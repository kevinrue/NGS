library(ggbio)

library(GenomicRanges)

df <- data.frame(
    seqnames = c(rep("chr1", 1200), rep("chr2", 1000), rep("chr3", 800)),
    start = c(1:1200, 1:1000, 1:800),
    distance = abs(rnorm(3000)) * 1E3
)

gr <- GRanges(
    seqnames = c(rep("chr1", 1200), rep("chr2", 1000), rep("chr3", 800)),
    ranges = IRanges(start = c(1:1200, 1:1000, 1:800), width = 5),
    distance = rnorm(3000)
)

plotGrandLinear(gr, aes(y = distance))
