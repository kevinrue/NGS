snp <- read.table(system.file("extdata", "plink.assoc.sub.txt", package = "biovizBase"),
                  header = TRUE)

gr.snp <- transformDfToGr(snp, seqnames = "CHR", start = "BP", width = 1)

require(biovizBase)

gr.snp <- transformDfToGr(snp, seqnames = "CHR", start = "BP", width = 1)


gr.kevin <- GRanges(
    seqnames = snp$CHR,
    ranges = IRanges(start = snp$BP, width = 1),
    OddRatio = snp$OR)

plotGrandLinear(gr.kevin, aes(y = -OddRatio, color = seqnames))

# autoplot(gr.kevin, type = "a")
