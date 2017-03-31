
stopifnot(
    requireNamespace('edgeR'),
    requireNamespace('DESeq2')
)

# Parameters ----

significance <- 0.01

# prepareFileInfo ----

# Create a data.frame that associates RNAaseq count files to samples
## Assumptions:
## The given root folder immediately contains one subfolder by sample name
## Structure of Rafik's pipeline
# Produces a data.frame with two columns: "files", "sampleNames"
prepareFileInfo.pattern <- function(
    folder, pattern = "%s_Treat.nodup.filtered.hits"){
    # list all immediate subfolders (expected to be sample names)
    sampleNames <- list.files(folder)
    # Derive expected full file paths !!! Hard-coded !!!
    countFiles <- file.path(
        folder, sampleNames, 'htseq', sprintf(pattern, sampleNames)
    )
    # Produce an edgeR-friendly data.frame to import counts
    fileInfo <- data.frame(
        files = countFiles,
        SampleName = sampleNames,
        row.names = sampleNames
    )
    return(fileInfo)
}

# Create a data.frame that associates RNAaseq count files to samples
## Assumptions:
## The given root folder immediately contains one subfolder by sample name
## Structure of Rafik's pipeline
# Produces a data.frame with two columns: "files", "sampleNames"
prepareFileInfo.fixed <- function(folder, fixed = 'treat_name.hits'){
    # list all immediate subfolders (expected to be sample names)
    sampleNames <- list.files(folder)
    # Derive expected full file paths !!! Hard-coded !!!
    countFiles <- file.path(folder, sampleNames, fixed)
    # Produce an edgeR-friendly data.frame to import counts
    fileInfo <- data.frame(
        files = countFiles,
        SampleName = sampleNames,
        Group = NA, Batch = NA,
        row.names = sampleNames
    )
    return(fileInfo)
}

# setGroupInfo ----

# Set the Group (i.e. Condition) information for each or all samples
# Return the input data.frame with Group column set to given value(s)
setGroupInfo <- function(x, condition){
    x[,'Group'] <- condition
    return(x)
}

# setBatchInfo ----

# Set the Batch information for each or all samples
# Return the input data.frame with Batch column set to given value(s)
setBatchInfo <- function(x, condition){
    x[,'Batch'] <- condition
    return(x)
}

# importCounts ----

# Read and Merge a Set of Files Containing Count Data
# Given a data.frame that links files and phenotype information
# Exclude metatags (e.g., count of unassigned reads)
# Return an object that stores both counts and samples information
importCounts <- function(x){
    # Parse all counts files
    RG <- edgeR::readDGE(x, group = x$Group, labels = x$SampleName)
    # Identify metatags if any
    metaTags <- grepl('^__', rownames(RG$counts))
    message('Meta tags dropped')
    # Exclude metatags if any
    RG <- RG[!metaTags,]
    # Recalculate the (true) library size after removing meta-tags
    RG$samples$lib.size <- colSums(RG$counts)
    return(RG)
}

# filterCounts ----

# Keep features (i.e., genes) expressed above a cutoff in a minimum fraction of samples
# Return the input object subsetted to the retained features
filterCounts <- function(x, cpm = 1, fraction = 0.8){
    cpmMatrix <- edgeR::cpm(x, normalized.lib.sizes = FALSE)
    keep <- (rowSums(cpmMatrix > cpm) > (ncol(x) * fraction))
    xSubset <- x[keep,]
    return(xSubset)
}

# getDE.DESeq2 ----

getDE.DESeq2 <- function(dge, groupTarget, groupRef, outFolder = 'DESeq2'){
    if (all(is.na(dge$samples$Batch))){
        colData <- dge$samples[,'Group', drop = FALSE]
        dds <- DESeqDataSetFromMatrix(
            dge$counts, colData, formula(~ Group)
        )
    } else {
        message('NOTE: Batch information supplied')
        colData <- dge$samples[,c('Group', 'Batch')]
        dds <- DESeqDataSetFromMatrix(
            countMat, colData,formula(~ condition + Batch)
        )
    }
    dds <- DESeq(dds, betaPrior = F)
    res <- results(dds, contrast = c("Group", groupTarget, groupRef))
    return(res)
}

# getDE.edgeR ----

getDE.edgeR <- function(dge, groupTarget, groupRef, outFolder = 'edgeR'){
    if (all(is.na(dge$samples$Batch))){
        design <- with(
            dge$samples,
            model.matrix(~ Group)
        )
    } else {
        message('NOTE: Batch information supplied')
        design <- with(
            dge$samples,
            model.matrix(~ Group + Batch)
        )
    }
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design)
    # et <- exactTest(deg)
    fit <- glmQLFit(dge, design)
    # Define the comparison to perform
    myContrast <- rep(0, ncol(design))
    idxGroupTarget <- (names(myContrast) == sprintf("Group%s", groupTarget))
    idxGroupRef <- (names(myContrast) == sprintf("Group%s", groupRef))
    myContrast[idxGroupTarget] <- 1
    myContrast[idxGroupRef] <- -1
    lrt <- glmQLFTest(fit, contrast = myContrast)
    print(names(design))
    print(names(design))
    et_tags = topTags(lrt, n = nrow(dge), sort.by = "none")[[1]]
}

# getLogCPM.edgeR ----

getLogCPM.edgeR <- function(dge, outFolder = 'edgeR'){
    dge <- calcNormFactors(dge)
    logCPMmatrix <- edgeR::cpm(dge, log = TRUE)
    return(logCPMmatrix)
}
