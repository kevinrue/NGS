
stopifnot(
    require('edgeR'),
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
    metaTags <- grepl('^__', rownames(RG[['counts']]))
    message('Meta tags dropped')
    # Exclude metatags if any; recalculate (true) library size
    RG <- RG[!metaTags, keep.lib.sizes = FALSE]
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
        dds <- DESeq2::DESeqDataSetFromMatrix(
            dge$counts, colData, formula(~ Group)
        )
    } else {
        message('NOTE: Batch information applied')
        colData <- dge$samples[,c('Group', 'Batch')]
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countMat, colData,formula(~ condition + Batch)
        )
    }
    dds <- DESeq2::DESeq(dds, betaPrior = F)
    res <- DESeq2::results(dds, contrast = c("Group", groupTarget, groupRef))
    # Write DE table to output file
    DEfile <- file.path(outFolder, 'DESeq2_table.csv')
    if (!dir.exists(outFolder)){stopifnot(dir.create(outFolder))}
    write.csv(res, DEfile)
    return(res)
}

# getDE.edgeR ----

getDE.edgeR <- function(dge, groupTarget, groupRef, outFolder = 'edgeR'){
    if (all(is.na(dge[['samples']][,'Batch']))){
        design <- with(
            dge[['samples']],
            model.matrix(~ Group)
        )
    } else {
        message('NOTE: Batch information applied')
        design <- with(
            dge[['samples']],
            model.matrix(~ Group + Batch)
        )
    }
    message('calcNormFactors ...')
    dge <- edgeR::calcNormFactors(dge)
    message('estimateDisp ...')
    dge <- edgeR::estimateDisp(dge, design)
    message('glmQLFit ...')
    fit <- edgeR::glmQLFit(dge, design)
    # Define the comparison to perform
    myContrast <- rep(0, ncol(design))
    idxGroupTarget <- (colnames(design) == sprintf("Group%s", groupTarget))
    idxGroupRef <- (colnames(design) == sprintf("Group%s", groupRef))
    myContrast[idxGroupTarget] <- 1
    myContrast[idxGroupRef] <- (-1)
    message('glmLRT ...')
    results <- glmLRT(fit, contrast = myContrast)
    # et <- edgeR::exactTest(deg, c(groupRef, groupTarget))
    # lrt <- glmQLFTest(fit, contrast = myContrast)
    # print(colnames(results))
    tp <- edgeR::topTags(results, n = nrow(dge), sort.by = "none")[[1]]
    # Write DE table to output file
    DEfile <- file.path(outFolder, 'edgeR_table.csv')
    if (!dir.exists(outFolder)){stopifnot(dir.create(outFolder))}
    write.csv(tp, DEfile)
    return(tp)
}

# getLogCPM.edgeR ----

getLogCPM.edgeR <- function(dge, outFolder = 'edgeR'){
    dge <- calcNormFactors(dge)
    logCPMmatrix <- edgeR::cpm(dge, log = TRUE)
    return(logCPMmatrix)
}
