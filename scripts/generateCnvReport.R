library(GenomicRanges, quietly=T)
library(VariantAnnotation, quietly=T)
library(stringr, quietly=T)

# Author: Christopher Medway
# Use: RScript generateCnvReport.R <seqId> <panel>
# Description: compiles a single run-wide report of CNVs, merging CNV calls
# from MANTA and ExomeDepth

# input arguments
roi <- commandArgs(trailingOnly = T)[1] #The bed file containing the roi
exomeDepthMatrics <- commandArgs(trailingOnly = T)[2]#Path to the exome metrics file
high_cov_bams <- commandArgs(trailingOnly = T)[3] #Path to the high coverage bam file
depth_folder <- commandArgs(trailingOnly = T)[4] #Path to the directory containing the depth summary files
exome_depth_dir <- commandArgs(trailingOnly = T)[5] #The folder containing the exome_depth results
manta_dir <- commandArgs(trailingOnly = T)[6] #The folder containing the manta results
out_name <- commandArgs(trailingOnly = T)[7] #Where to put the output csv


# load exome depth QC (r2 values)
exomeDepthMatricsDf <- read.table(exomeDepthMatrics, stringsAsFactors = F, header = T)

print (high_cov_bams)


# get list of samples that have CNV data
sampleId <- unlist(lapply(stringr::str_split(read.table(high_cov_bams, stringsAsFactors=F)[,1], pattern="/"), function(x) x[3]))
print(sampleId)
sampleId <-  unlist(lapply(stringr::str_split(sampleId, pattern="_"), function(x) paste0(x[1],"_", x[2])))
print (sampleId)


# generate initial QC flags where all samples = PASS
passQC <- data.frame(sampleId, "QC" = "PASS", stringsAsFactors = F)

# flag any sample with R@ < 0.98 as a FAIL
if (sum(exomeDepthMatricsDf$Correlation < 0.98) > 0) {
    fails <- sub(exomeDepthMatricsDf[exomeDepthMatricsDf$Correlation < 0.98,"BamPath"], perl = T, pattern = ".+_(\\w+).bam", replacement = "\\1")
    passQC[passQC$sampleId %in% fails, "QC"] <- "R2<0.98"
}

# add average depth for each sample
for (s in sampleId) {
	

    exomeDepthCoverage <- paste0(depth_folder,s,"_DepthOfCoverage.sample_summary")
    exomeDepthCoverageDf <- read.table(exomeDepthCoverage, stringsAsFactors = F, header = T, fill = T)
    passQC[passQC$sampleId %in% s, "Depth"] <- exomeDepthCoverageDf$mean[1]
}

# main function - loop over sample and compile MANTA / ExoneDepth data
mergedDf <- lapply(sampleId, function(sample) {
    
    qc    <- passQC[passQC$sampleId %in% sample,"QC"]
    depth <- passQC[passQC$sampleId %in% sample, "Depth"] 
    
    exomeDepth <- paste0(exome_depth_dir, sample, "_final_cnv_fixed.vcf.gz")
    manta      <- paste0(manta_dir, sample,"_diploidSV.vcf.gz")
    
    ## extract information from exomeDepth
    vcf    <- VariantAnnotation::readVcf(exomeDepth, "hg19")
    
    if (dim(vcf)[1] > 0) {
        chr     <- as.vector(vcf@rowRanges@seqnames)
        start   <- as.numeric(vcf@rowRanges@ranges@start)
        end     <- as.numeric(vcf@info$END)
        ref     <- as.character(vcf@fixed$REF)
        alt     <- as.character(vcf@fixed$ALT@unlistData)
        type    <- vcf@fixed$ALT@unlistData
        regions <- vcf@info$Regions
        
        if (any(is.na(end))) {end[is.na(end)] <- (start[is.na(end)] + 1)}
        
    } else {
        
        chr     <- NA
        start   <- NA
        end     <- NA
        ref     <- NA
        alt     <- NA
        type    <- NA
        regions <- NA
    }
    
    exomeDepthDF <- data.frame(
        sample,
        "method" = "exomeDepth",
        chr,
        start,
        end,
        ref,
        type,
        regions,
        qc,
        depth,
        stringsAsFactors = F
    )
    
    rm(vcf, chr, start, end, ref, alt, regions, type)
    

    
    ## extract information from MANTA
    vcf     <- VariantAnnotation::readVcf(manta, "hg19")
    
    if (dim(vcf)[1] > 0) {
        chr     <- as.vector(vcf@rowRanges@seqnames)
        start   <- as.numeric(vcf@rowRanges@ranges@start)
        end     <- as.numeric(vcf@info$END)
        ref     <- as.character(vcf@fixed$REF)
        alt     <- as.character(vcf@fixed$ALT@unlistData)
        type    <- vcf@info$SVTYPE
        
        if (any(is.na(end))) {end[is.na(end)] <- start[is.na(end)] + 1}
        
    } else {
        
        chr   <- NA
        start <- NA
        end   <- NA
        ref   <- NA
        alt   <- NA
        type  <- NA
        }
    
    
    mantaDF <- data.frame(
        sample,
        "method" = "manta",
        chr,
        start,
        end,
        ref,
        type,
        "regions" = "-",
        qc = "PASS",
        depth,
        stringsAsFactors = F
    )
    
    rm(vcf, chr, start, end, ref, alt, type)
    
    df <- rbind(exomeDepthDF, mantaDF)
    dfNonNull <- df[!is.na(df$chr),]
    dfNull    <- df[is.na(df$chr),]
    
    ## annotate with gene names
    roiDf <- read.table(roi, stringsAsFactors = F)
    roiDf$gene <- stringr::str_extract(string = roiDf$V4, pattern="(\\w+)")
    roiUnified <- do.call(rbind, lapply(split(roiDf, roiDf$gene), function(x) {data.frame("chr" = x$V1[1], "start" = min(x$V2), "end" = max(x$V3), "gene" = x$gene[1])}))

    roiGR <- GenomicRanges::GRanges(seqnames = roiUnified$chr, ranges = IRanges(start = roiUnified$start, end = roiUnified$end))
    
    gr <- GenomicRanges::GRanges(seqnames = dfNonNull$chr, ranges = IRanges(start = as.numeric(dfNonNull$start), end = as.numeric(dfNonNull$end)))
    ol <- findOverlaps(gr, roiGR)
    dfNonNull[queryHits(ol),"gene"] <- roiUnified[subjectHits(ol),"gene"]
    dfNonNull <- dfNonNull[order(as.numeric(dfNonNull$chr), as.numeric(dfNonNull$start)),]
    
    if (dim(dfNull)[1] == 0) {
        return(dfNonNull)
    } else {
        dfNull$gene <- NA
        return(rbind(dfNonNull, dfNull))
    }
})

mergedDf <- do.call(rbind, mergedDf)

mergedDf[mergedDf$depth < 160 & mergedDf$qc == "R2<0.98","qc"] <- "R2<0.98;Depth<160"
mergedDf[mergedDf$depth < 160 & mergedDf$qc == "PASS","qc"]    <- "Depth<160"

# write to current directory (run folder)
write.table(mergedDf, file = paste0(out_name), row.names = F, sep = ",", quote = F)
