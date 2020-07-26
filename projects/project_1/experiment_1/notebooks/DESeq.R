######################################################################## 
rm(list=ls())
gc()
setwd("~/projects/tmp/code/YC_combined")
source("../../../fuc/myfunc.R")
require(dplyr)
require(openxlsx)

setting = 1

for(setting in c(1)){
  # read settings
  source("settings.R") 
  
  # reload data if possible
  pA.df = read.csv(file.path(result_dir, "pA.df.csv"), as.is = T)
  # Extend 3'UTRs
  if(!exists("threeUTR_extended") && extend_threeUTR){
    source("../../../fuc/extend3UTRs.R")
    threeUTR_extended = T
  }
  
  # Map clusters to genomic features
  # using the following line in clusters2features_hg19.R
  # cluster = subset(cluster, rowSums(cluster[,-c(1:3)] >= min_cluster_read_num) >= 2)
  allow_multiple_regions = T 
  source("../../../fuc/clusters2features_hg19.R")
  # match with enhancer RNA
  eRNA = openxlsx::read.xlsx(file.path(result_dir, "result_500_encode_middle2.xlsx"))
  eRNA.gr = GRanges(seqnames=eRNA$chr,
                    IRanges(start=eRNA$start, end=eRNA$end),
                    strand="*",
                    eRNA_id = eRNA$eRNA_id)
  olp = findOverlaps(pA, eRNA.gr)
  pA$eRNA_id = ""
  pA[queryHits(olp)]$eRNA_id = eRNA.gr[subjectHits(olp)]$eRNA_id
  
  eRNA.pA = pA[queryHits(olp)]
  eRNA.pA$region
  
  eRNA.pA = as.data.frame(eRNA.pA)
  eRNA.pA$conserveLV = NULL
  eRNA.pA$end = NULL
  eRNA.pA$width = NULL
  eRNA.pA$signed_position = NULL
  eRNA.pA$seqnames = NULL
  names(eRNA.pA) = sub("^start$", "pA_position", names(eRNA.pA))
  
  combined = merge(eRNA, eRNA.pA, by="eRNA_id", all.x = T, sort = T)
  names(combined) = sub(".+\\.y", "pA_position", names(combined))
  names(combined) = sub("(.+)\\.x", "\\1", names(combined))
  openxlsx::write.xlsx(combined, file.path(result_dir, "updated4_result_500_encode_middle2.xlsx"))
  
  
  # Filter low expression pAs 
  source("../../../fuc/filter_pAs3.R")
  
  # Generate 3'READS report (Wencheng's style) 
  source("../../../fuc/3READS_Report_with_RPM.R")
  tmp = pA.df[rowSums(pA.df[,grep("count$", names(pA.df), value=T)] > 5) == length(grep("count$", names(pA.df))),]
  m = tmp[, grep("_rpm_ratio$", names(tmp))]
  require(pheatmap)
  pheatmap.2(m, scale = "column", clustering_distance_rows = "correlation")
  
  # Calculate RED
  source("../../../fuc/threeUTR_RED4.R")
  
  ##################### Differential gene expressiong using DESeq2
  names(pA.df)
  require(DESeq2)
  require(GenomicAlignments	)
  ?summarizeOverlaps
  ?SummarizedExperiment
  
  ## 4.2 Starting from count matrices
  # use read counts in 3'UTR to represent gene expression
  genes = pA.df[, c("gene_symbol", "gene_id", grep("[12]_3utr$", names(pA.df),value=T))]
  genes = unique(genes)
  genes = na.omit(genes)
  genes = genes[rowSums(genes[, grep("3utr", names(genes))]) > 1, ]
  counts <- as.matrix(genes[, grep("_3utr$", names(genes),value=T)])
  colnames(counts) = sub("_3utr$", "", colnames(counts))
  rownames(counts) = genes$gene_symbol
  
  # rowRanges <- GRanges(seqnames = pA.df$chr,
  #                      ranges = IRanges(start=pA.df$pA_pos, width=1),
  #                      strand=pA.df$strand,
  #                      feature_id=pA.df$pAid,
  #                      gene_symbol=pA.df$gene_symbol,
  #                      gene_id=pA.df$gene_id,
  #                      region=pA.df$region,
  #                      UTR_length=pA.df$UTR_length)
  # rowRanges <- GRanges(seqnames = genes$gene_symbol,
  #                      ranges = IRanges(start=rep(1, nrow(counts)), width=rep(1, nrow(counts))),
  #                      strand=rep("+", nrow(counts)),
  #                      #feature_id=genes$gene_symbol,
  #                      gene_symbol=genes$gene_symbol,
  #                      gene_id=genes$gene_id
  #                      )
  # 
  
  colData <- DataFrame(Treatment=factor(sub("(.+)_(.+)_(.+)", "\\1", colnames(counts)), levels=c("siCtrl", "siINTS4",  "siZFC3H1", "siZCCHC8")),
                       Batch= sub("(.+)_(.+)_(.+)", "\\3", colnames(counts)),
                       row.names=colnames(counts))

  dds <- DESeqDataSetFromMatrix(countData = counts,
                                   colData = colData,
                                   design = ~ Treatment)
  
  head(assay(dds))
  colData(dds)
  rowData(dds)
  
  nrow(dds)
  ## 5.1 Pre-filtering the dataset
  # Additional weighting/filtering to improve power is applied at a later step in the workflow.
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  nrow(dds)
  
  ## 5.2 The rlog and variance stabilizing transformations
  # Many common statistical methods for exploratory analysis of multidimensional data, for example clustering and principal components 
  # analysis (PCA), work best for data that generally has the same range of variance at different ranges of the mean values. When the 
  # expected amount of variance is approximately the same across different mean values, the data is said to be homoskedastic. For 
  # RNA-seq counts, however, the expected variance grows with the mean. For example, if one performs PCA directly on a matrix of counts 
  # or normalized counts (e.g. correcting for differences in sequencing depth), the resulting plot typically depends mostly on the genes 
  # with highest counts because they show the largest absolute differences between samples. A simple and often used strategy to avoid 
  # this is to take the logarithm of the normalized count values plus a pseudocount of 1; however, depending on the choice of 
  # pseudocount, now the genes with the very lowest counts will contribute a great deal of noise to the resulting plot, because 
  # taking the logarithm of small counts actually inflates their variance. We can quickly show this property of counts with some 
  # simulated data (here, Poisson counts with a range of lambda from 0.1 to 100). We plot the standard deviation of each row (genes) 
  # against the mean:
  
  lambda <- 10^seq(from = -1, to = 2, length = 1000)
  cts <- matrix(rpois(1000*100, lambda), ncol = 100)
  library("vsn")
  meanSdPlot(cts, ranks = FALSE)
  
  log.cts.one <- log2(cts + 1)
  meanSdPlot(log.cts.one, ranks = FALSE)
  
  #As a solution, DESeq2 offers two transformations for count data that stabilize the variance across the mean: the regularized-logarithm transformation or rlog (Love, Huber, and Anders 2014), and the variance stabilizing transformation (VST) for negative binomial data with a dispersion-mean trend (Anders and Huber 2010), implemented in the vst function.
  
  #For genes with high counts, the rlog and VST will give similar result to the ordinary log2 transformation of normalized counts. For genes with lower counts, however, the values are shrunken towards the genes’ averages across all samples. The rlog-transformed or VST data then becomes approximately homoskedastic, and can be used directly for computing distances between samples, making PCA plots, or as input to downstream methods which perform best with homoskedastic data.
  # The rlog tends to work well on small datasets (n < 30), sometimes outperforming the VST when there is a large range of sequencing depth across samples (an order of magnitude difference).
  # Note that the two transformations offered by DESeq2 are provided for applications other than differential testing. For differential testing we recommend the DESeq function applied to raw counts, as described later in this workflow, which also takes into account the dependence of the variance of counts on the mean value during the dispersion estimation step.
  rld <- rlog(dds, blind = FALSE)
  head(assay(rld), 3)
  vsd <- vst(dds, blind = FALSE)
  head(assay(vsd), 3)
  
  # the effect of the transformation
  library("dplyr")
  library("ggplot2")
  
  dds <- estimateSizeFactors(dds)
  df <- bind_rows(
    as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"),
    as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"),
    as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  colnames(df)[1:2] <- c("x", "y")
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation) 
  
  ## 5.3 Sample distances
  # To ensure we have a roughly equal contribution from all genes, we use it on the rlog-transformed data.
  sampleDists <- dist(t(assay(rld)))
  sampleDists
  library("pheatmap")
  library("RColorBrewer")
  sampleDistMatrix <- as.matrix( sampleDists )
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  
  # 5.4 PCA plot
  png(file.path(result_dir, "sample_distances_using_rlog_transformed_gene_expression.png"), 800, 800)
  plotPCA(rld, intgroup = "Treatment")
  dev.off()
  
  # 5.5 MDS plot
  # This is useful when we don’t have a matrix of data, but only a matrix of distances. 
  mds <- as.data.frame(colData(rld))  %>%
    cbind(cmdscale(sampleDistMatrix))
  ggplot(mds, aes(x = `1`, y = `2`, color = Treatment)) +
    geom_point(size = 3) + coord_fixed()
  
  ## 7.3 Gene clustering
  library("genefilter")
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 100)
  mat  <- assay(rld)[ topVarGenes, ]
  mat  <- mat - rowMeans(mat)
  #anno <- as.data.frame(colData(rld)[, c("Treatment")])
  pheatmap(mat, filename = file.path(result_dir, "top_100_DE_genes_across_samples.png"), fontsize = 8, width=8, height=11)
  
  ## 7.4 Independent filtering
  qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
  bins <- cut(resLFC1$baseMean, qs)
  levels(bins) <- paste0("~", round(signif((qs[-1] + qs[-length(qs)])/2, 2)))
  fractionSig <- tapply(resLFC1$pvalue, bins, function(p)
    mean(p < .05, na.rm = TRUE))
  barplot(fractionSig, xlab = "mean normalized count",
          ylab = "fraction of small p values")
  # The DESeq2 software automatically performs independent filtering that maximizes the number of genes with adjusted p value less than a critical value (by default, alpha is set to 0.1). This automatic independent filtering is performed by, and can be controlled by, the results function.
  # Such filtering is permissible only if the statistic that we filter on (here the mean of normalized counts across all samples) is independent of the actual test statistic (the p value) under the null hypothesis. 
  
  
  ## Running the differential expression pipeline
  dds <- DESeq(dds)
  ############## Other comparisons & reporting
  for(treatment in c("siINTS4", "siZCCHC8", "siZFC3H1")){
    ## 6.3 comparisons
    # the name of the variable, the name of the level for the numerator, and the name of the level for the denominator.
    #treatment = "siZFC3H1"
    contrast = c("Treatment", treatment, "siCtrl")
    res = results(dds, contrast = contrast)
    res[order(res$padj),]
    summary(res)
    
    resSig <- subset(res, padj < 0.1)
    head(resSig[ order(resSig$log2FoldChange), ]) # downregulated
    head(resSig[ order(resSig$log2FoldChange, decreasing=T), ]) # upregulated
    
    ## 7.2 MA-plot
    es <- lfcShrink(dds, contrast = contrast, res=res)
    png(file.path(result_dir, paste0(treatment, "_vs_", "siCtrl_MAplot.png")))
    plotMA(res, ylim = c(-5, 5), main = paste0(treatment, " vs ", "siCtrl"))
    dev.off()
    hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
         col = "grey50", border = "white")
    
    ## 8 Annotating and exporting results
    library("AnnotationDbi")
    library("org.Hs.eg.db")
    columns(org.Hs.eg.db)
    res$entrez_id <- mapIds(org.Hs.eg.db,
                            keys=row.names(res),
                            column="ENTREZID",
                            keytype="SYMBOL",
                            multiVals="first")
    res$gene_name <- mapIds(org.Hs.eg.db,
                            keys=row.names(res),
                            column="GENENAME",
                            keytype="SYMBOL",
                            multiVals="first")
    
    resOrdered <- res[order(res$pvalue),]
    head(resOrdered)
    resOrderedDF <- as.data.frame(resOrdered)
    write.csv(resOrderedDF, file = file.path(result_dir, paste0(treatment, "_vs_", "siCtrl_gene_FC_DESeq2.csv")))
  }
  
 
  
 
 
 
 
 
 
 
 
 ################################# APA analysis using DEXSeq
 # each row is like gene_id:exonic_part_number
 # gene_id will be aUTR id, defined by proximal and distal pA's ids.
 # exonice_part_number will be 1 or 2 for proximal and distal pAs, respectively
 
 library(DEXSeq)
 # use this function:
 ?DEXSeqDataSet
  names(pA.df) 
  

 n = 5
 APA.summary = list()
 for(treatment in c("siINTS4", "siZFC3H1", "siZCCHC8")){
   # treatment = "siINTS4"
   tmp = pA.df[, c(1:14, grep(paste0(paste0("(", treatment, "|siCtrl)"), ".+[12]_count$|pA_type|Num_pA|signed_pA_pos"), names(pA.df)))]
   tmp = tmp[rowSums(tmp[, grep(".+[12]_count$", names(tmp))] >= n) == length(grep(".+[12]_count$", names(tmp))), ]
   tmp = subset(tmp, pA_type != "Intergenic")
   # only keep genes with >1 pAs
   tmp = subset(tmp, Num_pA > 1)
   tmp = tmp[tmp$gene_symbol %in% tmp$gene_symbol[duplicated(tmp$gene_symbol)], ]
   str(tmp)
   # sort the data frame
   tmp = tmp[order(tmp$gene_symbol, tmp$signed_pA_pos),]
   
   # prepare input for ?DEXSeqDataSet
   counts = tmp[, grep("_count$", names(tmp))]
   groupID = tmp$gene_symbol
   featureID = with(tmp, paste0(chr, strand, pA_pos))  
   
   sampleData <- data.frame(Treatment=factor(sub("(.+)_(.+)_(.+)_(.+)", "\\1", colnames(counts)), levels=c("siCtrl",treatment)),
                            Batch= sub("(.+)_(.+)_(.+)_(.+)", "\\3", colnames(counts)),
                            row.names=colnames(counts))
   
   #sampleData = data.frame(condition = c("untreated", "treated", "untreated", "treated"))
   design <- formula( ~ sample + exon + Treatment:exon )
   
   dxd = DEXSeqDataSet( counts, sampleData, design, featureID, groupID )
   
   #geneIDs(dxd)
   colData(dxd)
   head( counts(dxd), 5 )
   split( seq_len(ncol(dxd)), colData(dxd)$exon )
   head( featureCounts(dxd), 5 )
   head( rowData(dxd), 3 )
   sampleAnnotation( dxd )
   dxd = estimateSizeFactors( dxd )
   dxd = estimateDispersions( dxd )
   plotDispEsts( dxd )
   dxd = testForDEU( dxd )
   dxd = estimateExonFoldChanges( dxd, fitExpToVar="Treatment")
   dxr1 = DEXSeqResults( dxd )
   table ( dxr1$padj < 0.05 )
   plotMA( dxr1, cex=0.8, alpha=0.05)
   
   dxr1 = subset(dxr1, padj < 0.05)
   #dxr1$gene_symbol = sub("(.+):(.+)", "\\1", rownames(dxr1))
   dxr1$pAid = sub("(.+):(.+)", "\\2", rownames(dxr1))
   dxr1 = split(as.data.frame(dxr1), dxr1$groupID)
   dxr1 = lapply(dxr1, function(gene) gene[which.min(gene$padj), ])
   dxr1 = do.call(rbind, dxr1)
   
   dxr1$pA_type = tmp$pA_type[match(dxr1$featureID, tmp$pAid)]
   dxr1$APA_direction = ifelse(dxr1[, grep("log2fold", names(dxr1))] > 0, "UP", "DWN")
   APA.summary[[paste0(treatment, "_siCtrl")]] = as.matrix(table(dxr1$APA_direction, dxr1$pA_type))
   
   write.csv(dxr1, file.path(result_dir, paste0("DEXSeq_APA_", treatment, "_siCtrl.csv")))
 }
 APA.summary[[3]] = cbind(E=c(0,0), APA.summary[[3]])
 cmp = rep(sub("_", "/", names(APA.summary)), each=2)
 
 APA.summary = do.call(rbind, APA.summary)
 direction = rownames(APA.summary)
 rownames(APA.summary) = NULL
 APA.summary = cbind(Sample_pair = cmp, Direction = direction, as.data.frame(APA.summary))
 write.csv(APA.summary, file.path(result_dir, "APA_count_summary_fdr_0.05.csv"), row.names = F)
 
  
  
  
  
  
 
  
}
