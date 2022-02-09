################################################################################
## This script contains the code to produce Figure 4 and Extended Data Figure 4
## Hetzel et al. The distinct DNA methylome of acute lymphoblastic leukemia 2022
################################################################################

library(pheatmap)
library(RColorBrewer)
library(data.table)
library(vioplot)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(gridExtra)
library(WebGestaltR)
library(RVAideMemoire)

################################################################################
## Figure 4a
################################################################################

avg_cgi_for_pca <- read.table("Figure4a_4b_4c_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")

dataForPca <- t(avg_cgi_for_pca)
pca <- prcomp(dataForPca)

samplesControlTALL <- c("SJNORM016314_G1", "SJNORM016314_G2", "SJNORM016314_G3", "SJNORM016314_G4", "SJNORM016314_G5", "SJNORM016315_G1", "SJNORM016315_G2", "SJNORM016315_G3", "SJNORM016315_G4", "SJNORM016315_G5")
samplesTALL <- c("SJALL030551_D1", "SJALL030933_D1", "SJALL048427_D3", "SJALL048442_D2", "SJALL064568_D1", "SJALL064571_D1", "SJALL064572_D1", "SJALL064573_D1", "SJTALL004_D", "SJTALL010_D", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012402_D1", "SJTALL022646_D1", "SJTALL187_D", "SJALL030493_D1", "SJALL030684_D1", "SJALL048428_D4", "SJALL048440_D2", "SJALL048449_D2", "SJALL064567_D1", "SJALL064569_D1", "SJTALL012401_D1", "SJTALL012403_D1", "SJTALL015715_D1", "SJTALL021687_D1", "SJTALL021688_D1", "SJTALL021912_D1", "SJTALL022648_D1", "SJTALL022649_D1", "SJTALL170_D", "SJTALL182_D", "SJTALL186_D", "SJALL030916_D1", "SJALL059529_D2", "SJALL064564_D1", "SJALL064565_D1", "SJALL064566_D1", "SJTALL009_D", "SJTALL021685_D1", "SJTALL021694_D1", "SJTALL021695_D1", "SJTALL022645_D1", "SJTALL022647_D1", "SJTALL169_D", "SJTALL174_D")

dfPCA <- data.frame(pca$x, type = factor(ifelse(rownames(pca$x) %in% samplesControlTALL, "control T cells", "T-ALL"), levels = c("control T cells", "T-ALL")), median_variable_cgi = apply(avg_cgi_for_pca[,c(samplesControlTALL, samplesTALL)], 2, median), fraction_variable_cgi = colSums(ifelse(avg_cgi_for_pca[,c(samplesControlTALL, samplesTALL)] > 0.2, 1, 0)) / nrow(avg_cgi_for_pca))

percentage <- summary(pca)$importance["Proportion of Variance",]
percentage <- paste(colnames(dfPCA), " (", as.character(round(percentage * 100, 2)), "%", ")", sep="")

pdf("Figure4a.pdf", width = 6.5, height = 5)
ggplot(dfPCA,aes(x=PC1,y=PC2, col = median_variable_cgi, shape = type))+
       geom_point(size=2) +
       xlab(percentage[1]) + ylab(percentage[2]) +
       scale_color_gradientn(colours = c("#0072B2", "seagreen3", "#F0E442"), limits=c(0,1)) +
       theme_classic()
dev.off()

################################################################################
## Figure 4b
################################################################################

avg_cgi_for_pca <- read.table("Figure4a_4b_4c_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")

dataForPca <- t(avg_cgi_for_pca)
dataForPcaBinary <- apply(dataForPca, 2, function(x) ifelse(x > 0.2, 1, 0))
pcaBinary <- prcomp(dataForPcaBinary)

dfPCABinary <- data.frame(pcaBinary$x, type = factor(ifelse(rownames(pcaBinary$x) %in% samplesControlTALL, "control T cells", "T-ALL"), levels = c("control T cells", "T-ALL")), median_variable_cgi = apply(avg_cgi_for_pca[,c(samplesControlTALL, samplesTALL)], 2, median), fraction_variable_cgi = colSums(ifelse(avg_cgi_for_pca[,c(samplesControlTALL, samplesTALL)] > 0.2, 1, 0)) / nrow(avg_cgi_for_pca))

percentageBinary <- summary(pcaBinary)$importance["Proportion of Variance",]
percentageBinary <- paste(colnames(dfPCABinary), " (", as.character(round(percentageBinary * 100, 2)), "%", ")", sep="")

pdf("Figure4b.pdf", width = 6.5, height = 5)
ggplot(dfPCABinary,aes(x=PC1,y=PC2, col = fraction_variable_cgi, shape = type))+
      geom_point(size=2) +
      xlab(percentageBinary[1]) + ylab(percentageBinary[2]) +
      scale_color_gradientn(colours = c("yellow1", "#FF8C00", "red3"), limits=c(0,1)) +
      theme_classic()
dev.off()

################################################################################
## Figure 4c
################################################################################

avg_cgi_for_heatmap <- read.table("Figure4a_4b_4c_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
cgi_annotation_complete <- read.table("Figure4c_Extended_Data_Figure4b_source_data_annotation_cgi.tsv", header = TRUE, row.names = 1, sep = "\t")
cluster_tall_df <- read.table("Figure4c_source_data_annotation_samples.tsv", header = TRUE, row.names = 1, sep = "\t")

annotation_colors <- list(Type = c("black", brewer.pal(3, "Reds")), Sex = c("gold", "royalblue4"), Age_Group = c(brewer.pal(n = 3, name = "YlGn")))
names(annotation_colors[[1]]) <- c("control T cells", "T-ALL LM", "T-ALL IM", "T-ALL HM")
names(annotation_colors[[2]]) <- c("female", "male")
names(annotation_colors[[3]]) <- c("Pediatric", "AYA", "Adult")

col_fun_median = colorRamp2(c(0, 0.5, 1), c("#0072B2", "seagreen3", "#F0E442"))
annotation_colors$median_variable_cgi <- col_fun_median

pdf("Figure4c.pdf", height = 8, width = 14)
Heatmap(as.matrix(avg_cgi_for_heatmap), row_split = cgi_annotation_complete[rownames(avg_cgi_for_heatmap),"cluster"], use_raster = TRUE, raster_quality = 2, show_row_names = FALSE, top_annotation =  columnAnnotation(df = cluster_tall_df[colnames(avg_cgi_for_heatmap),c("Type", "Age_Group", "Sex", "median_variable_cgi")], col = annotation_colors), cluster_row_slices = FALSE)
dev.off()

################################################################################
## Figure 4d
################################################################################

cgi_combined_feature_pct <- read.table(file = "Figure4d_source_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cgi_combined_feature_pct$cluster <- factor(cgi_combined_feature_pct$cluster, levels = c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high", "all_covered"))

pdf("Figure4d.pdf", width = 9, height = 5)
ggplot(cgi_combined_feature_pct, aes(x = cluster, y = fraction, color = feature, group = group)) + geom_line(size = 1) + geom_point(size=2) + ylab("Fraction") + xlab("") + theme_classic() + coord_cartesian(ylim = c(0,1))
dev.off()

################################################################################
## Figure 4e
################################################################################

df_states_per_cluster_hsc <- read.table(file = "Figure4e_source_data_hsc.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df_states_per_cluster_dnd41 <- read.table(file = "Figure4e_source_data_dnd41.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
df_states_per_cluster_hsc$cluster <- factor(df_states_per_cluster_hsc$cluster, levels = c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high", "all_covered"))
df_states_per_cluster_dnd41$cluster <- factor(df_states_per_cluster_dnd41$cluster, levels = c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high", "all_covered"))
df_states_per_cluster_hsc$state <- factor(df_states_per_cluster_hsc$state, levels = c("TssActive", "TssBivalent", "Transcript", "Enhancer", "Repressive", "Heterochromatin", "Quiescent"))
df_states_per_cluster_dnd41$state <- factor(df_states_per_cluster_dnd41$state, levels = c("TssActive", "TssBivalent", "Transcript", "Enhancer", "Repressive", "Heterochromatin", "Quiescent"))

states_col <- brewer.pal(7, "Set2")
names(states_col) <- c("TssActive", "Transcript", "Enhancer", "Heterochromatin", "TssBivalent", "Repressive", "Quiescent")

pdf("Figure4e.pdf", width = 9)
ggplot(df_states_per_cluster_hsc, aes(x = cluster, y = Freq, fill = state)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = states_col) + ylab("Fraction") + xlab("") + theme_classic() + ggtitle("HSC")
ggplot(df_states_per_cluster_dnd41, aes(x = cluster, y = Freq, fill = state)) + geom_bar(position="fill", stat="identity") + scale_fill_manual(values = states_col) + ylab("Fraction") + xlab("") + theme_classic() + ggtitle("DND41")
dev.off()

## Chi square test of homogeneity
# Null hypothesis: distributions of states per cluster are the same for HSC and DND41
# Alternative Hypothesis: distributions of states per cluster are not the same

mat_chisq_test <- lapply(c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high", "all_covered"), function(x)
{
    data <- data.frame("HSC" = subset(df_states_per_cluster_hsc, cluster == x)$Freq, "DND41" = subset(df_states_per_cluster_dnd41, cluster == x)$Freq)
    keep_rows <- !apply(data, 1, function(y) all(y == 0))
    data <- data[keep_rows,]
    return(t(data))
})
names(mat_chisq_test) <- c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high", "all_covered")

lapply(mat_chisq_test, function(x) chisq.test(x))
unlist(lapply(mat_chisq_test, function(x) chisq.test(x)$p.value))

lapply(mat_chisq_test, function(x) cramer.test(x))
unlist(lapply(mat_chisq_test, function(x) cramer.test(x)$estimate))

################################################################################
## Figure 4f
################################################################################

median_per_cluster_pancancer_df <- read.table(file = "Figure4f_source_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
median_per_cluster_pancancer_df$Var2 <- factor(median_per_cluster_pancancer_df$Var2, levels = rev(c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high")))

colors_pancancer_type <- c("grey35", "black", "firebrick", "deepskyblue1", "royalblue2", "#984EA3", brewer.pal(6, "Blues")[2:6], "grey65", brewer.pal(8, "Set3"))
names(colors_pancancer_type) <- c("precursor B cells", "precursor T cells", "T-ALL", "DUX4/ERG", "hypo-diploid", "Ph-like ALL", "B-ALL Blueprint", "AML", "TPLL", "CLL", "MCL", "control_solid", "BRCA", "UCEC", "LUSC", "LUAD", "COAD", "STAD", "READ", "BLCA")

median_per_cluster_pancancer_df$Var1 <- factor(median_per_cluster_pancancer_df$Var1, levels = names(colors_pancancer_type))

pdf("Figure4f.pdf", width = 9, height = 7)
ggplot(median_per_cluster_pancancer_df, aes(x = Var1, y = Var2, color = Var1, size = value)) + geom_point() + ylab("Median CGI methylation") + xlab("") + theme_classic() + scale_color_manual(values = colors_pancancer_type) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=50, hjust=1)) + theme(legend.position="bottom") + lims(size = c(0, 1))
dev.off()

################################################################################
## Extended Data Figure 4b
################################################################################

cgi_annotation_complete <- read.table("Figure4c_Extended_Data_Figure4b_source_data_annotation_cgi.tsv", header = TRUE, row.names = 1, sep = "\t")
cgi_annotation_complete$cluster <- factor(cgi_annotation_complete$cluster, levels = c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high"))

pdf("Extended_Data_Figure4b.pdf", height = 5)
ggplot(data = cgi_annotation_complete, aes(x = num_cpg, color = cluster)) + geom_line(stat="density", size = 1) + theme_classic() + xlab("Number of CpGs in CGI") + ylab("Density") + scale_color_manual(values = rev(brewer.pal(6, "RdYlBu"))) + coord_cartesian(xlim = c(0, 500))

ggplot(data = cgi_annotation_complete, aes(x = length, color = cluster)) + geom_line(stat="density", size = 1) + theme_classic() + xlab("Length CGI") + ylab("Density") + scale_color_manual(values = rev(brewer.pal(6, "RdYlBu"))) + coord_cartesian(xlim = c(0, 4000))

ggplot(data = cgi_annotation_complete, aes(x = contentGC, color = cluster)) + geom_line(stat="density", size = 1) + theme_classic() + xlab("GC content CGI") + ylab("Density") + scale_color_manual(values = rev(brewer.pal(6, "RdYlBu")))
dev.off()

################################################################################
## Extended Data Figure 4c
################################################################################

cgi_gene_tpm_df_promoter <- read.table(file = "Extended_Data_Figure4c_source_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cgi_gene_tpm_df_promoter$cluster <- factor(cgi_gene_tpm_df_promoter$cluster, levels = c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high"))

pdf("Extended_Data_Figure4c.pdf")
ggplot(cgi_gene_tpm_df_promoter, aes(x=cluster, y=log2_tpm)) + geom_boxplot(outlier.shape = NA) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=50, hjust=1)) + xlab("") + ylab("Log2 (TPM + 1)") + coord_cartesian(ylim = c(0, 10)) + geom_hline(yintercept = 1, lty = 2, color = "grey")
dev.off()

################################################################################
## Extended Data Figure 4d
################################################################################

cgi_gene_tpm_df_active_promoter <- read.table(file = "Extended_Data_Figure4d_source_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

active_genes_promoter_cgi_per_cluster <- lapply(sort(unique(as.character(cgi_gene_tpm_df_active_promoter$cluster))), function(x) sort(as.character(subset(cgi_gene_tpm_df_active_promoter, cluster == x)$gene_name)))
names(active_genes_promoter_cgi_per_cluster) <- sort(unique(as.character(cgi_gene_tpm_df_active_promoter$cluster)))

ora <- lapply(active_genes_promoter_cgi_per_cluster, function(x) WebGestaltR(interestGene = x, enrichMethod = "ORA", enrichDatabase = "geneontology_Biological_Process", organism = "hsapiens", referenceSet = "genome", isOutput = FALSE, interestGeneType="genesymbol", referenceGeneType="genesymbol",  maxNum = 500, sigMethod = "top", topThr = 20))

dotplot_enrichment <- function(data, top, selection = NULL, analysis = "ora")
{
    if (length(selection) > 0)
    {
        subset_data <- subset(data, geneSet %in% selection)
    }
    else
    {
        subset_data <- data
    }

    subset_data <- subset_data[order(subset_data$FDR),]
    subset_data <- head(subset_data, top)

    subset_data <- subset_data[order(subset_data$enrichmentRatio, decreasing = TRUE),]
    subset_data$description <- factor(subset_data$description, levels = rev(subset_data$description))
    subset_data$gene_ratio <- subset_data$overlap / subset_data$size

    p <- ggplot(subset_data, aes(x = enrichmentRatio, y = description, size = gene_ratio, color = FDR)) + geom_point() + scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.05)) + theme_light() + ylab("") + xlab("Enrichment ratio") + guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) + lims(size = c(0, 1))
    return(p)
}

dotplots_ora <- lapply(1:length(ora), function(x) dotplot_enrichment(ora[[x]], top = 20))
names(dotplots_ora) <- names(ora)

pdf("Extended_Data_Figure4d.pdf", width = 18)
p1 <- dotplots_ora[["low"]] + ggtitle("Cluster low")
p2 <- dotplots_ora[["cluster1"]] + ggtitle("Cluster 1")
grid.arrange(p1, p2, ncol = 2)
dev.off()

################################################################################
## Extended Data Figure 4e
################################################################################

avg_cgi_pancancer_df <- read.table(file = "Extended_Data_Figure4e_source_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

colors_pancancer_type <- c("grey35", "grey35", "firebrick", "deepskyblue1", "royalblue2", "#984EA3", brewer.pal(6, "Blues")[2:6], "grey35", brewer.pal(8, "Set3"))
names(colors_pancancer_type) <- c("precursor B cells", "precursor T cells", "T-ALL", "DUX4/ERG", "hypo-diploid", "Ph-like ALL", "B-ALL Blueprint", "AML", "TPLL", "CLL", "MCL", "control_solid", "BRCA", "UCEC", "LUSC", "LUAD", "COAD", "STAD", "READ", "BLCA")

avg_cgi_pancancer_df$cluster <- factor(avg_cgi_pancancer_df$cluster, levels = c("low", "cluster1", "cluster2", "cluster3", "cluster4", "high"))
avg_cgi_pancancer_df$variable <- factor(avg_cgi_pancancer_df$variable, levels = names(colors_pancancer_type))

pdf("Extended_Data_Figure4e.pdf", height = 10)
ggplot(avg_cgi_pancancer_df, aes(x=variable, y=value, fill = variable)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1)) + theme_classic() + scale_fill_manual(values = colors_pancancer_type) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=50, hjust=1)) + xlab("") + ylab("Mean CGI methylation") + facet_grid(cluster~.)
dev.off()
