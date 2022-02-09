################################################################################
## This script contains the code to produce Figure 7 and Extended Data Figure 7
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
library(ggridges)
library(ggrepel)
theme_set(theme_ridges())

################################################################################
## Figure 7a
################################################################################

most_variable_cpgs_with_ccl <- data.frame(fread("Figure7a_source_data.tsv", header = TRUE, stringsAsFactors = FALSE), check.names = FALSE)
metadata_ccl <- read.table("Figure7a_source_data_annotation.tsv", row.names = 1, header = TRUE, check.names = FALSE, sep = "\t")

most_variable_cpgs_with_ccl <- most_variable_cpgs_with_ccl[,-c(1:3)]
dist_samples_with_ccl <- dist(t(most_variable_cpgs_with_ccl))

annotation_colors_with_ccl <- list(Type = c("black", "gray50", brewer.pal(3, "Reds"), "deepskyblue1", "royalblue2", "#984EA3", "turquoise3", "darkorange"), Sex = c("gold", "royalblue4"), Age_Group = c(brewer.pal(n = 3, name = "YlGn")))
names(annotation_colors_with_ccl[[1]]) <- c("control T cells", "control B cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "DUX4/ERG", "hypo-diploid", "Ph-like ALL", "B-ALL cell line", "T-ALL cell line")
names(annotation_colors_with_ccl[[2]]) <- c("female", "male")
names(annotation_colors_with_ccl[[3]]) <- c("Pediatric", "AYA", "Adult")

breaksList <- seq(0, 800, by = 1)
pdf("Figure7a.pdf", width = 14, height = 12)
pheatmap::pheatmap(dist_samples_with_ccl, annotation_row = metadata_ccl[colnames(most_variable_cpgs_with_ccl),c("Type", "Sex", "Age_Group")], annotation_colors = annotation_colors_with_ccl, breaks = breaksList, color = colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(length(breaksList)), labels_row = metadata_ccl[colnames(most_variable_cpgs_with_ccl),"label"])
dev.off()

################################################################################
## Figure 7b
################################################################################

avg_cgi_avg_subtypes_cell_lines_tall_variable <- read.table("Figure7b_source_data.tsv", row.names = 1, header = TRUE, check.names = FALSE, sep = "\t")

avg_cgi_avg_subtypes_cell_lines_tall_variable_df <- melt(avg_cgi_avg_subtypes_cell_lines_tall_variable)
avg_cgi_avg_subtypes_cell_lines_tall_variable_df$variable <- factor(avg_cgi_avg_subtypes_cell_lines_tall_variable_df$variable, levels = c("precursor T cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "MOLT-16", "Jurkat", "PEER", "DND41", "PER-117", "RMPI-8402", "LOUCY", "TALL-1", "ALL-SIL"))

pdf("Figure7b.pdf", height = 5, width = 12)
par(mar=c(9, 4, 4, 2))
vioplot(value ~ variable, data = avg_cgi_avg_subtypes_cell_lines_tall_variable_df, ylim = c(0,1), las = 2, xlab = "", ylab = "Mean methylation variable CGIs", col = c("gray50", brewer.pal(4, "Reds")[2:4], rep("darkorange", 9)))
dev.off()

################################################################################
## Figure 7c
################################################################################

counts <- read.table("Figure7c_Extended_Data_Figure7h_source_data_counts.tsv", header = TRUE, row.names = 1)
tpm <- read.table("Figure7c_Extended_Data_Figure7h_source_data_tpm.tsv", header = TRUE, row.names = 1)
annotation_genes <- read.table("Figure7c_Extended_Data_Figure7h_source_data_annotation.tsv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

counts_comparison <- counts[,c("Jurkat_R1", "Jurkat_R2", "Jurkat_R3", "DND41_R1", "DND41_R2", "DND41_R3")]
tpm_comparison <- tpm[,c("Jurkat_R1", "Jurkat_R2", "Jurkat_R3", "DND41_R1", "DND41_R2", "DND41_R3")]
active_genes_comparison <- names(which(apply(tpm_comparison, 1, mean, na.rm = TRUE) >= 0.5))
tpm_active_genes_comparison <- tpm_comparison[active_genes_comparison,]
col_data <- data.frame(row.names = colnames(tpm), sample = c(rep("Jurkat", 3), rep("DND41", 3), rep("Jurkat_TET2KO", 3)))

## Differential expression

dds <- DESeqDataSetFromMatrix(countData = counts_comparison, colData = col_data[1:6,,drop=FALSE], design = ~ sample)
keep <- rowSums(counts(dds)) >= 10
table(keep)
dds <- dds[keep,]
dds <- DESeq(dds)

results <- results(dds, contrast = c("sample", "Jurkat", "DND41"))
results_significant <- subset(results, (padj < 0.05) & (abs(log2FoldChange) > 1))
results_significant_filtered <- results_significant[intersect(rownames(results_significant), rownames(tpm_active_genes_comparison)),]
results_filtered <- results[intersect(rownames(results), rownames(tpm_active_genes_comparison)),]

results_significant_filtered_up <- results_significant_filtered[results_significant_filtered$log2FoldChange > 1,]
results_significant_filtered_down <- results_significant_filtered[results_significant_filtered$log2FoldChange < -1,]

## Plot epigenetic regulator expression

epi_regulators <-c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3", "MYC", "IDH1", "IDH2", "WT1", "EED", "EZH2", "SUZ12", "RING1", "RNF2", "KDM2B", "BAP1", "SUV39H1", "SUV39H2", "EHMT2", "ARID1A", "ARID1B", "ARID2", "PBRM1", "SMARCA4", "SMARCB1", "QSER1", "HELLS")
epi_regulators_id <- sapply(epi_regulators, function(x) rownames(subset(annotation_genes, gene_name == x)))

intersect(epi_regulators_id, rownames(results_significant_filtered_up)) ## TET2
intersect(epi_regulators_id, rownames(results_significant_filtered_down)) ## DNMT1, DNMT3B, WT1

col_fun_rna_log2 <- colorRamp2(c(1,5,9), brewer.pal(9, "RdYlBu")[c(9,5,1)])

pdf("Figure7c.pdf", width = 4)
Heatmap(as.matrix(log2(tpm_active_genes_comparison + 1)[epi_regulators_id,]), show_row_names = TRUE, row_labels = epi_regulators, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_rna_log2)
dev.off()

################################################################################
## Figure 7d
################################################################################

methyl_data_ccl <- data.frame(fread("Figure7d_source_data.tsv", stringsAsFactors = FALSE, header = TRUE), check.names = FALSE)

pdf("Figure7d_densityplot.pdf", width = 15, height = 5)
layout(matrix(1:3, ncol = 3, byrow = TRUE))
smoothScatter(methyl_data_ccl[,"Jurkat"], methyl_data_ccl[,"Jurkat_TET2_KO"], xlab = "Jurkat_MPIMG", ylab = "Jurkat_TET2_KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
smoothScatter(methyl_data_ccl[,"PEER"], methyl_data_ccl[,"Jurkat_TET2_KO"], xlab = "PEER", ylab = "Jurkat_TET2_KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
smoothScatter(methyl_data_ccl[,"DND41"], methyl_data_ccl[,"Jurkat_TET2_KO"], xlab = "DND41", ylab = "Jurkat_TET2_KO", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
dev.off()

diff_per_cpg <- data.frame(Jurkat_TET2KO_vs_Jurkat = methyl_data_ccl[,"Jurkat_TET2_KO"] - methyl_data_ccl[,"Jurkat"],
                           Jurkat_TET2KO_vs_PEER = methyl_data_ccl[,"Jurkat_TET2_KO"] - methyl_data_ccl[,"PEER"],
                           Jurkat_TET2KO_vs_DND41 = methyl_data_ccl[,"Jurkat_TET2_KO"] - methyl_data_ccl[,"DND41"])

count_per_cpg_class <- matrix(NA, nrow = 3, ncol = 3)
colnames(count_per_cpg_class) <- c("Jurkat_TET2KO_vs_Jurkat","Jurkat_TET2KO_vs_PEER", "Jurkat_TET2KO_vs_DND41")
rownames(count_per_cpg_class) <- c("hyper", "stable", "hypo")
count_per_cpg_class["hyper",] <- apply(diff_per_cpg, 2, function(x) length(which(x > 0.1)))
count_per_cpg_class["hypo",] <- apply(diff_per_cpg, 2, function(x) length(which(x < -0.1)))
count_per_cpg_class["stable",] <- apply(diff_per_cpg, 2, function(x) length(which(x >= -0.1 & x <= 0.1)))
count_per_cpg_class <- data.frame(count_per_cpg_class)
count_per_cpg_class$feature <- rownames(count_per_cpg_class)
count_per_cpg_class_df <- melt(count_per_cpg_class)
count_per_cpg_class_df$feature <- factor(count_per_cpg_class_df$feature, levels = c("hypo", "stable", "hyper"))

pdf("Figure7d_barplot.pdf", width = 5)
ggplot(data = count_per_cpg_class_df, aes(x = variable, fill = feature, y = value)) + geom_bar(position = "fill", stat = "identity") + theme_classic() + xlab("") + ylab("Fraction CpGs") + theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=12)) + theme(axis.text.x=element_text(hjust = 1, angle = 45, size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=12)) + scale_fill_manual(values = c("darkblue", "grey", "darkred"))
dev.off()

################################################################################
## Figure 7e
################################################################################

diff_avg_sw_ccl_pmd_hmd <- data.frame(fread("Figure7e_source_data.tsv", header = TRUE), stringsAsFactors = FALSE, row.names = 1)
diff_avg_sw_ccl_pmd_hmd_df <- melt(diff_avg_sw_ccl_pmd_hmd)
diff_avg_sw_ccl_pmd_hmd_df$variable <- factor(diff_avg_sw_ccl_pmd_hmd_df$variable, levels = c(rev(c("Jurkat_TET2KO_vs_Jurkat", "Jurkat_vs_PEER", "Jurkat_TET2KO_vs_PEER", "Jurkat_vs_DND41", "Jurkat_TET2KO_vs_DND41"))))

diff_avg_sw_ccl_pmd_hmd_df$feature <- factor(diff_avg_sw_ccl_pmd_hmd_df$feature, levels = c("PMD", "HMD"))

pdf("Figure7e.pdf", width = 8, height = 6)
ggplot(diff_avg_sw_ccl_pmd_hmd_df, aes(x=variable, y=value, fill = feature)) + theme_classic() + geom_boxplot(outlier.shape = NA, alpha = 0.5) + theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=14)) + ylab("Difference sliding window") + xlab("") + coord_flip(ylim=c(-1, 0.5)) + geom_hline(yintercept = 0.1, lty = 2) + geom_hline(yintercept = -0.1, lty = 2) + geom_hline(yintercept = 0, lty = 2) + scale_fill_manual(values = c("darkblue", "forestgreen"))
dev.off()

################################################################################
## Figure 7f
################################################################################

avg_features_cell_line_ko <- data.frame(fread(file = "Figure7f_source_data.tsv", header = TRUE), stringsAsFactors = FALSE, row.names = 1)

avg_features_cell_line_ko_df <- melt(avg_features_cell_line_ko)
avg_features_cell_line_ko_df$variable <- factor(avg_features_cell_line_ko_df$variable, levels = c("DND41", "PEER", "Jurkat_TET2_KO", "Jurkat"))

pdf("Figure7f.pdf", height = 7, width = 7)
ggplot(avg_features_cell_line_ko_df, aes(y=variable, x=value, fill = feature)) + geom_density_ridges(aes(fill = feature), rel_min_height=.01, alpha = 0.5, panel_scaling = FALSE, scale = 1) + coord_cartesian(xlim = c(0, 1)) + theme_classic() + theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) + ylab("") + xlab("Mean methylation") + scale_fill_manual(values = c("forestgreen", "darkblue", "darkorange"))
dev.off()

################################################################################
## Extended Data Figure 7a
################################################################################

avg_cgi_common_tall_with_ccl_variable <- data.frame(fread("Extended_Data_Figure7a_source_data.tsv", header = TRUE), stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

samplesControlTALL <- c("SJNORM016314_G1", "SJNORM016314_G2", "SJNORM016314_G3", "SJNORM016314_G4", "SJNORM016314_G5", "SJNORM016315_G1", "SJNORM016315_G2", "SJNORM016315_G3", "SJNORM016315_G4", "SJNORM016315_G5")
samplesTALL <- c("SJALL030551_D1", "SJALL030933_D1", "SJALL048427_D3", "SJALL048442_D2", "SJALL064568_D1", "SJALL064571_D1", "SJALL064572_D1", "SJALL064573_D1", "SJTALL004_D", "SJTALL010_D", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012402_D1", "SJTALL022646_D1", "SJTALL187_D", "SJALL030493_D1", "SJALL030684_D1", "SJALL048428_D4", "SJALL048440_D2", "SJALL048449_D2", "SJALL064567_D1", "SJALL064569_D1", "SJTALL012401_D1", "SJTALL012403_D1", "SJTALL015715_D1", "SJTALL021687_D1", "SJTALL021688_D1", "SJTALL021912_D1", "SJTALL022648_D1", "SJTALL022649_D1", "SJTALL170_D", "SJTALL182_D", "SJTALL186_D", "SJALL030916_D1", "SJALL059529_D2", "SJALL064564_D1", "SJALL064565_D1", "SJALL064566_D1", "SJTALL009_D", "SJTALL021685_D1", "SJTALL021694_D1", "SJTALL021695_D1", "SJTALL022645_D1", "SJTALL022647_D1", "SJTALL169_D", "SJTALL174_D")

dataForPcaCCL <- t(avg_cgi_common_tall_with_ccl_variable[,c(samplesControlTALL, samplesTALL)])
pcaCCL <- prcomp(dataForPcaCCL)

dfPCACCL <- data.frame(pcaCCL$x, type = factor(ifelse(rownames(pcaCCL$x) %in% samplesControlTALL, "control T cells", "T-ALL"), levels = c("control T cells", "T-ALL")), median_variable_cgi = apply(avg_cgi_common_tall_with_ccl_variable[,c(samplesControlTALL, samplesTALL)], 2, median), fraction_variable_cgi = colSums(ifelse(avg_cgi_common_tall_with_ccl_variable[,c(samplesControlTALL, samplesTALL)] > 0.2, 1, 0)) / nrow(avg_cgi_common_tall_with_ccl_variable))

percentageCCL <- summary(pcaCCL)$importance["Proportion of Variance",]
percentageCCL <- paste(colnames(dfPCACCL), " (", as.character(round(percentageCCL * 100, 2)), "%", ")", sep="")

## Project cell line samples into PCA
prediction_tall_cell_lines <- predict(pcaCCL, t(avg_cgi_common_tall_with_ccl_variable[,!colnames(avg_cgi_common_tall_with_ccl_variable) %in% c(samplesControlTALL, samplesTALL)]))

dfPCACCL_2 <- data.frame(prediction_tall_cell_lines, type = rep("T-ALL cell line", nrow(prediction_tall_cell_lines)), median_variable_cgi = apply(avg_cgi_common_tall_with_ccl_variable[,!colnames(avg_cgi_common_tall_with_ccl_variable) %in% c(samplesControlTALL, samplesTALL)], 2, median), fraction_variable_cgi = colSums(ifelse(avg_cgi_common_tall_with_ccl_variable[,!colnames(avg_cgi_common_tall_with_ccl_variable) %in% c(samplesControlTALL, samplesTALL)] > 0.2, 1, 0)) / nrow(avg_cgi_common_tall_with_ccl_variable))

dfPCACCL$sample <- NA
dfPCACCL_2$sample <- rownames(dfPCACCL_2)
dfPCACCL <- rbind(dfPCACCL, dfPCACCL_2)

pdf("Extended_Data_Figure7a_mean.pdf", width = 6.5, height = 5)
ggplot(dfPCACCL,aes(x=PC1,y=PC2, col = median_variable_cgi, shape = type, label = sample))+
       geom_point(size=2) +
       xlab(percentageCCL[1]) + ylab(percentageCCL[2]) +
       scale_color_gradientn(colours = c("#0072B2", "seagreen3", "#F0E442"), limits=c(0,1)) +
       theme_classic() + geom_text_repel(vjust="inward",hjust="inward", size = 3)
dev.off()

## Repeat with binary

dataForPcaCCLBinary <- apply(dataForPcaCCL, 2, function(x) ifelse(x > 0.2, 1, 0))
pcaCCLBinary <- prcomp(dataForPcaCCLBinary)

dfPCACCLBinary <- data.frame(pcaCCLBinary$x, type = factor(ifelse(rownames(pcaCCL$x) %in% samplesControlTALL, "control T cells", "T-ALL"), levels = c("control T cells", "T-ALL")), median_variable_cgi = apply(avg_cgi_common_tall_with_ccl_variable[,c(samplesControlTALL, samplesTALL)], 2, median), fraction_variable_cgi = colSums(ifelse(avg_cgi_common_tall_with_ccl_variable[,c(samplesControlTALL, samplesTALL)] > 0.2, 1, 0)) / nrow(avg_cgi_common_tall_with_ccl_variable))

percentageCCLBinary <- summary(pcaCCLBinary)$importance["Proportion of Variance",]
percentageCCLBinary <- paste(colnames(dfPCACCLBinary), " (", as.character(round(percentageCCLBinary * 100, 2)), "%", ")", sep="")

## Project cell line samples into PCA
prediction_tall_cell_lines_binary <- predict(pcaCCLBinary, apply(t(avg_cgi_common_tall_with_ccl_variable[,!colnames(avg_cgi_common_tall_with_ccl_variable) %in% c(samplesControlTALL, samplesTALL)]), 2, function(x) ifelse(x > 0.2, 1, 0)))

dfPCACCLBinary_2 <- data.frame(prediction_tall_cell_lines_binary, type = rep("T-ALL cell line", nrow(prediction_tall_cell_lines)), median_variable_cgi = apply(avg_cgi_common_tall_with_ccl_variable[,!colnames(avg_cgi_common_tall_with_ccl_variable) %in% c(samplesControlTALL, samplesTALL)], 2, median), fraction_variable_cgi = colSums(ifelse(avg_cgi_common_tall_with_ccl_variable[,!colnames(avg_cgi_common_tall_with_ccl_variable) %in% c(samplesControlTALL, samplesTALL)] > 0.2, 1, 0)) / nrow(avg_cgi_common_tall_with_ccl_variable))

dfPCACCLBinary$sample <- NA
dfPCACCLBinary_2$sample <- rownames(dfPCACCLBinary_2)
dfPCACCLBinary <- rbind(dfPCACCLBinary, dfPCACCLBinary_2)

pdf("Extended_Data_Figure7a_binary.pdf", width = 6.5, height = 5)
ggplot(dfPCACCLBinary,aes(x=PC1,y=PC2, col = fraction_variable_cgi, shape = type, label = sample))+
      geom_point(size=2) +
      xlab(percentageCCLBinary[1]) + ylab(percentageCCLBinary[2]) +
      scale_color_gradientn(colours = c("yellow1", "#FF8C00", "red3"), limits=c(0,1)) +
      theme_classic() + geom_text_repel(vjust="inward",hjust="inward", size = 3)
dev.off()

################################################################################
## Extended Data Figure 7b
################################################################################

avg_tile_avg_subtypes_cell_lines <- data.frame(fread("Extended_Data_Figure7b_source_data.tsv", header = TRUE), stringsAsFactors = FALSE, row.names = 1, check.names = FALSE)

avg_tile_avg_subtypes_cell_lines_df <- melt(avg_tile_avg_subtypes_cell_lines)
avg_tile_avg_subtypes_cell_lines_df$variable <- factor(avg_tile_avg_subtypes_cell_lines_df$variable, levels = c("precursor T cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "MOLT-16", "Jurkat", "PEER", "DND41", "PER-117", "RMPI-8402", "LOUCY", "TALL-1", "ALL-SIL", "precursor B cells", "DUX4/ERG", "NALM6", "hypo-diploid", "MHH-CALL-2", "NALM16", "Ph-like ALL", "MHH-CALL-4", "MUTZ5"))

pdf("Extended_Data_Figure7b.pdf", height = 5, width = 16)
par(mar=c(9, 4, 4, 2))
vioplot(value ~ variable, data = subset(avg_tile_avg_subtypes_cell_lines_df, feature == "HMD"), ylim = c(0,1), las = 2, xlab = "", ylab = "Mean methylation 1kb tiles", col = c("gray50", brewer.pal(4, "Reds")[-1], rep("darkorange", 9), "gray50", "deepskyblue1", "turquoise3", "royalblue2", "turquoise3", "turquoise3", "#984EA3", "turquoise3", "turquoise3"), side = "left", plotCentre = "line")
vioplot(value ~ variable, data = subset(avg_tile_avg_subtypes_cell_lines_df, feature == "PMD"), ylim = c(0,1), las = 2, xlab = "", ylab = "Mean methylation 1kb tiles", col = c("gray50", brewer.pal(4, "Reds")[-1], rep("darkorange", 9), "gray50", "deepskyblue1", "turquoise3", "royalblue2", "turquoise3", "turquoise3", "#984EA3", "turquoise3", "turquoise3"), side = "right", plotCentre = "line", add = TRUE)
dev.off()

################################################################################
## Extended Data Figure 7c
################################################################################

avg_cgi_avg_subtypes_cell_lines_ball_variable <- read.table("Extended_Data_Figure7c_source_data.tsv", row.names = 1, header = TRUE, check.names = FALSE, sep = "\t")

avg_cgi_avg_subtypes_cell_lines_ball_variable_df <- melt(avg_cgi_avg_subtypes_cell_lines_ball_variable)
avg_cgi_avg_subtypes_cell_lines_ball_variable_df$variable <- factor(avg_cgi_avg_subtypes_cell_lines_ball_variable_df$variable, levels = c("precursor B cells", "DUX4/ERG", "NALM6", "hypo-diploid", "MHH-CALL-2", "NALM16", "Ph-like ALL", "MHH-CALL-4", "MUTZ5"))

pdf("Extended_Data_Figure7c.pdf", height = 5, width = 10)
par(mar=c(9, 4, 4, 2))
vioplot(value ~ variable, data = avg_cgi_avg_subtypes_cell_lines_ball_variable_df, ylim = c(0,1), las = 2, xlab = "", ylab = "Mean methylation variable CGIs", col = c("gray50", "deepskyblue1", "turquoise3", "royalblue2", "turquoise3", "turquoise3", "#984EA3", "turquoise3", "turquoise3"))
dev.off()

################################################################################
## Extended Data Figure 7d
################################################################################

avg_promoter_5x_epi_reg_cell_lines <- read.table("Extended_Data_Figure7d_source_data.tsv", row.names = 1, header = TRUE, check.names = FALSE, sep = "\t")

col_fun_dname <- colorRamp2(c(0.1, 0.5, 0.9), brewer.pal(9, "RdBu")[c(9,5,1)])

pdf("Extended_Data_Figure7d.pdf", width = 6, height = 8)
Heatmap(as.matrix(avg_promoter_5x_epi_reg_cell_lines), show_row_names = TRUE, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_dname)
dev.off()

################################################################################
## Extended Data Figure 7h
################################################################################

counts <- read.table("Figure7c_Extended_Data_Figure7h_source_data_counts.tsv", header = TRUE, row.names = 1)
tpm <- read.table("Figure7c_Extended_Data_Figure7h_source_data_tpm.tsv", header = TRUE, row.names = 1)
annotation_genes <- read.table("Figure7c_Extended_Data_Figure7h_source_data_annotation.tsv", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

counts_ko <- counts[,c("Jurkat_R1", "Jurkat_R2", "Jurkat_R3", "Jurkat_TET2KO_R1", "Jurkat_TET2KO_R2", "Jurkat_TET2KO_R3")]
tpm_ko <- tpm[,c("Jurkat_R1", "Jurkat_R2", "Jurkat_R3", "Jurkat_TET2KO_R1", "Jurkat_TET2KO_R2", "Jurkat_TET2KO_R3")]
active_genes_ko <- names(which(apply(tpm_ko, 1, mean, na.rm = TRUE) >= 0.5))
tpm_active_genes_ko <- tpm_ko[active_genes_ko,]

col_data <- data.frame(row.names = colnames(tpm), sample = c(rep("Jurkat", 3), rep("DND41", 3), rep("Jurkat_TET2KO", 3)))

## Differential expression

dds_ko <- DESeqDataSetFromMatrix(countData = counts_ko, colData = col_data[c(1:3,7:9),,drop=FALSE], design = ~ sample)
keep <- rowSums(counts(dds_ko)) >= 10
table(keep)
dds_ko <- dds_ko[keep,]
dds_ko <- DESeq(dds_ko)

results_ko <- results(dds_ko, contrast = c("sample", "Jurkat_TET2KO", "Jurkat"))
results_ko_significant <- subset(results_ko, (padj < 0.05) & (abs(log2FoldChange) > 1))
results_ko_significant_filtered <- results_ko_significant[intersect(rownames(results_ko_significant), rownames(tpm_active_genes_ko)),]
results_ko_filtered <- results[intersect(rownames(results), rownames(tpm_active_genes_ko)),]

results_ko_significant_filtered_up <- results_ko_significant_filtered[results_ko_significant_filtered$log2FoldChange > 1,]
results_ko_significant_filtered_down <- results_ko_significant_filtered[results_ko_significant_filtered$log2FoldChange < -1,]

## Plot epigenetic regulator expression

epi_regulators <-c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3", "MYC", "IDH1", "IDH2", "WT1", "EED", "EZH2", "SUZ12", "RING1", "RNF2", "KDM2B", "BAP1", "SUV39H1", "SUV39H2", "EHMT2", "ARID1A", "ARID1B", "ARID2", "PBRM1", "SMARCA4", "SMARCB1", "QSER1", "HELLS")
epi_regulators_id <- sapply(epi_regulators, function(x) rownames(subset(annotation_genes, gene_name == x)))

intersect(epi_regulators_id, rownames(results_ko_significant_filtered_up)) ## DNMT3B
intersect(epi_regulators_id, rownames(results_ko_significant_filtered_down)) ## TET2

col_fun_rna_log2 <- colorRamp2(c(1,5,9), brewer.pal(9, "RdYlBu")[c(9,5,1)])

pdf("Extended_Data_Figure7h.pdf", width = 4)
Heatmap(as.matrix(log2(tpm_active_genes_ko + 1)[epi_regulators_id,]), show_row_names = TRUE, row_labels = epi_regulators, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_rna_log2)
dev.off()
