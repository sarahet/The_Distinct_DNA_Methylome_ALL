################################################################################
## This script contains the code to produce Figure 6 and Extended Data Figure 6
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
library(WebGestaltR)

################################################################################
## Figure 6a + Extended Data Figure 6c
################################################################################

tpm_active_genes_for_correlation_test <- read.table(file = "Figure6a_source_data.tsv", sep = "\t", row.names = 1, header = TRUE)
cor_cgi_dname_genes_with_significance <- read.table(file = "Figure6a_Extended_Data_Figure_6a_source_data_cor_test_cgi.tsv", sep = "\t", row.names = 1, header = TRUE)
cor_global_dname_genes_with_significance <- read.table(file = "Figure6a_Extended_Data_Figure_6a_source_data_cor_test_global.tsv", sep = "\t", row.names = 1, header = TRUE)
df_avg_methylation_per_sample_cor <- read.table(file = "Figure6a_source_data_annotation.tsv", sep = "\t", row.names = 1, header = TRUE)

## Epigenetic regulators
epi_reg <- c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3", "MYC", "IDH1", "IDH2", "WT1", "EED", "EZH2", "SUZ12", "RING1", "RNF2", "KDM2B", "BAP1", "SUV39H1", "SUV39H2", "EHMT2", "ARID1A", "ARID1B", "ARID2", "PBRM1", "SMARCA4", "SMARCB1", "QSER1", "HELLS")

annotation_colors <- c("black", "grey35", brewer.pal(4, "Reds")[2:4], "deepskyblue1", "royalblue2", "#984EA3")
names(annotation_colors) <- c("control T cells", "control B cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

pdf("Figure6a_Extended_Data_Figure6c_cgi.pdf")
for (i in sapply(epi_reg, function(x) rownames(subset(cor_cgi_dname_genes_with_significance, gene_name == x))))
{
    plot(log2(do.call(c, tpm_active_genes_for_correlation_test[i,]) + 1), df_avg_methylation_per_sample_cor[,"avg_cgi"], pch = 19, xlab = "Log2 (TPM + 1)", ylab = "Mean CGI DNAme", col = annotation_colors[as.character(df_avg_methylation_per_sample_cor[,"Type"])], main = paste(cor_cgi_dname_genes_with_significance[i,"gene_name"], "cor =", round(cor_cgi_dname_genes_with_significance[i,"rho"], 2), "padj =", round(cor_cgi_dname_genes_with_significance[i,"padj"], 6)))
}
dev.off()

pdf("Figure6a_Extended_Data_Figure6c_global.pdf")
for (i in sapply(epi_reg, function(x) rownames(subset(cor_global_dname_genes_with_significance, gene_name == x))))
{
    plot(log2(do.call(c, tpm_active_genes_for_correlation_test[i,]) + 1), df_avg_methylation_per_sample_cor[,"avg_global_no_cgi"], pch = 19, xlab = "Log2 (TPM + 1)", ylab = "Mean global no CGI DNAme", col = annotation_colors[as.character(df_avg_methylation_per_sample_cor[,"Type"])], main = paste(cor_global_dname_genes_with_significance[i,"gene_name"], "cor =", round(cor_global_dname_genes_with_significance[i,"rho"], 2), "padj =", round(cor_global_dname_genes_with_significance[i,"padj"], 6)))
}
dev.off()

################################################################################
## Figure 6b
################################################################################

avg_promoter_epi_reg <- read.table("Figure6b_6c_Extended_Data_Figure_6e_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata_rnaseq_with_avg_cgi_global <- read.table("Figure6b_Extended_Data_Figure6b_6d_6e_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

annotation_colors <- list(Type = c("black", "grey35", brewer.pal(4, "Reds")[2:4], "deepskyblue1", "royalblue2", "#984EA3"))
names(annotation_colors[[1]]) <- c("control T cells", "control B cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

col_fun_cgi <- colorRamp2(c(0.1, 0.5), c("#F7FCFD", "#810F7C"))
col_fun_global <- colorRamp2(c(0.8, 1), c("#F7FCFD", "#006D2C"))

annotation_colors$avg_cgi <- col_fun_cgi
annotation_colors$avg_global_no_cgi <- col_fun_global

col_fun_dname <- colorRamp2(c(0.1, 0.5, 0.9), brewer.pal(9, "RdBu")[c(9,5,1)])

order_samples_heatmap <- c("SJBALL218_D", "SJBALL230_D", "SJERG016_D", "SJERG003_D", "SJERG011_D", "SJERG020299_D1", "SJERG000020_D", "SJERG000009_D", "SJERG014_D", "SJERG006_D", "SJERG000030_D", "SJERG020306_D1", "SJBALL020852_D1", "SJBALL021530_D1", "SJBALL021644_D1", "SJBALL015_D", "SJBALL021080_D1", "SJBALL021109_D1", "SJBALL000010_D", "SJTALL004_D", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012402_D1", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL187_D", "SJTALL022646_D1", "SJTALL186_D", "SJTALL012403_D1", "SJTALL022648_D1", "SJTALL021687_D1", "SJTALL022649_D1", "SJTALL170_D", "SJTALL021912_D1", "SJTALL015715_D1", "SJTALL021688_D1", "SJTALL012401_D1", "SJTALL182_D", "SJTALL021685_D1", "SJTALL021695_D1", "SJTALL169_D", "SJTALL021694_D1", "SJTALL022647_D1", "SJTALL009_D", "SJTALL022645_D1", "SJTALL174_D")
order_control_samples_heatmap <- c("SJNORM016096_G1", "SJNORM016096_G2", "SJNORM016096_G3", "SJNORM016096_G4", "SJNORM016314_G4", "SJNORM016314_G5")

pdf("Figure6b.pdf", width = 12, height = 8)
h1 <- Heatmap(as.matrix(avg_promoter_epi_reg[,order_control_samples_heatmap]), show_row_names = TRUE, top_annotation = columnAnnotation(df = metadata_rnaseq_with_avg_cgi_global[order_control_samples_heatmap,c("Type", "avg_cgi", "avg_global_no_cgi")], col = annotation_colors),cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_dname)

h2 <- Heatmap(as.matrix(avg_promoter_epi_reg[,order_samples_heatmap]), show_row_names = TRUE, top_annotation = columnAnnotation(df = metadata_rnaseq_with_avg_cgi_global[order_samples_heatmap,c("Type", "avg_cgi", "avg_global_no_cgi")], col = annotation_colors), cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_dname)
draw(h1 + h2)
dev.off()

################################################################################
## Figure 6c
################################################################################

df_avg_methylation_per_sample <- read.table(file = "Figure6c_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")
avg_promoter_epi_reg <- read.table("Figure6b_6c_Extended_Data_Figure_6e_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata_rnaseq_with_avg_cgi_global <- read.table("Figure6b_Extended_Data_Figure6b_6d_6e_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

common_samples <- c("SJBALL000010_D", "SJBALL015_D", "SJBALL020852_D1", "SJBALL021080_D1", "SJBALL021109_D1", "SJBALL021530_D1", "SJBALL021644_D1", "SJBALL218_D", "SJBALL230_D", "SJERG000009_D", "SJERG000020_D", "SJERG000030_D", "SJERG003_D", "SJERG006_D", "SJERG011_D", "SJERG014_D", "SJERG016_D", "SJERG020299_D1", "SJERG020306_D1", "SJNORM016096_G1", "SJNORM016096_G2", "SJNORM016096_G3", "SJNORM016096_G4", "SJNORM016314_G4", "SJNORM016314_G5", "SJTALL004_D", "SJTALL009_D", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012401_D1", "SJTALL012402_D1", "SJTALL012403_D1", "SJTALL015715_D1", "SJTALL021685_D1", "SJTALL021687_D1", "SJTALL021688_D1", "SJTALL021694_D1", "SJTALL021695_D1", "SJTALL021912_D1", "SJTALL022645_D1", "SJTALL022646_D1", "SJTALL022647_D1", "SJTALL022648_D1", "SJTALL022649_D1", "SJTALL169_D", "SJTALL170_D", "SJTALL174_D", "SJTALL182_D", "SJTALL186_D", "SJTALL187_D")

wt1_mutations_samples <- c("SJTALL012397_D1", "SJTALL012399_D1", "SJTALL012401_D1", "SJTALL021694_D1", "SJTALL022649_D1", "SJTALL170_D")

annotation_colors <- list(Type = c("black", "grey35", brewer.pal(4, "Reds")[2:4], "deepskyblue1", "royalblue2", "#984EA3"))
names(annotation_colors[[1]]) <- c("control T cells", "control B cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

pdf("Figure6c.pdf")
plot(do.call(c, avg_promoter_epi_reg["TET2",common_samples]), df_avg_methylation_per_sample[common_samples,"avg_cgi"], pch = 19, ylab = "Mean CGI methylation", xlab = "Mean TET2 promoter methylation", col = annotation_colors$Type[as.character(df_avg_methylation_per_sample[common_samples,"Type"])], xlim = c(0,1))

plot(do.call(c, avg_promoter_epi_reg["TET2",common_samples]), df_avg_methylation_per_sample[common_samples,"avg_global_no_cgi"], pch = 19, ylab = "Mean global no CGI methylation", xlab = "Mean TET2 promoter methylation", col = annotation_colors$Type[as.character(df_avg_methylation_per_sample[common_samples,"Type"])], xlim = c(0,1))

plot(do.call(c, avg_promoter_epi_reg["WT1",common_samples]), df_avg_methylation_per_sample[common_samples,"avg_cgi"], pch = ifelse(common_samples %in% wt1_mutations_samples, 17, 19), ylab = "Mean CGI methylation", xlab = "Mean WT1 promoter methylation", col = annotation_colors$Type[as.character(df_avg_methylation_per_sample[common_samples,"Type"])], xlim = c(0,1))

plot(do.call(c, avg_promoter_epi_reg["WT1",common_samples]), df_avg_methylation_per_sample[common_samples,"avg_global_no_cgi"], pch = ifelse(common_samples %in% wt1_mutations_samples, 17, 19), ylab = "Mean global no CGI methylation", xlab = "Mean WT1 promoter methylation", col = annotation_colors$Type[as.character(df_avg_methylation_per_sample[common_samples,"Type"])], xlim = c(0,1))
dev.off()

################################################################################
## Extended Data Figure 6a
################################################################################

cor_cgi_dname_genes_with_significance <- read.table(file = "Figure6a_Extended_Data_Figure_6a_source_data_cor_test_cgi.tsv", sep = "\t", row.names = 1, header = TRUE)
cor_global_dname_genes_with_significance <- read.table(file = "Figure6a_Extended_Data_Figure_6a_source_data_cor_test_global.tsv", sep = "\t", row.names = 1, header = TRUE)

cor_genes_cgi <- as.character(subset(cor_cgi_dname_genes_with_significance, padj < 0.01)$gene_name)
cor_genes_global <- as.character(subset(cor_global_dname_genes_with_significance, padj < 0.01)$gene_name)

gene_sets <- list(cor_cgi = cor_genes_cgi, cor_global = cor_genes_global)

ora <- lapply(gene_sets, function(x) WebGestaltR(interestGene = x, enrichMethod = "ORA", enrichDatabase = "geneontology_Biological_Process", organism = "hsapiens", referenceSet = "genome", isOutput = FALSE, interestGeneType="genesymbol", referenceGeneType="genesymbol", minNum = 10, maxNum = 500, sigMethod = "top", topThr = 20))

ora[["cor_cgi"]]$sample <- "CGI"
ora[["cor_global"]]$sample <- "global"

candidates_cor_combined <- rbind(ora[["cor_cgi"]][,c("geneSet", "description", "FDR", "overlap", "size", "sample")], ora[["cor_global"]][,c("geneSet", "description", "FDR", "overlap", "size", "sample")])
candidates_cor_combined$sample <- factor(candidates_cor_combined$sample, levels = c("CGI", "global"))
candidates_cor_combined$description <- factor(candidates_cor_combined$description, sort(unique(candidates_cor_combined$description), decreasing = TRUE))
candidates_cor_combined$gene_ratio <- candidates_cor_combined$overlap / candidates_cor_combined$size

pdf("Extended_Data_Figure6a.pdf", width = 8)
ggplot(candidates_cor_combined, aes(x = sample, y = description, size = gene_ratio, color = FDR)) + geom_point() + scale_color_gradient(low = "red", high = "blue", limits = c(0, 0.05)) + theme_light() + ylab("") + xlab("") + guides(colour = guide_colourbar(order = 1), size = guide_legend(order = 2)) + lims(size = c(0, 1))
dev.off()

################################################################################
## Extended Data Figure 6b
################################################################################

tpm_epi_reg <- read.table("Extended_Data_Figure6b_6e_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata_rnaseq_with_avg_cgi_global <- read.table("Figure6b_Extended_Data_Figure6b_6d_6e_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

annotation_colors <- list(Type = c("black", "grey35", brewer.pal(4, "Reds")[2:4], "deepskyblue1", "royalblue2", "#984EA3"))
names(annotation_colors[[1]]) <- c("control T cells", "control B cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

col_fun_cgi <- colorRamp2(c(0.1, 0.5), c("#F7FCFD", "#810F7C"))
col_fun_global <- colorRamp2(c(0.8, 1), c("#F7FCFD", "#006D2C"))

annotation_colors$avg_cgi <- col_fun_cgi
annotation_colors$avg_global_no_cgi <- col_fun_global

order_samples_heatmap <- c("SJBALL218_D", "SJBALL230_D", "SJERG016_D", "SJERG003_D", "SJERG011_D", "SJERG020299_D1", "SJERG000020_D", "SJERG000009_D", "SJERG014_D", "SJERG006_D", "SJERG000030_D", "SJERG020306_D1", "SJBALL020852_D1", "SJBALL021530_D1", "SJBALL021644_D1", "SJBALL015_D", "SJBALL021080_D1", "SJBALL021109_D1", "SJBALL000010_D", "SJTALL004_D", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012402_D1", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL187_D", "SJTALL022646_D1", "SJTALL186_D", "SJTALL012403_D1", "SJTALL022648_D1", "SJTALL021687_D1", "SJTALL022649_D1", "SJTALL170_D", "SJTALL021912_D1", "SJTALL015715_D1", "SJTALL021688_D1", "SJTALL012401_D1", "SJTALL182_D", "SJTALL021685_D1", "SJTALL021695_D1", "SJTALL169_D", "SJTALL021694_D1", "SJTALL022647_D1", "SJTALL009_D", "SJTALL022645_D1", "SJTALL174_D")
order_control_samples_heatmap <- c("SJNORM016096_G1", "SJNORM016096_G2", "SJNORM016096_G3", "SJNORM016096_G4", "SJNORM016314_G4", "SJNORM016314_G5")

col_fun_rna_log2 <- colorRamp2(c(1,4,7), brewer.pal(9, "RdYlBu")[c(9,5,1)])

epi_reg <- c("DNMT1", "DNMT3A", "DNMT3B", "TET1", "TET2", "TET3", "MYC", "IDH1", "IDH2", "WT1", "EED", "EZH2", "SUZ12", "RING1", "RNF2", "KDM2B", "BAP1", "SUV39H1", "SUV39H2", "EHMT2", "ARID1A", "ARID1B", "ARID2", "PBRM1", "SMARCA4", "SMARCB1", "QSER1", "HELLS")

pdf("Extended_Data_Figure6b.pdf", width = 12, height = 8)
h1 <- Heatmap(as.matrix(log2(tpm_epi_reg[,order_control_samples_heatmap] + 1)), show_row_names = TRUE, top_annotation = columnAnnotation(df = metadata_rnaseq_with_avg_cgi_global[order_control_samples_heatmap,c("Type", "avg_cgi", "avg_global_no_cgi")], col = annotation_colors), row_labels = epi_reg, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_rna_log2)

h2 <- Heatmap(as.matrix(log2(tpm_epi_reg[,order_samples_heatmap] + 1)), show_row_names = TRUE, top_annotation = columnAnnotation(df = metadata_rnaseq_with_avg_cgi_global[order_samples_heatmap,c("Type", "avg_cgi", "avg_global_no_cgi")], col = annotation_colors), row_labels = epi_reg, cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_rna_log2)
draw(h1 + h2)
dev.off()

################################################################################
## Extended Data Figure 6d
################################################################################

mat_for_heatmap <- read.table("Extended_Data_Figure6d_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata_rnaseq_with_avg_cgi_global <- read.table("Figure6b_Extended_Data_Figure6b_6d_6e_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

annotation_colors <- list(Type = c("black", "grey35", brewer.pal(4, "Reds")[2:4], "deepskyblue1", "royalblue2", "#984EA3"))
names(annotation_colors[[1]]) <- c("control T cells", "control B cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

col_fun_cgi <- colorRamp2(c(0.1, 0.5), c("#F7FCFD", "#810F7C"))
col_fun_global <- colorRamp2(c(0.8, 1), c("#F7FCFD", "#006D2C"))

annotation_colors$avg_cgi <- col_fun_cgi
annotation_colors$avg_global_no_cgi <- col_fun_global

order_samples_heatmap <- c("SJBALL218_D", "SJBALL230_D", "SJERG016_D", "SJERG003_D", "SJERG011_D", "SJERG020299_D1", "SJERG000020_D", "SJERG000009_D", "SJERG014_D", "SJERG006_D", "SJERG000030_D", "SJERG020306_D1", "SJBALL020852_D1", "SJBALL021530_D1", "SJBALL021644_D1", "SJBALL015_D", "SJBALL021080_D1", "SJBALL021109_D1", "SJBALL000010_D", "SJTALL004_D", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012402_D1", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL187_D", "SJTALL022646_D1", "SJTALL186_D", "SJTALL012403_D1", "SJTALL022648_D1", "SJTALL021687_D1", "SJTALL022649_D1", "SJTALL170_D", "SJTALL021912_D1", "SJTALL015715_D1", "SJTALL021688_D1", "SJTALL012401_D1", "SJTALL182_D", "SJTALL021685_D1", "SJTALL021695_D1", "SJTALL169_D", "SJTALL021694_D1", "SJTALL022647_D1", "SJTALL009_D", "SJTALL022645_D1", "SJTALL174_D")
order_control_samples_heatmap <- c("SJNORM016096_G1", "SJNORM016096_G2", "SJNORM016096_G3", "SJNORM016096_G4", "SJNORM016314_G4", "SJNORM016314_G5")

col_fun_rna_log2 <- colorRamp2(c(1,2.5,4), brewer.pal(9, "RdYlBu")[c(9,5,1)])

pdf("Extended_Data_Figure6d.pdf", height = 5, width = 14)
h1 <- Heatmap(as.matrix(log2(mat_for_heatmap[,order_control_samples_heatmap] + 1)), show_row_names = TRUE, top_annotation = columnAnnotation(df = metadata_rnaseq_with_avg_cgi_global[order_control_samples_heatmap,c("Type", "avg_cgi", "avg_global_no_cgi")], col = annotation_colors), cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_rna_log2)

h2 <- Heatmap(as.matrix(log2(mat_for_heatmap[,order_samples_heatmap] + 1)), show_row_names = TRUE, top_annotation = columnAnnotation(df = metadata_rnaseq_with_avg_cgi_global[order_samples_heatmap,c("Type", "avg_cgi", "avg_global_no_cgi")], col = annotation_colors), cluster_rows = FALSE, cluster_columns = FALSE, col = col_fun_rna_log2)
draw(h1 + h2)
dev.off()

################################################################################
## Extended Data Figure 6e
################################################################################

tpm_epi_reg <- read.table("Extended_Data_Figure6b_6e_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
avg_promoter_epi_reg <- read.table("Figure6b_6c_Extended_Data_Figure_6e_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata_rnaseq_with_avg_cgi_global <- read.table("Figure6b_Extended_Data_Figure6b_6d_6e_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")
cor_cgi_dname_genes_with_significance <- read.table(file = "Figure6a_Extended_Data_Figure_6a_source_data_cor_test_cgi.tsv", sep = "\t", row.names = 1, header = TRUE)

common_samples <- c("SJBALL000010_D", "SJBALL015_D", "SJBALL020852_D1", "SJBALL021080_D1", "SJBALL021109_D1", "SJBALL021530_D1", "SJBALL021644_D1", "SJBALL218_D", "SJBALL230_D", "SJERG000009_D", "SJERG000020_D", "SJERG000030_D", "SJERG003_D", "SJERG006_D", "SJERG011_D", "SJERG014_D", "SJERG016_D", "SJERG020299_D1", "SJERG020306_D1", "SJNORM016096_G1", "SJNORM016096_G2", "SJNORM016096_G3", "SJNORM016096_G4", "SJNORM016314_G4", "SJNORM016314_G5", "SJTALL004_D", "SJTALL009_D", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012401_D1", "SJTALL012402_D1", "SJTALL012403_D1", "SJTALL015715_D1", "SJTALL021685_D1", "SJTALL021687_D1", "SJTALL021688_D1", "SJTALL021694_D1", "SJTALL021695_D1", "SJTALL021912_D1", "SJTALL022645_D1", "SJTALL022646_D1", "SJTALL022647_D1", "SJTALL022648_D1", "SJTALL022649_D1", "SJTALL169_D", "SJTALL170_D", "SJTALL174_D", "SJTALL182_D", "SJTALL186_D", "SJTALL187_D")

wt1_mutations_samples <- c("SJTALL012397_D1", "SJTALL012399_D1", "SJTALL012401_D1", "SJTALL021694_D1", "SJTALL022649_D1", "SJTALL170_D")

annotation_colors <- list(Type = c("black", "grey35", brewer.pal(4, "Reds")[2:4], "deepskyblue1", "royalblue2", "#984EA3"))
names(annotation_colors[[1]]) <- c("control T cells", "control B cells", "T-ALL LM", "T-ALL IM", "T-ALL HM", "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

pdf("Extended_Data_Figure6e.pdf")
for (i in c("TET1", "TET2", "WT1", "KDM2B"))
{
    gene_id <- rownames(subset(cor_cgi_dname_genes_with_significance, gene_name == i))

    if (i != "WT1")
    {
        plot(log2(do.call(c, tpm_epi_reg[gene_id,common_samples]) + 1), avg_promoter_epi_reg[i,common_samples], pch = 19, xlab = "Log2 (TPM + 1)", ylab = "Mean promoter methylation", col = annotation_colors$Type[as.character(metadata_rnaseq_with_avg_cgi_global[common_samples,"Type"])], main = i, ylim = c(0,1))
    } else
    {
        plot(log2(do.call(c, tpm_epi_reg[gene_id,common_samples]) + 1), avg_promoter_epi_reg[i,common_samples],xlab = "Log2 (TPM + 1)", ylab = "Mean promoter methylation", col = annotation_colors$Type[as.character(metadata_rnaseq_with_avg_cgi_global[common_samples,"Type"])], main = i, ylim = c(0,1), pch = ifelse(common_samples %in% wt1_mutations_samples, 17, 19))
    }
}
dev.off()

################################################################################
## Extended Data Figure 6f
################################################################################

## Data can be found in Extended_Data_Figure5c_6f_source_data.tsv
