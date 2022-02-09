################################################################################
## This script contains the code to produce Figure 5 and Extended Data Figure 5
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

################################################################################
## Figure 5a
################################################################################

df_oncoprint_tall <- read.table("Figure5a_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
annotation_column_tall_variants <- read.table("Figure5a_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

col <- c(brewer.pal(11, "Set3")[-9])
names(col) <- c("utr_3", "utr_5", "frameshift", "missense", "nonsense", "proteindel", "proteinins", "splice", "splice_region", "exon")

alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    utr_3 = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["utr_3"], col = NA))
    },
    utr_5 = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["utr_5"], col = NA))
    },
    frameshift = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["frameshift"], col = NA))
    },
    missense = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["missense"], col = NA))
    },
    nonsense = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["nonsense"], col = NA))
    },
    proteindel = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["proteindel"], col = NA))
    },
    proteinins = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["proteinins"], col = NA))
    },
    splice_region = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["splice_region"], col = NA))
    },
    exon = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["exon"], col = NA))
    },
    splice = function(x, y, w, h) {
        grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
            gp = gpar(fill = col["splice"], col = NA))
    }
)

heatmap_legend_param = list(title = "Alterations", at = c("exon", "frameshift", "missense", "nonsense", "proteindel", "proteinins", "splice", "splice_region", "utr_3", "utr_5")  , labels = c("exon", "frameshift", "missense", "nonsense", "proteindel", "proteinins", "splice", "splice_region", "utr_3", "utr_5")  )
sample_order_tall <- rownames(annotation_column_tall_variants)[order(annotation_column_tall_variants$median_variable_cgi)]

col_fun_median = colorRamp2(c(0, 0.5, 1), c("#0072B2", "seagreen3", "#F0E442"))

annotation_colors <- list(Type = c(brewer.pal(3, "Reds")), Sex = c("gold", "royalblue4"), Age_Group = c(brewer.pal(n = 3, name = "YlGn")))
names(annotation_colors[[1]]) <- c("T-ALL LM", "T-ALL IM", "T-ALL HM")
names(annotation_colors[[2]]) <- c("female", "male")
names(annotation_colors[[3]]) <- c("Pediatric", "AYA", "Adult")

annotation_colors$median_variable_cgi <- col_fun_median
annotation_colors$ALL.subtype <- brewer.pal(6, "Dark2")
names(annotation_colors$ALL.subtype) <- sort(unique(annotation_column_tall_variants$ALL.subtype))

pdf("Figure5a.pdf", width = 12, height = 9)
oncoPrint(df_oncoprint_tall, alter_fun = alter_fun, col = col, heatmap_legend_param = heatmap_legend_param, show_column_names = TRUE, column_order = sample_order_tall, bottom_annotation = HeatmapAnnotation(df = annotation_column_tall_variants, col = annotation_colors))
dev.off()

################################################################################
## Figure 5b
################################################################################

zscore_tpm_tall_top_500_genes_sd <- read.table("Figure5b_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
annotation_rna_clustering <- read.table("Figure5b_Extended_Data_Figure5b_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

col_fun_median = colorRamp2(c(0, 0.5, 1), c("#0072B2", "seagreen3", "#F0E442"))

annotation_colors <- list(Type = c(brewer.pal(3, "Reds")), Sex = c("gold", "royalblue4"), Age_Group = c(brewer.pal(n = 3, name = "YlGn")))
names(annotation_colors[[1]]) <- c("T-ALL LM", "T-ALL IM", "T-ALL HM")
names(annotation_colors[[2]]) <- c("female", "male")
names(annotation_colors[[3]]) <- c("Pediatric", "AYA", "Adult")

annotation_colors$median_variable_cgi <- col_fun_median
annotation_colors$ALL.subtype <- brewer.pal(6, "Dark2")
names(annotation_colors$ALL.subtype) <- sort(unique(annotation_rna_clustering$ALL.subtype))

annotation_colors$rna_cluster <- brewer.pal(3, "Accent")
names(annotation_colors$rna_cluster) <- c("1", "2", "3")
col_fun_rna <- colorRamp2(c(-3, 0, 3), c("#1B7837", "white", "#762A83"))

pdf("Figure5b.pdf", width = 8, height = 9)
Heatmap(as.matrix(zscore_tpm_tall_top_500_genes_sd), show_row_names = FALSE, top_annotation = columnAnnotation(df = annotation_rna_clustering, col = annotation_colors), col = col_fun_rna, clustering_method_columns = "ward.D2", clustering_method_rows = "ward.D2", clustering_distance_rows = "pearson", clustering_distance_columns = "pearson")
dev.off()

################################################################################
## Extended Data Figure 5a
################################################################################

avg_cgi_common_tall_variable_tall_samples <- read.table("Extended_Data_Figure5a_source_data_heatmap.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata_tall_clustering <- read.table("Extended_Data_Figure5a_source_data_heatmap_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")
avg_cgi_avg_subtypes_tall_variable <- read.table("Extended_Data_Figure5a_source_data_violin_plot.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

dist_tall_samples_cgis <- dist(t(avg_cgi_common_tall_variable_tall_samples))
breaksList <- seq(0, 62, by = 0.1)

annotation_colors <- list(Type = c(brewer.pal(3, "Reds")), Sex = c("gold", "royalblue4"), Age_Group = c(brewer.pal(n = 3, name = "YlGn")))
names(annotation_colors[[1]]) <- c("T-ALL LM", "T-ALL IM", "T-ALL HM")
names(annotation_colors[[2]]) <- c("female", "male")
names(annotation_colors[[3]]) <- c("Pediatric", "AYA", "Adult")

annotation_colors$ALL.subtype <- brewer.pal(6, "Dark2")
names(annotation_colors$ALL.subtype) <- sort(unique(metadata_tall_clustering$ALL.subtype))

samplesTALL <- c("SJALL030551_D1", "SJALL030933_D1", "SJALL048427_D3", "SJALL048442_D2", "SJALL064568_D1", "SJALL064571_D1", "SJALL064572_D1", "SJALL064573_D1", "SJTALL004_D", "SJTALL010_D", "SJTALL012397_D1", "SJTALL012398_D1", "SJTALL012399_D1", "SJTALL012400_D1", "SJTALL012402_D1", "SJTALL022646_D1", "SJTALL187_D", "SJALL030493_D1", "SJALL030684_D1", "SJALL048428_D4", "SJALL048440_D2", "SJALL048449_D2", "SJALL064567_D1", "SJALL064569_D1", "SJTALL012401_D1", "SJTALL012403_D1", "SJTALL015715_D1", "SJTALL021687_D1", "SJTALL021688_D1", "SJTALL021912_D1", "SJTALL022648_D1", "SJTALL022649_D1", "SJTALL170_D", "SJTALL182_D", "SJTALL186_D", "SJALL030916_D1", "SJALL059529_D2", "SJALL064564_D1", "SJALL064565_D1", "SJALL064566_D1", "SJTALL009_D", "SJTALL021685_D1", "SJTALL021694_D1", "SJTALL021695_D1", "SJTALL022645_D1", "SJTALL022647_D1", "SJTALL169_D", "SJTALL174_D")

pdf("Extended_Data_Figure_5a_heatmap.pdf", width = 17, height = 14)
pheatmap::pheatmap(dist_tall_samples_cgis, show_colnames = FALSE, labels_row = samplesTALL, annotation_row = metadata_tall_clustering, annotation_colors = annotation_colors, breaks = breaksList, color = colorRampPalette(brewer.pal(n = 11, name = "RdYlBu"))(length(breaksList)))
dev.off()

avg_cgi_avg_subtypes_tall_variable_df <- melt(avg_cgi_avg_subtypes_tall_variable)
avg_cgi_avg_subtypes_tall_variable_df$variable <- factor(avg_cgi_avg_subtypes_tall_variable_df$variable, levels = c("precursor T cells", "T-ALL LM", "T-ALL IM", "T-ALL HM"))

pdf("Extended_Data_Figure_5a_violin_plot.pdf", height = 5, width = 5)
par(mar=c(9, 4, 4, 2))
vioplot(value ~ variable, data = avg_cgi_avg_subtypes_tall_variable_df, ylim = c(0,1), las = 2, xlab = "", ylab = "Mean methylation variable CGIs", col = c("gray50", brewer.pal(4, "Reds")[-1]))
dev.off()

################################################################################
## Extended Data Figure 5b
################################################################################

zscore_tpm_tall_marker <- read.table("Extended_Data_Figure5b_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
annotation_rna_clustering <- read.table("Figure5b_Extended_Data_Figure5b_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

col_fun_median = colorRamp2(c(0, 0.5, 1), c("#0072B2", "seagreen3", "#F0E442"))

annotation_colors <- list(Type = c(brewer.pal(3, "Reds")), Sex = c("gold", "royalblue4"), Age_Group = c(brewer.pal(n = 3, name = "YlGn")))
names(annotation_colors[[1]]) <- c("T-ALL LM", "T-ALL IM", "T-ALL HM")
names(annotation_colors[[2]]) <- c("female", "male")
names(annotation_colors[[3]]) <- c("Pediatric", "AYA", "Adult")

annotation_colors$median_variable_cgi <- col_fun_median
annotation_colors$ALL.subtype <- brewer.pal(6, "Dark2")
names(annotation_colors$ALL.subtype) <- sort(unique(annotation_rna_clustering$ALL.subtype))

annotation_colors$rna_cluster <- brewer.pal(3, "Accent")
names(annotation_colors$rna_cluster) <- c("1", "2", "3")
col_fun_rna <- colorRamp2(c(-3, 0, 3), c("#1B7837", "white", "#762A83"))

tall_marker_genes <- c("TAL1", "TAL2", "LMO1", "LMO2", "LYL1", "TLX1", "TLX2", "TLX3", "BCL11B", "HOXA1", "HOXA2", "HOXA3", "HOXA4", "HOXA5", "HOXA6", "HOXA7", "HOXA9", "HOXA10", "HOXA11", "HOXA13", "NKX2-1", "CD1A", "CD1B", "CD1C", "CD1D", "CD1E", "CD4", "CD5", "CD8A", "CD8B", "ANPEP", "IL2RA", "CD33", "CD34", "KIT")
tall_marker_genes_labels <- gsub("ANPEP", "CD13", tall_marker_genes)
tall_marker_genes_labels <- gsub("KIT", "CD117", tall_marker_genes_labels)
tall_marker_genes_labels <- gsub("IL2RA", "CD25", tall_marker_genes_labels)

pdf("Extended_Data_Figure_5b.pdf", width = 8, height = 10)
Heatmap(as.matrix(zscore_tpm_tall_marker), show_row_names = TRUE, top_annotation = columnAnnotation(df = annotation_rna_clustering[colnames(zscore_tpm_tall_marker),c("Age_Group", "Sex", "Type", "ALL.subtype", "median_variable_cgi", "rna_cluster")], col = annotation_colors), col = col_fun_rna, row_labels = tall_marker_genes_labels, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

################################################################################
## Extended Data Figure 5c
################################################################################

## Data can be found in Extended_Data_Figure5c_6f_source_data.tsv

################################################################################
## Extended Data Figure 5d
################################################################################

avg_cgi_entropy_variable <- read.table("Extended_Data_Figure5d_source_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.table("Extended_Data_Figure5d_source_data_annotation.tsv", header = TRUE, row.names = 1, sep = "\t")

avg_cgi_entropy_dname_variable <- avg_cgi_entropy_variable[,grep("methylation", colnames(avg_cgi_entropy_variable))]
colnames(avg_cgi_entropy_dname_variable) <- gsub("methylation_", "", colnames(avg_cgi_entropy_dname_variable))

avg_cgi_entropy_variable <- avg_cgi_entropy_variable[,grep("entropy", colnames(avg_cgi_entropy_variable))]
colnames(avg_cgi_entropy_variable) <- gsub("entropy_", "", colnames(avg_cgi_entropy_variable))

median_cgi_per_tall_sample <- data.frame(row.names = colnames(avg_cgi_entropy_dname_variable), median_val = apply(avg_cgi_entropy_dname_variable, 2, median, na.rm = TRUE), variable = colnames(avg_cgi_entropy_dname_variable))
median_cgi_per_tall_sample$group <- "DNAme"

avg_cgi_entropy_variable_df <- melt(avg_cgi_entropy_variable)
avg_cgi_entropy_variable_df$variable <- factor(avg_cgi_entropy_variable_df$variable, levels = colnames(avg_cgi_entropy_dname_variable)[order(median_cgi_per_tall_sample[,1])])
avg_cgi_entropy_variable_df$type_broad <- ifelse(grepl("NORM", as.character(avg_cgi_entropy_variable_df$variable)), "control T cells", "T-ALL")
avg_cgi_entropy_variable_df$type <- metadata[as.character(avg_cgi_entropy_variable_df$variable),"Type"]

annotation_colors <- list(Type = c("grey", brewer.pal(3, "Reds")))
names(annotation_colors[[1]]) <- c("control T cells", "T-ALL LM", "T-ALL IM", "T-ALL HM")

pdf("Extended_Data_Figure_5d.pdf", height = 7, width = 18)
ggplot(avg_cgi_entropy_variable_df, aes(x=variable, y=value)) + geom_boxplot(outlier.shape = NA, alpha = 0.5, aes(fill = type)) + theme_classic()  + scale_fill_manual(values = annotation_colors$Type) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ylab("Mean entropy variable CGIs") + geom_point(data = median_cgi_per_tall_sample, aes(x = variable, y = median_val)) + geom_line(data = median_cgi_per_tall_sample, aes(x = variable, y = median_val, group = group), lty = 2)
dev.off()

avg_cgi_entropy_variable_mean <- apply(avg_cgi_entropy_variable, 2, mean, na.rm = TRUE)

wilcox.test(avg_cgi_entropy_variable_mean[rownames(subset(metadata, Type == "T-ALL LM"))], avg_cgi_entropy_variable_mean[rownames(subset(metadata, Type == "T-ALL HM"))]) ## 2.321e-06
wilcox.test(avg_cgi_entropy_variable_mean[rownames(subset(metadata, Type == "T-ALL IM"))], avg_cgi_entropy_variable_mean[rownames(subset(metadata, Type == "T-ALL HM"))]) ## 0.0001591
wilcox.test(avg_cgi_entropy_variable_mean[rownames(subset(metadata, Type == "T-ALL LM"))], avg_cgi_entropy_variable_mean[rownames(subset(metadata, Type == "T-ALL IM"))]) ## 0.1626
