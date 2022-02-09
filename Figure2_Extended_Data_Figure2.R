################################################################################
## This script contains the code to produce Figure 2 and Extended Data Figure 2
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
## Figure 2a
################################################################################

diff_avg_sw_pmd_hmd <- data.frame(fread("Figure2a_Extended_Data_Figure2b_source_data.tsv", header = TRUE, sep = "\t"), stringsAsFactors = FALSE, row.names = 1)

diff_avg_sw_pmd_hmd_df <- melt(diff_avg_sw_pmd_hmd)
diff_avg_sw_pmd_hmd_df$variable <- factor(diff_avg_sw_pmd_hmd_df$variable, levels = c(rev(c("BRCA", "UCEC", "LUSC", "LUAD", "COAD", "STAD", "READ", "BLCA")), "MCL", "CLL","TPLL", "AML", "ball_blueprint", "ball_phlike", "ball_hypo", "ball_erg", "tall", "differentiation_control", "precursor_control"))
diff_avg_sw_pmd_hmd_df$feature <- factor(diff_avg_sw_pmd_hmd_df$feature, levels = c("PMD", "HMD"))

pdf("Figure2a.pdf", width = 7, height = 9)
ggplot(diff_avg_sw_pmd_hmd_df, aes(x=variable, y=value, fill = feature)) + theme_classic() + geom_boxplot(outlier.shape = NA, alpha = 0.5) + theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=14)) + ylab("Difference sliding window") + xlab("") + coord_flip(ylim=c(-0.75, 0.5)) + geom_hline(yintercept = 0.1, lty = 2) + geom_hline(yintercept = -0.1, lty = 2) + geom_hline(yintercept = 0, lty = 2) + scale_fill_manual(values = c("darkblue", "forestgreen"))
dev.off()

################################################################################
## Figure 2b
################################################################################

mean_solo_wcgw_pancancer_chr16 <- read.table(file = "Figure2b_source_data.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
annotation_df_regions_solo_wcgw_pancancer_chr16 <- read.table(file = "Figure2b_source_data_annotation_regions.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
annotation_df_samples_solo_wcgw_pancancer_chr16 <- read.table(file = "Figure2b_source_data_annotation_samples.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")

annotation_colors <- NULL
annotation_colors$Region <- c("PMD" = "darkblue", "HMD" = "forestgreen")
annotation_colors$Type <- c("grey40", brewer.pal(4, "YlGnBu")[2:4], brewer.pal(4, "PuRd")[2:4], "grey60", "firebrick", "deepskyblue1", "royalblue2", "#984EA3", brewer.pal(6, "Blues")[2:6], brewer.pal(8, "Set3"))
names(annotation_colors$Type) <- c("hematopoietic progenitor", "precursor T cell", "T cell", "memory T cell", "precursor B cell", "naive B cell", "memory B cell", "control_solid", "T-ALL", "DUX4/ERG", "hypo-diploid", "Ph-like ALL", "B-ALL Blueprint", "AML", "TPLL", "CLL", "MCL", "BRCA", "UCEC", "LUSC", "LUAD", "COAD", "STAD", "READ", "BLCA")

breaksList <- seq(0, 1, by = 0.01)

pdf("Figure2b.pdf", width = 12)
pheatmap::pheatmap(mean_solo_wcgw_pancancer_chr16, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(length(breaksList)), breaks = breaksList, annotation_row = annotation_df_samples_solo_wcgw_pancancer_chr16,  annotation_colors = annotation_colors, show_rownames = FALSE, annotation_col = annotation_df_regions_solo_wcgw_pancancer_chr16, main = c("chr16p"), gaps_row = cumsum(table(annotation_df_samples_solo_wcgw_pancancer_chr16$condition)))
dev.off()

################################################################################
## Figure 2c
################################################################################

dmrs_overlap_combined_features_table_normalized <- read.table(file = "Figure2c_source_data.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
dmrs_overlap_combined_features_table_normalized$feature <- factor(dmrs_overlap_combined_features_table_normalized$feature, levels = c("CGI", "shore", "shelf", "DMV", "PMD", "TssActive", "TssBivalent", "Transcript", "Enhancer", "Repressive", "Heterochromatin", "Quiescent"))
dmrs_overlap_combined_features_table_normalized$subtype <- factor(dmrs_overlap_combined_features_table_normalized$subtype, levels = c("T-ALL", "ERG", "Hypo", "PhLike"))

pdf("Figure2c.pdf", width = 10)
ggplot(dmrs_overlap_combined_features_table_normalized, aes(x=feature, y=norm_pct, color=subtype)) + geom_point(shape = 4, size = 3, stroke = 1) + ylab("Enrichment of DMRs in features") + xlab ("") + theme_classic() + theme(axis.text.x=element_text(hjust = 1, angle = 45, size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=14)) + facet_grid(~type_dmr) + scale_color_manual(values = c("darkred", "deepskyblue1", "royalblue2", "#984EA3")) + geom_hline(yintercept = 1, lty = 2)
dev.off()

################################################################################
## Extended Data Figure 2a
################################################################################

avg_cgi_common_all_no_ccl_variable <- read.table("Extended_Data_Figure3b_source_data.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
metadata <- read.table("Extended_Data_Figure3b_source_data_annotation.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")

avg_cgi_common_all_no_ccl_variable_df <- melt(avg_cgi_common_all_no_ccl_variable)
avg_cgi_common_all_no_ccl_variable_df$type <- factor(metadata[as.character(avg_cgi_common_all_no_ccl_variable_df$variable),"type"], levels = c("control T cells", "control B cells", "T-ALL",     "DUX4/ERG", "hypo-diploid", "Ph-like ALL"))

type_for_vioplot <- data.frame(unique(avg_cgi_common_all_no_ccl_variable_df[c("variable", "type")]), row.names = 1)

avg_cgi_common_all_no_ccl_variable_df$variable <- factor(avg_cgi_common_all_no_ccl_variable_df$variable, levels = colnames(avg_cgi_common_all_no_ccl_variable)[order(type_for_vioplot[levels(avg_cgi_common_all_no_ccl_variable_df$variable),1], apply(avg_cgi_common_all_no_ccl_variable, 2, median))])

avg_tile_1kb_all_no_ccl <- data.frame(fread("Extended_Data_Figure2a_source_data.tsv", sep = "\t", header = TRUE), stringsAsFactors = FALSE, row.names = 1)
metadata <- read.table("Extended_Data_Figure2a_source_data_annotation.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")

avg_tile_1kb_all_no_ccl_df <- melt(avg_tile_1kb_all_no_ccl)
avg_tile_1kb_all_no_ccl_df$variable <- factor(avg_tile_1kb_all_no_ccl_df$variable, levels = levels(avg_cgi_common_all_no_ccl_variable_df$variable))

annotation_colors <- c("grey20", "grey35", "firebrick", "deepskyblue1", "royalblue2", "#984EA3")
names(annotation_colors) <- c("control T cells", "control B cells", "T-ALL",     "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

pdf("Extended_Data_Figure2a.pdf", height = 5, width = 55)
par(mar=c(9, 4, 4, 2))
vioplot(value ~ variable, data = subset(avg_tile_1kb_all_no_ccl_df, feature == "HMD"), ylim = c(0,1), col = annotation_colors[as.character(metadata[levels(avg_tile_1kb_all_no_ccl_df$variable),1])], las = 2, side = "left", plotCentre = "line", xlab = "", ylab = "Mean methylation 1kb tiles")
vioplot(value ~ variable, data = subset(avg_tile_1kb_all_no_ccl_df, feature == "PMD"), ylim = c(0,1), col = annotation_colors[as.character(metadata[levels(avg_tile_1kb_all_no_ccl_df$variable),1])], las = 2, side = "right", plotCentre = "line", add = TRUE)
dev.off()

################################################################################
## Extended Data Figure 2b
################################################################################

diff_avg_sw_pmd_hmd <- data.frame(fread("Figure2a_Extended_Data_Figure2b_source_data.tsv", header = TRUE, sep = "\t"), stringsAsFactors = FALSE, row.names = 1)

diff_avg_sw_pmd_hmd_hyper_freq <- data.frame(sapply(c("HMD", "PMD"), function(x) apply(subset(diff_avg_sw_pmd_hmd, feature == x)[,1:19], 2, function(y) length(which(y > 0.1)))))
diff_avg_sw_pmd_hmd_hypo_freq <- data.frame(sapply(c("HMD", "PMD"), function(x) apply(subset(diff_avg_sw_pmd_hmd, feature == x)[,1:19], 2, function(y) length(which(y < -0.1)))))

diff_avg_sw_pmd_hmd_hyper_freq$condition <- "hyper"
diff_avg_sw_pmd_hmd_hypo_freq$condition <- "hypo"

diff_avg_sw_pmd_hmd_hyper_freq$type <- rownames(diff_avg_sw_pmd_hmd_hyper_freq)
diff_avg_sw_pmd_hmd_hypo_freq$type <- rownames(diff_avg_sw_pmd_hmd_hypo_freq)

diff_avg_sw_pmd_hmd_hyper_freq$PMD <- diff_avg_sw_pmd_hmd_hyper_freq$PMD / apply(diff_avg_sw_pmd_hmd[,1:19], 2, function(x) length(na.omit(x)))
diff_avg_sw_pmd_hmd_hyper_freq$HMD <- diff_avg_sw_pmd_hmd_hyper_freq$HMD / apply(diff_avg_sw_pmd_hmd[,1:19], 2, function(x) length(na.omit(x)))
diff_avg_sw_pmd_hmd_hypo_freq$PMD <- diff_avg_sw_pmd_hmd_hypo_freq$PMD / apply(diff_avg_sw_pmd_hmd[,1:19], 2, function(x) length(na.omit(x)))
diff_avg_sw_pmd_hmd_hypo_freq$HMD <- diff_avg_sw_pmd_hmd_hypo_freq$HMD / apply(diff_avg_sw_pmd_hmd[,1:19], 2, function(x) length(na.omit(x)))

num_hyper_hypo_sw_all_subtypes <- rbind(diff_avg_sw_pmd_hmd_hyper_freq, diff_avg_sw_pmd_hmd_hypo_freq)
num_hyper_hypo_sw_all_subtypes_df <- melt(num_hyper_hypo_sw_all_subtypes)

num_hyper_hypo_sw_all_subtypes_df$type <- factor(num_hyper_hypo_sw_all_subtypes$type, levels = c(rev(c("BRCA", "UCEC", "LUSC", "LUAD", "COAD", "STAD", "READ", "BLCA")), "MCL", "CLL", "TPLL", "AML", "ball_blueprint", "ball_phlike", "ball_hypo", "ball_erg", "tall", "differentiation_control", "precursor_control"))
num_hyper_hypo_sw_all_subtypes_df[which(num_hyper_hypo_sw_all_subtypes_df$condition == "hypo"),"value"] <- -num_hyper_hypo_sw_all_subtypes_df[which(num_hyper_hypo_sw_all_subtypes_df$condition == "hypo"),"value"]
num_hyper_hypo_sw_all_subtypes_df$variable <- paste(num_hyper_hypo_sw_all_subtypes_df$condition, num_hyper_hypo_sw_all_subtypes_df$variable, sep = "_")

pdf("Extended_Data_Figure2b.pdf", width = 9)
ggplot(data = num_hyper_hypo_sw_all_subtypes_df, aes(x = type, y = value, fill = variable)) + geom_bar(position = "stack", stat = "identity") + theme_classic() + xlab("") + ylab("Fraction of sliding windows") + theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=12)) + coord_flip(ylim = c(-0.65, 0.25)) + geom_hline(yintercept = 0.1, lty = 2) + geom_hline(yintercept = -0.1, lty = 2) + geom_hline(yintercept = 0, lty = 2) + scale_fill_manual(values = c(brewer.pal(4, "Paired"))[c(3,1,4,2)])
dev.off()

################################################################################
## Extended Data Figure 2c
################################################################################

## Data can be found in Extended_Data_Figure2c_source_data.tsv
