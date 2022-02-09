################################################################################
## This script contains the code to produce Figure 1 and Extended Data Figure 1
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
## Figure 1b
################################################################################

## Figure 1b was generated using the processed methylation rates of samples
## MBC_S017RE, 182CLL, SJNORM016314_G1, SJTALL186_D, COAD_N3158, COAD_3158
## visualized in IGV

################################################################################
## Figure 1c
################################################################################

example_tracks <- data.frame(fread("Figure1c_Extended_Data_Figure1b_source_data.tsv", header = TRUE, sep = "\t"), stringsAsFactors = FALSE, check.names = FALSE)

pdf("Figure1c_densityplot.pdf", width = 15, height = 5)
layout(matrix(1:3, ncol = 3, byrow = TRUE))
smoothScatter(example_tracks[,"memory B cells"], example_tracks[,"CLL"], xlab = "Memory B cells", ylab = "CLL", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
smoothScatter(example_tracks[,"precursor T cells"], example_tracks[,"T-ALL"], xlab = "Precursor T cells", ylab = "T-ALL", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
smoothScatter(example_tracks[,"colon"], example_tracks[,"COAD"], xlab = "Colon", ylab = "COAD", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
dev.off()

diff_per_cpg <- data.frame(cll = example_tracks[,"CLL"] - example_tracks[,"memory B cells"],
                           tall = example_tracks[,"T-ALL"] - example_tracks[,"precursor T cells"],
                           coad = example_tracks[,"COAD"] - example_tracks[,"colon"])

count_per_cpg_class <- matrix(NA, nrow = 3, ncol = 3)
colnames(count_per_cpg_class) <- c("cll", "tall", "coad")
rownames(count_per_cpg_class) <- c("hyper", "stable", "hypo")
count_per_cpg_class["hyper",] <- apply(diff_per_cpg, 2, function(x) length(which(x > 0.1)))
count_per_cpg_class["hypo",] <- apply(diff_per_cpg, 2, function(x) length(which(x < -0.1)))
count_per_cpg_class["stable",] <- apply(diff_per_cpg, 2, function(x) length(which(x >= -0.1 & x <= 0.1)))
count_per_cpg_class <- data.frame(count_per_cpg_class)
count_per_cpg_class$feature <- rownames(count_per_cpg_class)
count_per_cpg_class_df <- melt(count_per_cpg_class)
count_per_cpg_class_df$feature <- factor(count_per_cpg_class_df$feature, levels = c("hypo", "stable", "hyper"))

pdf("Figure1c_barplot.pdf", width = 6, height = 6)
ggplot(data = count_per_cpg_class_df, aes(x = variable, fill = feature, y = value)) + geom_bar(position = "fill", stat = "identity") + theme_classic() + xlab("") + ylab("Fraction CpGs") + theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=12)) + scale_fill_manual(values = c("darkblue", "grey", "darkred"))
dev.off()

################################################################################
## Figure 1d
################################################################################

df_avg_methylation_per_sample <- read.table("Figure1d_source_data.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")

df_avg_methylation_per_sample$type <- factor(df_avg_methylation_per_sample$type, levels = c("control_hematopoietic", "T-ALL", "DUX4/ERG", "hypo-diploid", "Ph-like ALL", "B-ALL Blueprint", "AML", "TPLL", "CLL", "MCL", "control_solid", "BRCA", "UCEC", "LUSC", "LUAD", "COAD", "STAD", "READ", "BLCA"))
df_avg_methylation_per_sample$origin <- factor(df_avg_methylation_per_sample$origin, levels = c("This study", "Blueprint", "TCGA"))

pdf("Figure1d_pancancer.pdf", width = 12, height = 6)
ggplot(df_avg_methylation_per_sample, aes(x=type, y=avg_global_no_cgi, fill = type)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1)) + theme_classic() + scale_fill_manual(values = c("grey35", "firebrick", "deepskyblue1", "royalblue2", "#984EA3", brewer.pal(6, "Blues")[2:6], "grey35", brewer.pal(8, "Set3"))) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ylab("CGI mean methylation") + geom_jitter(aes(shape=origin), position=position_jitter(width = 0.3, height = 0, seed = 42), alpha = 0.6)
dev.off()

## Separate boxplot hematopoietic controls
df_avg_methylation_per_sample_hematopoietic_controls <- subset(df_avg_methylation_per_sample, type == "control_hematopoietic")
df_avg_methylation_per_sample_hematopoietic_controls$type <- c("B_cell_progenitor_CD34p_CD19p", "pro_B_cell_CD34n_CD19p_CD10p", "new_B_cell_CD34n_CD19p_CD10n_KLpos", "pre_B_cell_CD34n_CD19p_CD10n_KLneg", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p", "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p" , "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "CD4p_alpha_beta_Tcell", "CD4p_alpha_beta_Tcell", "CD4p_alpha_beta_Tcell", "CD4p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "HPC_CD34p_CD19n", "MBC", "MBC", "NBC", "NBC", "NBC", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell")
df_avg_methylation_per_sample_hematopoietic_controls$type_broad <- "control_hematopoietic"

df_avg_methylation_per_sample_hematopoietic_controls$type <- factor(df_avg_methylation_per_sample_hematopoietic_controls$type, levels = c("HPC_CD34p_CD19n", "B_cell_progenitor_CD34p_CD19p", "pro_B_cell_CD34n_CD19p_CD10p", "new_B_cell_CD34n_CD19p_CD10n_KLpos", "pre_B_cell_CD34n_CD19p_CD10n_KLneg", "NBC", "MBC", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p", "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "CD4p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell"))

colors_controls_hematopoietic <- c("grey40", brewer.pal(7, "PuRd")[2:7], brewer.pal(9, "YlGnBu"))
names(colors_controls_hematopoietic) <- c("HPC_CD34p_CD19n", "B_cell_progenitor_CD34p_CD19p", "pro_B_cell_CD34n_CD19p_CD10p", "new_B_cell_CD34n_CD19p_CD10n_KLpos", "pre_B_cell_CD34n_CD19p_CD10n_KLneg", "NBC", "MBC", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p", "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "CD4p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell")

pdf("Figure1d_controls.pdf", width = 5, height = 6)
ggplot(df_avg_methylation_per_sample_hematopoietic_controls, aes(x=type_broad, y=avg_global_no_cgi)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0.8, 0.95)) + theme_classic() + scale_fill_manual(values = c("grey35")) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ylab("CGI mean methylation") + geom_jitter(aes(shape=origin, color = type), position=position_jitter(width = 0.3, height = 0, seed = 42)) + scale_color_manual(values = colors_controls_hematopoietic)
dev.off()

################################################################################
## Extended Data Figure 1a
################################################################################

## Extended Data Figure 1a was generated using the processed methylation rates of samples
## SJNORM016096_G1, SJERG003_D, SJHYPO003046_D, SJBALL020488_D1
## visualized in IGV

################################################################################
## Extended Data Figure 1b
################################################################################

example_tracks <- data.frame(fread("Figure1c_Extended_Data_Figure1b_source_data.tsv", header = TRUE, sep = "\t"), stringsAsFactors = FALSE, check.names = FALSE)

pdf("Extended_Data_Figure1b_densityplot.pdf", width = 15, height = 5)
layout(matrix(1:3, ncol = 3, byrow = TRUE))
smoothScatter(example_tracks[,"precursor B cells"], example_tracks[,"DUX4/ERG"], xlab = "Precursor B cells", ylab = "DUX4/ERG", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
smoothScatter(example_tracks[,"precursor B cells"], example_tracks[,"hypodiploid"], xlab = "Precursor B cells", ylab = "Hypodiploid", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
smoothScatter(example_tracks[,"precursor B cells"], example_tracks[,"Ph-like"], xlab = "Precursor B cells", ylab = "Ph-like", colramp = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu"))), cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0)
abline(0.1, 1, lty = 2)
abline(-0.1, 1, lty = 2)
dev.off()

diff_per_cpg <- data.frame(erg = example_tracks[,"DUX4/ERG"] - example_tracks[,"precursor B cells"],
                           hypo = example_tracks[,"hypodiploid"] - example_tracks[,"precursor B cells"],
                           phlike = example_tracks[,"Ph-like"] - example_tracks[,"precursor B cells"])

count_per_cpg_class <- matrix(NA, nrow = 3, ncol = 3)
colnames(count_per_cpg_class) <- c("erg", "hypo", "phlike")
rownames(count_per_cpg_class) <- c("hyper", "stable", "hypo")
count_per_cpg_class["hyper",] <- apply(diff_per_cpg, 2, function(x) length(which(x > 0.1)))
count_per_cpg_class["hypo",] <- apply(diff_per_cpg, 2, function(x) length(which(x < -0.1)))
count_per_cpg_class["stable",] <- apply(diff_per_cpg, 2, function(x) length(which(x >= -0.1 & x <= 0.1)))
count_per_cpg_class <- data.frame(count_per_cpg_class)
count_per_cpg_class$feature <- rownames(count_per_cpg_class)
count_per_cpg_class_df <- melt(count_per_cpg_class)
count_per_cpg_class_df$feature <- factor(count_per_cpg_class_df$feature, levels = c("hypo", "stable", "hyper"))

pdf("Extended_Data_Figure1b_barplot.pdf", width = 6, height = 6)
ggplot(data = count_per_cpg_class_df, aes(x = variable, fill = feature, y = value)) + geom_bar(position = "fill", stat = "identity") + theme_classic() + xlab("") + ylab("Fraction CpGs") + theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=12)) + scale_fill_manual(values = c("darkblue", "grey", "darkred"))
dev.off()

################################################################################
## Extended Data Figure 1c
################################################################################

df_avg_methylation_per_sample_tall <- read.table("Extended_Data_Figure1c_source_data.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
df_avg_methylation_per_sample_tall$Age_Group[is.na(df_avg_methylation_per_sample_tall$Age_Group)] <- "Unknown"
df_avg_methylation_per_sample_tall$Age_Group <- factor(df_avg_methylation_per_sample_tall$Age_Group, levels = c("Pediatric", "AYA", "Adult", "Unknown"))

pdf("Extended_Data_Figure1c.pdf", width = 4, height = 6)
ggplot(df_avg_methylation_per_sample_tall, aes(x=Age_Group, y=avg_global_no_cgi, fill = Age_Group)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1)) + theme_classic() + scale_fill_manual(values = c(brewer.pal(n = 3, name = "YlGn"), "grey")) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ylab("Global no CGI mean methylation") + geom_jitter(shape=20, position=position_jitter(width = 0.3, height = 0, seed = 42), alpha = 0.6)
ggplot(df_avg_methylation_per_sample_tall, aes(x=Sex, y=avg_global_no_cgi, fill = Sex)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1)) + theme_classic() + scale_fill_manual(values = c("gold", "royalblue4")) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ylab("Global no CGI mean methylation") + geom_jitter(shape=20, position=position_jitter(width = 0.3, height = 0, seed = 42), alpha = 0.6)
dev.off()

wilcox.test(subset(df_avg_methylation_per_sample_tall, Age_Group == "Pediatric")$avg_global_no_cgi, subset(df_avg_methylation_per_sample_tall, Age_Group %in% c("AYA", "Adult"))$avg_global_no_cgi)
wilcox.test(subset(df_avg_methylation_per_sample_tall, Sex == "female")$avg_global_no_cgi, subset(df_avg_methylation_per_sample_tall, Sex == "male")$avg_global_no_cgi)
