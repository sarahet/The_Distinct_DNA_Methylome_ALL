################################################################################
## This script contains the code to produce Figure 3 and Extended Data Figure 3
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
## Figure 3a
################################################################################

df_avg_methylation_per_sample <- read.table("Figure3a_source_data.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")

df_avg_methylation_per_sample$type <- factor(df_avg_methylation_per_sample$type, levels = c("control_hematopoietic", "T-ALL", "DUX4/ERG", "hypo-diploid", "Ph-like ALL", "B-ALL Blueprint", "AML", "TPLL", "CLL", "MCL", "control_solid", "BRCA", "UCEC", "LUSC", "LUAD", "COAD", "STAD", "READ", "BLCA"))
df_avg_methylation_per_sample$origin <- factor(df_avg_methylation_per_sample$origin, levels = c("This study", "Blueprint", "TCGA"))

pdf("Figure3a_pancancer.pdf", width = 12, height = 6)
ggplot(df_avg_methylation_per_sample, aes(x=type, y=avg_cgi, fill = type)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0, 1)) + theme_classic() + scale_fill_manual(values = c("grey35", "firebrick", "deepskyblue1", "royalblue2", "#984EA3", brewer.pal(6, "Blues")[2:6], "grey35", brewer.pal(8, "Set3"))) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ylab("CGI mean methylation") + geom_jitter(aes(shape=origin), position=position_jitter(width = 0.3, height = 0, seed = 42), alpha = 0.6)
dev.off()

## Separate boxplot hematopoietic controls
df_avg_methylation_per_sample_hematopoietic_controls <- subset(df_avg_methylation_per_sample, type == "control_hematopoietic")
df_avg_methylation_per_sample_hematopoietic_controls$type <- c("B_cell_progenitor_CD34p_CD19p", "pro_B_cell_CD34n_CD19p_CD10p", "new_B_cell_CD34n_CD19p_CD10n_KLpos", "pre_B_cell_CD34n_CD19p_CD10n_KLneg", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p", "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p" , "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "CD4p_alpha_beta_Tcell", "CD4p_alpha_beta_Tcell", "CD4p_alpha_beta_Tcell", "CD4p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "HPC_CD34p_CD19n", "MBC", "MBC", "NBC", "NBC", "NBC", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell")
df_avg_methylation_per_sample_hematopoietic_controls$type_broad <- "control_hematopoietic"

df_avg_methylation_per_sample_hematopoietic_controls$type <- factor(df_avg_methylation_per_sample_hematopoietic_controls$type, levels = c("HPC_CD34p_CD19n", "B_cell_progenitor_CD34p_CD19p", "pro_B_cell_CD34n_CD19p_CD10p", "new_B_cell_CD34n_CD19p_CD10n_KLpos", "pre_B_cell_CD34n_CD19p_CD10n_KLneg", "NBC", "MBC", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p", "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "CD4p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell"))

colors_controls_hematopoietic <- c("grey40", brewer.pal(7, "PuRd")[2:7], brewer.pal(9, "YlGnBu"))
names(colors_controls_hematopoietic) <- c("HPC_CD34p_CD19n", "B_cell_progenitor_CD34p_CD19p", "pro_B_cell_CD34n_CD19p_CD10p", "new_B_cell_CD34n_CD19p_CD10n_KLpos", "pre_B_cell_CD34n_CD19p_CD10n_KLneg", "NBC", "MBC", "ETP_CD34p_CD1an_CD3n_CD4n_CD8n", "cortical_CD34n_CD1ap_CD4p_CD8p", "CD1ap_CD3p_CD4_CD8_double_positive", "CD1ap_CD3p_CD4_CD8_single_positive", "CD1an_CD3p_CD4_CD8_single_positive", "CD4p_alpha_beta_Tcell", "CD8p_alpha_beta_Tcell", "central_memory_CD4p_alpha_beta_Tcell", "central_memory_CD8p_alpha_beta_Tcell")

pdf("Figure3a_controls.pdf", width = 5, height = 6)
ggplot(df_avg_methylation_per_sample_hematopoietic_controls, aes(x=type_broad, y=avg_cgi)) + geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(0.15, 0.20)) + theme_classic() + scale_fill_manual(values = c("grey35")) + theme(axis.text=element_text(size=10), axis.title=element_text(size=12), axis.text.x=element_text(angle=45, hjust=1)) + xlab("") + ylab("CGI mean methylation") + geom_jitter(aes(shape=origin, color = type), position=position_jitter(width = 0.3, height = 0, seed = 42)) + scale_color_manual(values = colors_controls_hematopoietic)
dev.off()

################################################################################
## Figure 3b
################################################################################

## Figure 3b was generated using the processed methylation rates of samples
## MBC_S017RE, 182CLL, SJNORM016314_G1, SJTALL186_D, COAD_N3158, COAD_3158
## visualized in IGV

################################################################################
## Figure 3c
################################################################################

df_avg_methylation_per_sample <- read.table("Figure3c_source_data.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
df_avg_methylation_per_sample$origin <- factor(df_avg_methylation_per_sample$origin, levels = c("This study", "Blueprint"))
df_avg_methylation_per_sample$type <- factor(df_avg_methylation_per_sample$type, levels = c("control_hematopoietic", "T-ALL", "DUX4/ERG", "hypo-diploid", "Ph-like ALL", "B-ALL Blueprint"))

pdf("Figure3c.pdf", height = 6, width = 8)
ggplot(df_avg_methylation_per_sample, aes(x=avg_global_no_cgi, y=avg_cgi, group = 1, color = type, shape = origin)) + geom_point() + theme_classic() + scale_color_manual(values = c("grey35", "firebrick", "deepskyblue1", "royalblue2", "#984EA3", brewer.pal(6, "Blues")[2], "grey75")) + coord_cartesian(ylim = c(0.1, 0.5), xlim = c(0.75, 0.95)) + xlab("Global mean methylation (outside of CGIs)") + ylab("CGI mean methylation")
dev.off()

################################################################################
## Extended Data Figure 3a
################################################################################

## Extended Data Figure 3a was generated using the processed methylation rates of samples
## SJNORM016096_G1, SJERG003_D, SJHYPO003046_D, SJBALL020488_D1
## visualized in IGV

################################################################################
## Extended Data Figure 3b
################################################################################

avg_cgi_common_all_no_ccl_variable <- read.table("Extended_Data_Figure3b_source_data.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
metadata <- read.table("Extended_Data_Figure3b_source_data_annotation.tsv", header = TRUE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")

avg_cgi_common_all_no_ccl_variable_df <- melt(avg_cgi_common_all_no_ccl_variable)
avg_cgi_common_all_no_ccl_variable_df$type <- factor(metadata[as.character(avg_cgi_common_all_no_ccl_variable_df$variable),"type"], levels = c("control T cells", "control B cells", "T-ALL",     "DUX4/ERG", "hypo-diploid", "Ph-like ALL"))

type_for_vioplot <- data.frame(unique(avg_cgi_common_all_no_ccl_variable_df[c("variable", "type")]), row.names = 1)

avg_cgi_common_all_no_ccl_variable_df$variable <- factor(avg_cgi_common_all_no_ccl_variable_df$variable, levels = colnames(avg_cgi_common_all_no_ccl_variable)[order(type_for_vioplot[levels(avg_cgi_common_all_no_ccl_variable_df$variable),1], apply(avg_cgi_common_all_no_ccl_variable, 2, median))])

annotation_colors <- c("black", "grey35", "firebrick", "deepskyblue1", "royalblue2", "#984EA3")
names(annotation_colors) <- c("control T cells", "control B cells", "T-ALL", "DUX4/ERG", "hypo-diploid", "Ph-like ALL")

pdf("Extended_Data_Figure3b.pdf", height = 5, width = 55)
par(mar=c(9, 4, 4, 2))
vioplot(value ~ variable, data = avg_cgi_common_all_no_ccl_variable_df, ylim = c(0,1), col = annotation_colors[as.character(type_for_vioplot[levels(avg_cgi_common_all_no_ccl_variable_df$variable),1])], las = 2, drawRect = FALSE, xlab = "", ylab = "Mean methylation variable CGIs")
points(1:length(levels(avg_cgi_common_all_no_ccl_variable_df$variable)), apply(avg_cgi_common_all_no_ccl_variable[,levels(avg_cgi_common_all_no_ccl_variable_df$variable)], 2, median), pch = 19, col = "white")
dev.off()
