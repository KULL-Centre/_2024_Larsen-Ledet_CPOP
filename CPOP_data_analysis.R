#Load in required packages
library(tidyverse)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)
library(pROC)
library(forstringr)
library(ggpubr)
library(ggpmisc)

#Change the path to your directory
path <- "your_path"

#Read in dataframes
cpop_data <- read.csv(file.path(path, "cpop_data.csv"))
cpop_data_ROC_ins <- read.csv(file.path(path, "cpop_data_ROC_ins.csv"))
cpop_data_ROC_del <- read.csv(file.path(path, "cpop_data_ROC_del.csv")) 
dssp <- read.csv(file.path(path, "1u72_dssp.csv"), sep = ",")
ins_dplddt_ddg <- read.csv(file.path(path, "ins_dplddt_ddg.csv"))
del_dplddt_ddg <- read.csv(file.path(path, "del_dplddt_ddg.csv"))

#Create dataframes for each type of mutation
cpop_ins_gather <- cpop_data %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_ins_syn_gather <- cpop_data %>%
  filter(mut_type == "syn") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_ins_stop_gather <- cpop_data %>%
  filter(mut_type == "stop") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_del_gather <- cpop_data %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_del_syn_gather <- cpop_data %>%
  filter(mut_type == "syn") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_del_stop_gather <- cpop_data %>%
  filter(mut_type == "stop") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))


#Fig. 3A

#CPOP score distributions for insertions
ggplot() + 
  geom_histogram(cpop_ins_syn_gather, mapping = aes(score, fill = "syn", color = "syn"), binwidth = 0.035, alpha = 1) +
  geom_histogram(cpop_ins_stop_gather, mapping = aes(score, fill = "stop", color = "stop"), binwidth = 0.035, alpha = 1) +
  geom_histogram(cpop_ins_gather, mapping = aes(score, fill = "ins", color = "ins"), binwidth = 0.035, alpha = 0.8) +
  geom_segment(cpop_ins_gather, mapping = aes(x = 1, y = 0, xend = 1, yend = 113), color = "black", linetype = "longdash") +
  annotate("label", x = 1, y = 117, label = "WT" , color = "black", alpha = 0.3) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100)) +
  scale_colour_manual(name = NULL, values = c("ins" = "black", "syn"="black", "stop" = "black"), labels = c("ins" = "Insertion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("ins","syn", "stop")) +
  scale_fill_manual(name = NULL, values = adjustcolor(c("darkorange","white", "darkred"), alpha.f = 0.3), labels = c("ins" = "Insertion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("ins","syn", "stop")) +
  guides(fill = guide_legend(override.aes = list(alpha = c(0.8, 1, 1)))) +
  xlab("CPOP score") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14), panel.spacing = unit(1, "lines"), text = element_text(size = 14))

#CPOP score distributions for deletions
ggplot() + 
  geom_histogram(cpop_del_syn_gather, mapping = aes(score, fill = "syn", color = "syn"), binwidth = 0.035, alpha = 1) +
  geom_histogram(cpop_del_stop_gather, mapping = aes(score, fill = "stop", color = "stop"), binwidth = 0.035, alpha = 1) +
  geom_histogram(cpop_del_gather, mapping = aes(score, fill = "del", color = "del"), binwidth = 0.035, alpha = 0.8) +
  geom_segment(cpop_del_gather, mapping = aes(x = 1, y = 0, xend = 1, yend = 113), color = "black", linetype = "longdash") +
  annotate("label", x = 1, y = 117, label = "WT" , color = "black", alpha = 0.3) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(limits = c(0, 120), breaks = c(0, 25, 50, 75, 100)) +
  scale_colour_manual(name = NULL, values = c("del" = "black", "syn"="black", "stop" = "black"), labels = c("del" = "Deletion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("del","syn", "stop")) +
  scale_fill_manual(name = NULL, values = adjustcolor(c("darkorange","white", "darkred"), alpha.f = 0.3), labels = c("del" = "Deletion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("del","syn", "stop")) +
  guides(fill = guide_legend(override.aes = list(alpha = c(0.8, 1, 1)))) +
  xlab("CPOP score") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14), panel.spacing = unit(1, "lines"), text = element_text(size = 14))


#Fig. 3B

#CPOP heatmaps for insertions
cpop_ins <- cpop_data %>%
  filter(mut_type == "ins") %>%
  select(pos, score_30, score_30_MTX, score_37, score_37_MTX)
cpop_ins_compl <- t(cpop_ins)
cpop_ins_compl_heatmap <- cpop_ins_compl[-1,]
colnames(cpop_ins_compl_heatmap) <- cpop_ins_compl[1,]
draw(Heatmap(cpop_ins_compl_heatmap,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(-0.228, 1 , 2.228), c("darkred", "white", "darkblue")),
             name = ("CPOP score"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(0, 1, 2),
                                         labels = c("0", "1", "2")),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1.5),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                          legend_gp = gpar(fill = c("dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 12),
                                          row_gap = unit(2, "mm"))))


#CPOP heatmaps for deletions
cpop_del <- cpop_data %>%
  filter(mut_type == "del") %>%
  select(pos, score_30, score_30_MTX, score_37, score_37_MTX)
cpop_del_compl <- t(cpop_del)  
cpop_del_compl_heatmap <- cpop_del_compl[-1,]
colnames(cpop_del_compl_heatmap) <- cpop_del_compl[1,]
draw(Heatmap(cpop_del_compl_heatmap,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(-0.288, 1 , 2.288), c("darkred", "white", "darkblue")),
             name = ("CPOP score"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(0, 1, 2),
                                         labels = c("0", "1", "2")),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1.5),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                          legend_gp = gpar(fill = c("dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 12),
                                          row_gap = unit(2, "mm"))))


#Fig. 4B

#Group secondary structure elements
dssp_group <- dssp %>%
  mutate(pos = DSSP.index + 1,
         sec = ifelse(Secondary.structure == "G" | 
                      Secondary.structure == "H" | 
                      Secondary.structure == "I", "helices",
               ifelse(Secondary.structure == "E" | 
                      Secondary.structure == "B", "sheets",
               ifelse(Secondary.structure == "T" |
                      Secondary.structure == "S" |
                      Secondary.structure == "-", "loops", NA))))

#Dot plots for secondary structure elements for insertions
dssp_cpop_ins <- merge(cpop_data, dssp_group, by = "pos") %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
dssp_cpop_ins$sec <- factor(dssp_cpop_ins$sec, levels = c("loops", "helices", "sheets"))
ggplot(dssp_cpop_ins, aes(experiment, score, fill = sec), na.rm = F) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = "dodge", dotsize = 6, binwidth = 0.0078, stackratio = 0.9) +
  scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                    labels=c("Loops", "Helices", "Sheets"),
                    name = "Secondary structure") +
  scale_x_discrete(labels = c("30 °C", "30 °C + MTX", "37 °C", "37 °C + MTX")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0)

#Dot plots for secondary structure elements for deletions
dssp_cpop_del <- merge(cpop_data, dssp_group, by = "pos") %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
n_fun <- function(x){
  return(data.frame(y = 1.2, label = length(x)))
}
dssp_cpop_del$sec <- factor(dssp_cpop_del$sec, levels = c("loops", "helices", "sheets"))
ggplot(dssp_cpop_del, aes(experiment, score, fill = sec), na.rm = F) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', position = "dodge", dotsize = 6, binwidth = 0.0078, stackratio = 0.9) +
  scale_fill_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                    labels=c("Loops", "Helices", "Sheets"),
                    name = "Secondary structure") +
  scale_x_discrete(labels = c("30 °C", "30 °C + MTX", "37 °C", "37 °C + MTX")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0)


#Fig. 5A

#CPOP delta heatmaps for insertions
cpop_data_delta_ins <- cpop_data %>%
  filter(mut_type == "ins") %>%
  mutate(dtemp = score_30 - score_37,
         dmtx_30 = score_30_MTX - score_30,
         dmtx_37 = score_37_MTX - score_37) %>%
  select(pos, dtemp, dmtx_30, dmtx_37)
cpop_data_delta_ins_t <- t(cpop_data_delta_ins)
cpop_data_delta_ins_heat <- cpop_data_delta_ins_t[-1, , drop = FALSE]
colnames(cpop_data_delta_ins_heat) <- cpop_data_delta_ins_t[1,]
draw(Heatmap(cpop_data_delta_ins_heat,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(0, 1.113), c("white", "magenta")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*"CPOP score"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(0, 1),
                                         labels = c(0, 1)),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                          legend_gp = gpar(fill = c("dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 12),
                                          row_gap = unit(2, "mm"))))


#CPOP delta heatmaps for deletions
cpop_data_delta_del <- cpop_data %>%
  filter(mut_type == "del") %>%
  mutate(dtemp = score_30 - score_37,
         dmtx_30 = score_30_MTX - score_30,
         dmtx_37 = score_37_MTX - score_37) %>%
  select(pos, dtemp, dmtx_30, dmtx_37)
cpop_data_delta_del_t <- t(cpop_data_delta_del)
cpop_data_delta_del_heat <- cpop_data_delta_del_t[-1, , drop = FALSE]
colnames(cpop_data_delta_del_heat) <- cpop_data_delta_del_t[1,]
draw(Heatmap(cpop_data_delta_del_heat,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(0, 1.054), c("white", "magenta")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*"CPOP score"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(0, 1),
                                         labels = c(0, 1)),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom", 
     annotation_legend_side = "bottom",
     annotation_legend_list = list(Legend(labels = c("Missing variant"),
                                          legend_gp = gpar(fill = c("dimgrey")),
                                          border = c("black", "black"),
                                          labels_gp = gpar(fontsize = 12),
                                          row_gap = unit(2, "mm"))))


#Fig. 6A

#Scatter plots for dplddt and ddg for insertions
cpop_ins_gather <- cpop_data %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_ins_dplddt_ddg <- cpop_data %>%
  filter(mut_type == "ins") %>%
  merge(ins_dplddt_ddg, by = "pos") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cor(cpop_ins_dplddt_ddg$dplddt, cpop_ins_dplddt_ddg$ddg, method = "spearman")
ggplot() +
  geom_point(cpop_ins_dplddt_ddg, mapping = aes(dplddt, ddg), color = "black", size = 3) +
  geom_point(cpop_ins_dplddt_ddg, mapping = aes(dplddt, ddg, color = score), size = 2) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "longdash") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "longdash") +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  scale_colour_gradient2(limits = c(-0.228, 2.228), 
                         breaks = c(0, 1, 2), 
                         midpoint = 1, 
                         low = 'darkred', 
                         mid = 'white', 
                         high = 'darkblue', 
                         labels = c(0, 1, 2)) +
  labs(x = expression(Delta*pLDDT), y = expression(Delta*Delta*G), color = "CPOP score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14), panel.spacing = unit(1, "lines"), text = element_text(size = 14))

#Scatter plots for dplddt and ddg for deletions
cpop_del_gather <- cpop_data %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_del_dplddt_ddg <- cpop_data %>%
  filter(mut_type == "del") %>%
  merge(del_dplddt_ddg, by = "pos") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cor(cpop_del_dplddt_ddg$dplddt, cpop_del_dplddt_ddg$ddg, method = "spearman")
ggplot() +
  geom_point(cpop_del_dplddt_ddg, mapping = aes(dplddt, ddg), color = "black", size = 3) +
  geom_point(cpop_del_dplddt_ddg, mapping = aes(dplddt, ddg, color = score), size = 2) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5, linetype = "longdash") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5, linetype = "longdash") +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  scale_colour_gradient2(limits = c(-0.288, 2.288), 
                         breaks = c(0, 1, 2), 
                         midpoint = 1, 
                         low = 'darkred', 
                         mid = 'white', high = 'darkblue', 
                         labels = c(0, 1, 2)) +
  labs(x = expression(Delta*pLDDT), y = expression(Delta*Delta*G), color = "CPOP score") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14), panel.spacing = unit(1, "lines"), text = element_text(size = 14))


#Fig. 6B

#ROC curve for insertions
ins_thres <- cpop_data_ROC_ins %>%
  filter(mut_type == "ins") %>%
  merge(ins_dplddt_ddg, by = "pos") %>%
  mutate(cutoff_30 = ifelse(score_30 > 0.5, 1, 0),
         cutoff_30_MTX = ifelse(score_30_MTX > 0.5, 1, 0),
         cutoff_37 = ifelse(score_37 > 0.5, 1, 0),
         cutoff_37_MTX = ifelse(score_37_MTX > 0.5, 1, 0))
ggroc(list(dplddt_30 = roc(ins_thres$cutoff_30, ins_thres$dplddt),
           ddG_30 = roc(ins_thres$cutoff_30, ins_thres$ddg),
           dplddt_30_MTX = roc(ins_thres$cutoff_30_MTX, ins_thres$dplddt),
           ddG_30_MTX = roc(ins_thres$cutoff_30_MTX, ins_thres$ddg),
           dplddt_37 = roc(ins_thres$cutoff_37, ins_thres$dplddt),
           ddG_37 = roc(ins_thres$cutoff_37, ins_thres$ddg),
           dplddt_37_MTX = roc(ins_thres$cutoff_37_MTX, ins_thres$dplddt),
           ddG_37_MTX = roc(ins_thres$cutoff_37_MTX, ins_thres$ddg)),
      legacy.axes = TRUE, linewidth = 2, aes = c("color", "linetype")) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", labels = c(bquote(Delta*pLDDT: "30 °C " (AUC == .(round(auc(ins_thres$cutoff_30, ins_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "30 °C " (AUC == .(round(auc(ins_thres$cutoff_30, ins_thres$ddg), 2)))),
                                       bquote(Delta*pLDDT: "30 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_30_MTX, ins_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "30 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_30_MTX, ins_thres$ddg), 2)))),
                                       bquote(Delta*pLDDT: "37 °C " (AUC == .(round(auc(ins_thres$cutoff_37, ins_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "37 °C " (AUC == .(round(auc(ins_thres$cutoff_37, ins_thres$ddg), 2)))),
                                       bquote(Delta*pLDDT: "37 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_37_MTX, ins_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "37 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_37_MTX, ins_thres$ddg), 2))))),
                            values = c(1, 2, 1, 2, 1, 2, 1, 2)) +
  scale_color_manual("", labels = c(bquote(Delta*pLDDT: "30 °C " (AUC == .(round(auc(ins_thres$cutoff_30, ins_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "30 °C " (AUC == .(round(auc(ins_thres$cutoff_30, ins_thres$ddg), 2)))),
                                    bquote(Delta*pLDDT: "30 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_30_MTX, ins_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "30 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_30_MTX, ins_thres$ddg), 2)))),
                                    bquote(Delta*pLDDT: "37 °C " (AUC == .(round(auc(ins_thres$cutoff_37, ins_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "37 °C " (AUC == .(round(auc(ins_thres$cutoff_37, ins_thres$ddg), 2)))),
                                    bquote(Delta*pLDDT: "37 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_37_MTX, ins_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "37 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_37_MTX, ins_thres$ddg), 2))))),
                         values = c("orange", "orange", "limegreen", "limegreen", "magenta3", "magenta3", "dodgerblue2", "dodgerblue2")) +
  geom_abline(linetype = "dashed", linewidth = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0, aspect.ratio = 1)

#ROC curve for deletions
del_thres <- cpop_data_ROC_del %>%
  filter(mut_type == "del") %>%
  merge(del_dplddt_ddg, by = "pos") %>%
  mutate(cutoff_30 = ifelse(score_30 > 0.5, 1, 0),
         cutoff_30_MTX = ifelse(score_30_MTX > 0.5, 1, 0),
         cutoff_37 = ifelse(score_37 > 0.5, 1, 0),
         cutoff_37_MTX = ifelse(score_37_MTX > 0.5, 1, 0))
ggroc(list(dplddt_30 = roc(del_thres$cutoff_30, del_thres$dplddt),
           ddG_30 = roc(del_thres$cutoff_30, del_thres$ddg),
           dplddt_30_MTX = roc(del_thres$cutoff_30_MTX, del_thres$dplddt),
           ddG_30_MTX = roc(del_thres$cutoff_30_MTX, del_thres$ddg),
           dplddt_37 = roc(del_thres$cutoff_37, del_thres$dplddt),
           ddG_37 = roc(del_thres$cutoff_37, del_thres$ddg),
           dplddt_37_MTX = roc(del_thres$cutoff_37_MTX, del_thres$dplddt),
           ddG_37_MTX = roc(del_thres$cutoff_37_MTX, del_thres$ddg)),
      legacy.axes = TRUE, linewidth = 2, aes = c("color", "linetype")) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", labels = c(bquote(Delta*pLDDT: "30 °C " (AUC == .(round(auc(del_thres$cutoff_30, del_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "30 °C " (AUC == .(round(auc(del_thres$cutoff_30, del_thres$ddg), 2)))),
                                       bquote(Delta*pLDDT: "30 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_30_MTX, del_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "30 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_30_MTX, del_thres$ddg), 2)))),
                                       bquote(Delta*pLDDT: "37 °C " (AUC == .(round(auc(del_thres$cutoff_37, del_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "37 °C " (AUC == .(round(auc(del_thres$cutoff_37, del_thres$ddg), 2)))),
                                       bquote(Delta*pLDDT: "37 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_37_MTX, del_thres$dplddt), 2)))),
                                       bquote(Delta*Delta * G: "37 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_37_MTX, del_thres$ddg), 2))))),
                        values = c(1, 2, 1, 2, 1, 2, 1, 2)) +
  scale_color_manual("", labels = c(bquote(Delta*pLDDT: "30 °C " (AUC == .(round(auc(del_thres$cutoff_30, del_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "30 °C " (AUC == .(round(auc(del_thres$cutoff_30, del_thres$ddg), 2)))),
                                    bquote(Delta*pLDDT: "30 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_30_MTX, del_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "30 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_30_MTX, del_thres$ddg), 2)))),
                                    bquote(Delta*pLDDT: "37 °C " (AUC == .(round(auc(del_thres$cutoff_37, del_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "37 °C " (AUC == .(round(auc(del_thres$cutoff_37, del_thres$ddg), 2)))),
                                    bquote(Delta*pLDDT: "37 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_37_MTX, del_thres$dplddt), 2)))),
                                    bquote(Delta*Delta * G: "37 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_37_MTX, del_thres$ddg), 2))))),
                     values = c("orange", "orange", "limegreen", "limegreen", "magenta3", "magenta3", "dodgerblue2", "dodgerblue2")) +
  geom_abline(linetype = "dashed", linewidth = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0, aspect.ratio = 1)
