#Load in required packages
library(tidyverse)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)

#Change the path to your directory
path <- "your_path_goes_here"

#Read in dataframes
cpop_data <- read.csv(file.path(path, "cpop_data.csv"))
dssp <- read.csv(file.path(path, "1u72_dssp.csv"), sep = ",")
ins_dplddt_ddg <- read.csv(file.path(path, "ins_dplddt_ddg.csv"))
del_dplddt_ddg <- read.csv(file.path(path, "del_dplddt_ddg.csv"))
mtx_dist <- read.csv(file.path(path, "mtx_dist.csv"))

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


#Fig. S2 (CPOP delta insertion/deletion heatmaps)

#Calculate CPOP delta insertion/deletion scores
cpop_data_delta_indel <- tibble(pos = 1:187) %>%
  tibble(dindel_30 = cpop_data$score_30[cpop_data$mut_type == "ins"] - cpop_data$score_30[cpop_data$mut_type == "del"],
         dindel_30_MTX = cpop_data$score_30_MTX[cpop_data$mut_type == "ins"] - cpop_data$score_30_MTX[cpop_data$mut_type == "del"],
         dindel_37 = cpop_data$score_37[cpop_data$mut_type == "ins"] - cpop_data$score_37[cpop_data$mut_type == "del"],
         dindel_37_MTX = cpop_data$score_37_MTX[cpop_data$mut_type == "ins"] - cpop_data$score_37_MTX[cpop_data$mut_type == "del"])

#CPOP delta insertion/deletion heatmaps 
cpop_data_delta_indel_t <- t(cpop_data_delta_indel)
cpop_data_delta_indel_heat <- cpop_data_delta_indel_t[-1, , drop = FALSE]
colnames(cpop_data_delta_indel_heat) <- cpop_data_delta_indel_t[1,]
draw(Heatmap(cpop_data_delta_indel_heat,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(-0.9513439, 0, 0.9513439), c("cyan", "white", "magenta")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*"indel score"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-0.9513439, 0, 0.9513439),
                                         labels = c("Insertion\nmore unstable", "Equally stable", "Deletion\nmore unstable")),
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


#Fig. S3 (rSASA)

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

#Scatter plots for relative solvent accessible surface area (rSASA) for insertions
dssp_cpop_ins <- merge(cpop_data, dssp_group, by = "pos") %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
dssp_cpop_ins$sec <- factor(dssp_cpop_ins$sec , levels = c("loops", "helices", "sheets"))
ggplot() +
  geom_point(dssp_cpop_ins, mapping = aes(score, Relative.ASA), color = "black", size = 3) +
  geom_point(dssp_cpop_ins, mapping = aes(score, Relative.ASA, color = sec), size = 2) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  labs(x = "CPOP score", y = "rSASA") +
  scale_color_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                    labels=c("Loops", "Helices", "Sheets"),
                    name = "Secondary structure") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))

#Scatter plots for relative solvent accessible surface area (rSASA) for deletions
dssp_cpop_del <- merge(cpop_data, dssp_group, by = "pos") %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
dssp_cpop_del$sec <- factor(dssp_cpop_del$sec , levels = c("loops", "helices", "sheets"))
ggplot() +
  geom_point(dssp_cpop_del, mapping = aes(score, Relative.ASA), color = "black", size = 3) +
  geom_point(dssp_cpop_del, mapping = aes(score, Relative.ASA, color = sec), size = 2) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  labs(x = "CPOP score", y = "rSASA") +
  scale_color_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                     labels=c("Loops", "Helices", "Sheets"),
                     name = "Secondary structure") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))


#Fig. S4 (MTX distance)

#Scatter plots for distance to MTX for dMTX scores for insertions
cpop_data_delta_ins <- cpop_data %>%
  filter(mut_type == "ins") %>%
  mutate(dmtx_30 = score_30_MTX - score_30,
         dmtx_37 = score_37_MTX - score_37) %>%
  select(pos, dmtx_30, dmtx_37) %>%
  merge(mtx_dist, by = "pos") %>%
  gather(experiment, score, c(dmtx_30, dmtx_37))
ggplot() +
  geom_point(cpop_data_delta_ins, mapping = aes(score, MTXdist), size = 2) +
  facet_grid(~ experiment, labeller = as_labeller(c(dmtx_30 = "30 °C", dmtx_37 = "37 °C"))) +
  labs(x = expression(Delta*MTX~score), y = "Distance to MTX (Å)") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))

#Scatter plots for distance to MTX for dMTX scores for deletions
cpop_data_delta_del <- cpop_data %>%
  filter(mut_type == "del") %>%
  mutate(dmtx_30 = score_30_MTX - score_30,
         dmtx_37 = score_37_MTX - score_37) %>%
  select(pos, dmtx_30, dmtx_37) %>%
  merge(mtx_dist, by = "pos") %>%
  gather(experiment, score, c(dmtx_30, dmtx_37))
ggplot() +
  geom_point(cpop_data_delta_del, mapping = aes(score, MTXdist), size = 2) +
  facet_grid(~ experiment, labeller = as_labeller(c(dmtx_30 = "30 °C", dmtx_37 = "37 °C"))) +
  labs(x = expression(Delta*MTX~score), y = "Distance to MTX (Å)") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))


#Fig. S5 (pLDDT)

#Average pLDDT for insertion and deletion variants
ggplot() +
  geom_histogram(ins_dplddt_ddg, mapping = aes(plddt, fill = "Insertion", color = "Insertion"), binwidth = 0.5, alpha = 0.5) +
  geom_histogram(del_dplddt_ddg, mapping = aes(plddt, fill = "Deletion", color = "Deletion"), binwidth = 0.5, alpha = 0.5) +
  scale_x_continuous(breaks = c(86, 88, 90, 92, 94, 96)) +
  scale_colour_manual(name = NULL, values = c("Insertion" = "black", "Deletion"="black"), breaks = c("Insertion","Deletion")) +
  scale_fill_manual(name = NULL, values = adjustcolor(c("lightblue","salmon"), alpha.f = 0.6), breaks = c("Insertion","Deletion")) +
  guides(fill = guide_legend(override.aes = list(alpha = c(0.4, 0.4)))) +
  xlab("AlphaFold pLDDT") +
  ylab("Count") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))
