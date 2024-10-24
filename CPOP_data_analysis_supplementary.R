#Load in required packages
library(tidyverse)
library(forstringr)
library(ggpubr)
library(ggpmisc)
library(ComplexHeatmap) #install from Bioconductor
library(circlize)
library(patchwork)

#Change the path to your directory
path <- "your_path"

#Read in dataframes
cpop_data <- read.csv(file.path(path, "cpop_data.csv"))
dssp <- read.csv(file.path(path, "rSASA.csv"), sep = ",")
ins_dplddt_ddg <- read.csv(file.path(path, "ins_dplddt_ddg.csv"))
del_dplddt_ddg <- read.csv(file.path(path, "del_dplddt_ddg.csv"))
mtx_dist <- read.csv(file.path(path, "mtx_dist.csv"))
cpop_data_pre_rescale <- read.csv(file.path(path, "cpop_data_pre_rescale.csv"))
tile1_count <- read.csv(file.path(path, "tile1.csv"))
tile2_count <- read.csv(file.path(path, "tile2.csv"))
tile3_count <- read.csv(file.path(path, "tile3.csv"))
tile4_count <- read.csv(file.path(path, "tile4.csv"))
tile5_count <- read.csv(file.path(path, "tile5.csv"))
del_sequence_alignment <- read.csv(file.path(path, "del_sequence_alignment.csv"))
cpop_data_ROC_ins <- read.csv(file.path(path, "cpop_data_ROC_ins.csv"))
cpop_data_ROC_del <- read.csv(file.path(path, "cpop_data_ROC_del.csv")) 
ins_esm1b <- read.csv(file.path(path, "ins_esm1b.csv"))
del_esm1b <- read.csv(file.path(path, "del_esm1b.csv"))


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

#Create dataframes for each type of mutation prior to rescaling
cpop_pre_ins_gather <- cpop_data_pre_rescale %>%
  filter(mut_type == "ins") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_pre_ins_syn_gather <- cpop_data_pre_rescale %>%
  filter(mut_type == "syn") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_pre_ins_stop_gather <- cpop_data_pre_rescale %>%
  filter(mut_type == "stop") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_pre_del_gather <- cpop_data_pre_rescale %>%
  filter(mut_type == "del") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_pre_del_syn_gather <- cpop_data_pre_rescale %>%
  filter(mut_type == "syn") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
cpop_pre_del_stop_gather <- cpop_data_pre_rescale %>%
  filter(mut_type == "stop") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))



#Fig. S2 (Score distributions prior to rescaling)

#Score distributions of raw CPOP scores for insertions
ggplot() + 
  geom_histogram(cpop_pre_ins_syn_gather, mapping = aes(score, fill = "syn", color = "syn"), binwidth = 0.1, alpha = 1) +
  geom_histogram(cpop_pre_ins_stop_gather, mapping = aes(score, fill = "stop", color = "stop"), binwidth = 0.1, alpha = 1) +
  geom_histogram(cpop_pre_ins_gather, mapping = aes(score, fill = "ins", color = "ins"), binwidth = 0.1, alpha = 0.8) +
  geom_segment(cpop_ins_gather, mapping = aes(x = 0, y = 0, xend = 0, yend = 113), color = "black", linetype = "longdash") +
  annotate("label", x = 0, y = 115, label = "WT" , color = "black", alpha = 0.3) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  scale_colour_manual(name = NULL, values = c("ins" = "black", "syn"="black", "stop" = "black"), labels = c("ins" = "Insertion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("ins","syn", "stop")) +
  scale_fill_manual(name = NULL, values = adjustcolor(c("darkorange","white", "darkred"), alpha.f = 0.3), labels = c("ins" = "Insertion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("ins","syn", "stop")) +
  guides(fill = guide_legend(override.aes = list(alpha = c(0.8, 1, 1)))) +
  xlab("CPOP score") +
  ylab("Count") +
  ylim(0,118) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14), panel.spacing = unit(1, "lines"), text = element_text(size = 14))

#Score distributions of raw CPOP scores for deletions
ggplot() + 
  geom_histogram(cpop_pre_del_syn_gather, mapping = aes(score, fill = "syn", color = "syn"), binwidth = 0.1, alpha = 1) +
  geom_histogram(cpop_pre_del_stop_gather, mapping = aes(score, fill = "stop", color = "stop"), binwidth = 0.1, alpha = 1) +
  geom_histogram(cpop_pre_del_gather, mapping = aes(score, fill = "del", color = "del"), binwidth = 0.1, alpha = 0.8) +
  geom_segment(cpop_del_gather, mapping = aes(x = 0, y = 0, xend = 0, yend = 113), color = "black", linetype = "longdash") +
  annotate("label", x = 0, y = 115, label = "WT" , color = "black", alpha = 0.3) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  scale_colour_manual(name = NULL, values = c("del" = "black", "syn"="black", "stop" = "black"), labels = c("del" = "Deletion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("del","syn", "stop")) +
  scale_fill_manual(name = NULL, values = adjustcolor(c("darkorange","white", "darkred"), alpha.f = 0.3), labels = c("del" = "Deletion", "syn" = "Synonymous", "stop" = "Nonsense"), breaks = c("del","syn", "stop")) +
  guides(fill = guide_legend(override.aes = list(alpha = c(0.8, 1, 1)))) +
  xlab("CPOP score") +
  ylab("Count") +
  ylim(0,118) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 14), legend.text = element_text(size = 14), panel.spacing = unit(1, "lines"), text = element_text(size = 14))


#Fig. S3 (Count correlation)

#Sequencing count correlation between triplicates across tiles and conditions after selection

#Tile 1
tile1_rep1 <- tile1_count%>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no.sele, count_rep1_30, count_rep1_30.MTX, count_rep1_37, count_rep1_37.MTX))
tile1_rep2 <- tile1_count%>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no.sele, count_rep2_30, count_rep2_30.MTX, count_rep2_37, count_rep2_37.MTX))
tile1_rep3 <- tile1_count%>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no.sele, count_rep3_30, count_rep3_30.MTX, count_rep3_37, count_rep3_37.MTX))
tile1_merged <- cbind(tile1_rep1, tile1_rep2, tile1_rep3)
colnames(tile1_merged) <- make.unique(names(tile1_merged)) 
tile1_merged_plot <- tile1_merged %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no.sele") &
                              str_detect(experiment_rep2, "no.sele") &
                              str_detect(experiment_rep3, "no.sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_detect(experiment_rep1, "_30.MTX") & 
                                            str_detect(experiment_rep2, "_30.MTX") &
                                            str_detect(experiment_rep3, "_30.MTX"), "30 °C + MTX",
                                          ifelse(str_right(experiment_rep1, 2) == "37" & 
                                                   str_right(experiment_rep2, 2) == "37" &
                                                   str_right(experiment_rep3, 2) == "37", "37 °C",
                                                 ifelse(str_detect(experiment_rep1, "_37.MTX") &
                                                          str_detect(experiment_rep2, "_37.MTX") &
                                                          str_detect(experiment_rep3, "_37.MTX"), "37 °C + MTX", NA))))))
ggplot(tile1_merged_plot, aes(count_rep1, count_rep2, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) + 
  labs(x = "Replicate 1", y = "Replicate 2") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile1_merged_plot, aes(count_rep1, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 1", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile1_merged_plot, aes(count_rep2, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 2", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)

#Tile 2
tile2_rep1 <- tile2_count %>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no.sele, count_rep1_30, count_rep1_30.MTX, count_rep1_37, count_rep1_37.MTX))
tile2_rep2 <- tile2_count %>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no.sele, count_rep2_30, count_rep2_30.MTX, count_rep2_37, count_rep2_37.MTX))
tile2_rep3 <- tile2_count %>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no.sele, count_rep3_30, count_rep3_30.MTX, count_rep3_37, count_rep3_37.MTX))
tile2_merged <- cbind(tile2_rep1, tile2_rep2, tile2_rep3)
colnames(tile2_merged) <- make.unique(names(tile2_merged)) 
tile2_merged_plot <- tile2_merged %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no.sele") &
                              str_detect(experiment_rep2, "no.sele") &
                              str_detect(experiment_rep3, "no.sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_detect(experiment_rep1, "_30.MTX") & 
                                            str_detect(experiment_rep2, "_30.MTX") &
                                            str_detect(experiment_rep3, "_30.MTX"), "30 °C + MTX",
                                          ifelse(str_right(experiment_rep1, 2) == "37" & 
                                                   str_right(experiment_rep2, 2) == "37" &
                                                   str_right(experiment_rep3, 2) == "37", "37 °C",
                                                 ifelse(str_detect(experiment_rep1, "_37.MTX") &
                                                          str_detect(experiment_rep2, "_37.MTX") &
                                                          str_detect(experiment_rep3, "_37.MTX"), "37 °C + MTX", NA))))))
ggplot(tile2_merged_plot, aes(count_rep1, count_rep2, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) + 
  labs(x = "Replicate 1", y = "Replicate 2") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile2_merged_plot, aes(count_rep1, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 1", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile2_merged_plot, aes(count_rep2, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 2", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)

#Tile 3
tile3_rep1 <- tile3_count %>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no.sele, count_rep1_30, count_rep1_30.MTX, count_rep1_37, count_rep1_37.MTX))
tile3_rep2 <- tile3_count %>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no.sele, count_rep2_30, count_rep2_30.MTX, count_rep2_37, count_rep2_37.MTX))
tile3_rep3 <- tile3_count %>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no.sele, count_rep3_30, count_rep3_30.MTX, count_rep3_37, count_rep3_37.MTX))
tile3_merged <- cbind(tile3_rep1, tile3_rep2, tile3_rep3)
colnames(tile3_merged) <- make.unique(names(tile3_merged)) 
tile3_merged_plot <- tile3_merged %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no.sele") &
                              str_detect(experiment_rep2, "no.sele") &
                              str_detect(experiment_rep3, "no.sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_detect(experiment_rep1, "_30.MTX") & 
                                            str_detect(experiment_rep2, "_30.MTX") &
                                            str_detect(experiment_rep3, "_30.MTX"), "30 °C + MTX",
                                          ifelse(str_right(experiment_rep1, 2) == "37" & 
                                                   str_right(experiment_rep2, 2) == "37" &
                                                   str_right(experiment_rep3, 2) == "37", "37 °C",
                                                 ifelse(str_detect(experiment_rep1, "_37.MTX") &
                                                          str_detect(experiment_rep2, "_37.MTX") &
                                                          str_detect(experiment_rep3, "_37.MTX"), "37 °C + MTX", NA))))))
ggplot(tile3_merged_plot, aes(count_rep1, count_rep2, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) + 
  labs(x = "Replicate 1", y = "Replicate 2") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile3_merged_plot, aes(count_rep1, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 1", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile3_merged_plot, aes(count_rep2, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 2", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)

#Tile 4
tile4_rep1 <- tile4_count %>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no.sele, count_rep1_30, count_rep1_30.MTX, count_rep1_37, count_rep1_37.MTX))
tile4_rep2 <- tile4_count %>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no.sele, count_rep2_30, count_rep2_30.MTX, count_rep2_37, count_rep2_37.MTX))
tile4_rep3 <- tile4_count %>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no.sele, count_rep3_30, count_rep3_30.MTX, count_rep3_37, count_rep3_37.MTX))
tile4_merged <- cbind(tile4_rep1, tile4_rep2, tile4_rep3)
colnames(tile4_merged) <- make.unique(names(tile4_merged)) 
tile4_merged_plot <- tile4_merged %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no.sele") &
                              str_detect(experiment_rep2, "no.sele") &
                              str_detect(experiment_rep3, "no.sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_detect(experiment_rep1, "_30.MTX") & 
                                            str_detect(experiment_rep2, "_30.MTX") &
                                            str_detect(experiment_rep3, "_30.MTX"), "30 °C + MTX",
                                          ifelse(str_right(experiment_rep1, 2) == "37" & 
                                                   str_right(experiment_rep2, 2) == "37" &
                                                   str_right(experiment_rep3, 2) == "37", "37 °C",
                                                 ifelse(str_detect(experiment_rep1, "_37.MTX") &
                                                          str_detect(experiment_rep2, "_37.MTX") &
                                                          str_detect(experiment_rep3, "_37.MTX"), "37 °C + MTX", NA))))))
ggplot(tile4_merged_plot, aes(count_rep1, count_rep2, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) + 
  labs(x = "Replicate 1", y = "Replicate 2") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile4_merged_plot, aes(count_rep1, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 1", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile4_merged_plot, aes(count_rep2, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 2", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)

#Tile 5
tile5_rep1 <- tile5_count %>%
  gather(experiment_rep1, count_rep1, c(count_rep1_no.sele, count_rep1_30, count_rep1_30.MTX, count_rep1_37, count_rep1_37.MTX))
tile5_rep2 <- tile5_count %>%
  gather(experiment_rep2, count_rep2, c(count_rep2_no.sele, count_rep2_30, count_rep2_30.MTX, count_rep2_37, count_rep2_37.MTX))
tile5_rep3 <- tile5_count %>%
  gather(experiment_rep3, count_rep3, c(count_rep3_no.sele, count_rep3_30, count_rep3_30.MTX, count_rep3_37, count_rep3_37.MTX))
tile5_merged <- cbind(tile5_rep1, tile5_rep2, tile5_rep3)
colnames(tile5_merged) <- make.unique(names(tile5_merged)) 
tile5_merged_plot <- tile5_merged %>%
  mutate(Condition = ifelse(str_detect(experiment_rep1, "no.sele") &
                              str_detect(experiment_rep2, "no.sele") &
                              str_detect(experiment_rep3, "no.sele"), "no selection",
                            ifelse(str_right(experiment_rep1, 2) == "30" &
                                     str_right(experiment_rep2, 2) == "30" &
                                     str_right(experiment_rep3, 2) == "30", "30 °C", 
                                   ifelse(str_detect(experiment_rep1, "_30.MTX") & 
                                            str_detect(experiment_rep2, "_30.MTX") &
                                            str_detect(experiment_rep3, "_30.MTX"), "30 °C + MTX",
                                          ifelse(str_right(experiment_rep1, 2) == "37" & 
                                                   str_right(experiment_rep2, 2) == "37" &
                                                   str_right(experiment_rep3, 2) == "37", "37 °C",
                                                 ifelse(str_detect(experiment_rep1, "_37.MTX") &
                                                          str_detect(experiment_rep2, "_37.MTX") &
                                                          str_detect(experiment_rep3, "_37.MTX"), "37 °C + MTX", NA))))))
ggplot(tile5_merged_plot, aes(count_rep1, count_rep2, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) + 
  labs(x = "Replicate 1", y = "Replicate 2") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile5_merged_plot, aes(count_rep1, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 1", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)
ggplot(tile5_merged_plot, aes(count_rep2, count_rep3, color = Condition)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  stat_cor(method = "pearson", aes(label = after_stat(r.label)), show.legend = F) +
  stat_poly_line(show.legend = F) +
  labs(x = "Replicate 2", y = "Replicate 3") + 
  theme_bw() +
  theme(aspect.ratio = 1)


#Fig. S4 (SD heatmaps)

#Standard deviations for indel variants

#SD heatmaps for insertions
cpop_ins_SD <- cpop_data %>%
  filter(mut_type == "ins") %>%
  select(pos, SD_30, SD_30_MTX, SD_37, SD_37_MTX)
cpop_ins_SD_compl <- t(cpop_ins_SD)
cpop_ins_SD_compl_heatmap <- cpop_ins_SD_compl[-1,]
colnames(cpop_ins_SD_compl_heatmap) <- cpop_ins_SD_compl[1,]
draw(Heatmap(cpop_ins_SD_compl_heatmap,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(0.001, 0.547), c("white", "purple4")),
             name = ("SD"),
             heatmap_legend_param = list(direction = "horizontal",
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(0.001, 0.547),
                                         labels = c("0.001", "0.547")),
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

#SD heatmaps for deletions
cpop_del_SD <- cpop_data %>%
  filter(mut_type == "del") %>%
  select(pos, SD_30, SD_30_MTX, SD_37, SD_37_MTX)
cpop_del_SD_compl <- t(cpop_del_SD)
cpop_del_SD_compl_heatmap <- cpop_del_SD_compl[-1,]
colnames(cpop_del_SD_compl_heatmap) <- cpop_del_SD_compl[1,]
draw(Heatmap(cpop_del_SD_compl_heatmap,
                     cluster_rows = FALSE, 
                     cluster_columns = FALSE,
                     colorRamp2(c(0.002, 0.565), c("white", "purple4")),
                     name = ("SD"),
                     heatmap_legend_param = list(direction = "horizontal",
                                                 title_position = "topcenter", 
                                                 legend_width = unit(6, "cm"),
                                                 labels_gp = gpar(fontsize = 12),
                                                 title_gp = gpar(fontsize = 14, fontface = "bold"),
                                                 border = "black", 
                                                 at = c(0.002, 0.565),
                                                 labels = c("0.002", "0.565")),
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


#Fig. S6 (CPOP delta insertion/deletion heatmaps)

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
             colorRamp2(c(-0.9335, 0, 0.9335), c("cyan", "white", "magenta")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*"indel score"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-0.9335, 0, 0.9335),
                                         labels = c("-1", "0", "1")),
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


#Fig. S7 (rSASA)

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
ggplot(dssp_cpop_ins, mapping = aes(score, Relative.ASA), color = "black") +
  geom_point(size = 3) +
  geom_point(dssp_cpop_ins, mapping = aes(score, Relative.ASA, color = sec), size = 2) +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), label.y = 1.05) +
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
ggplot(dssp_cpop_del, mapping = aes(score, Relative.ASA), color = "black") +
  geom_point(size = 3) +
  geom_point(dssp_cpop_del, mapping = aes(score, Relative.ASA, color = sec), size = 2) +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), label.y = 1.05) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  labs(x = "CPOP score", y = "rSASA") +
  scale_color_manual(values = c("#73A222", "#FFF670", "#CF47FF"),
                     labels=c("Loops", "Helices", "Sheets"),
                     name = "Secondary structure") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))


#Fig. S8 (MTX distance)

#Scatter plots for distance to MTX for dMTX scores for insertions
cpop_data_delta_ins <- cpop_data %>%
  filter(mut_type == "ins") %>%
  mutate(dmtx_30 = score_30_MTX - score_30,
         dmtx_37 = score_37_MTX - score_37) %>%
  select(pos, dmtx_30, dmtx_37) %>%
  merge(mtx_dist, by = "pos") %>%
  gather(experiment, score, c(dmtx_30, dmtx_37))
ggplot(cpop_data_delta_ins, mapping = aes(score, MTXdist), size = 2) +
  geom_point() +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), label.y = 25.5) + 
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
ggplot(cpop_data_delta_del, mapping = aes(score, MTXdist), size = 2) +
  geom_point() +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), label.y = 25.5) + 
  facet_grid(~ experiment, labeller = as_labeller(c(dmtx_30 = "30 °C", dmtx_37 = "37 °C"))) +
  labs(x = expression(Delta*MTX~score), y = "Distance to MTX (Å)") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))


#Fig. S9 (dplddt and ddg vs. CPOP scores)

#Scatter plots for dplddt and ddg vs. CPOP scores for insertions
cpop_ins_dplddt_ddg <- cpop_data %>%
  filter(mut_type == "ins") %>%
  merge(ins_dplddt_ddg, by = "pos") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
ggplot(cpop_ins_dplddt_ddg, mapping = aes(score, ddg)) +
         geom_point(size = 1) +
         stat_cor(method = "spearman", aes(label = after_stat(r.label)), label.y = 37) +
         facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
         labs(x = expression("CPOP score"), y = expression(Delta*Delta*G)) +
         scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))

ggplot(cpop_ins_dplddt_ddg, mapping = aes(score, dplddt)) +
         geom_point(size = 1) +
         stat_cor(method = "spearman", aes(label = after_stat(r.label)), label.y = 1) +
         facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
         labs(x = expression("CPOP score"), y = expression(Delta*pLDDT)) +
         scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))


#Fig. S10 (dplddt and ddg vs. CPOP scores)

#Scatter plots for dplddt and ddg vs. CPOP scores for deletions
cpop_del_dplddt_ddg <- cpop_data %>%
  filter(mut_type == "del") %>%
  merge(del_dplddt_ddg, by = "pos") %>%
  gather(experiment, score, c(score_30, score_30_MTX, score_37, score_37_MTX))
ggplot(cpop_del_dplddt_ddg, mapping = aes(score, ddg)) +
  geom_point(size = 1) +
  stat_cor(method = "spearman", aes(label = after_stat(r.label)), , label.y = 22) +
  facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
  labs(x = expression("CPOP score"), y = expression(Delta*Delta*G)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))

ggplot(cpop_del_dplddt_ddg, mapping = aes(score, dplddt)) +
         geom_point(size = 1) +
         stat_cor(method = "spearman", aes(label = after_stat(r.label)), label.y = 1) +
         facet_grid(~ experiment, labeller = as_labeller(c(score_30 = "30 °C", score_30_MTX = "30 °C + MTX", score_37 = "37 °C", score_37_MTX = "37 °C + MTX"))) +
         labs(x = expression("CPOP score"), y = expression(Delta*pLDDT)) +
         scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
         theme_bw() +
         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14))


#Fig. S11 (pLDDT)

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


#Fig. S12 (dpLDDT and ddG heatmaps)

#dpLDDT and ddG heatmaps for insertions
ins_dplddt <- ins_dplddt_ddg %>%
  select(pos, dplddt)
ins_dplddt_compl <- t(ins_dplddt)
ins_dplddt_heatmap <- matrix(ins_dplddt_compl[-1,], ncol = 187)
colnames(ins_dplddt_heatmap) <- ins_dplddt_compl[1,]
draw(Heatmap(ins_dplddt_heatmap,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(-9.451868782, 0), c("purple", "white")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*"pLDDT"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-10, -7.5, -5, -2.5, 0),
                                         labels = c("-10", "-7.5", "-5", "-2.5", "0")),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom")

ins_ddg <- ins_dplddt_ddg %>%
  select(pos, ddg)
ins_ddg_compl <- t(ins_ddg)
ins_ddg_heatmap <- matrix(ins_ddg_compl[-1,], ncol = 187)
colnames(ins_ddg_heatmap) <- ins_ddg_compl[1,]
draw(Heatmap(ins_ddg_heatmap,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(-10, 0, 10), c("darkgreen", "white", "darkorange")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*Delta*"G"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-10, -5, 0, 5, 10),
                                         labels = c("-10", "-5", "0", "5", "10")),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom")

#dpLDDT and ddG heatmaps for deletions
del_dplddt <- del_dplddt_ddg %>%
  select(pos, dplddt)
del_dplddt_compl <- t(del_dplddt)
del_dplddt_heatmap <- matrix(del_dplddt_compl[-1,], ncol = 187)
colnames(del_dplddt_heatmap) <- del_dplddt_compl[1,]
draw(Heatmap(del_dplddt_heatmap,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(-2.822004747, 0), c("purple", "white")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*"pLDDT"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-3, -2, -1, 0),
                                         labels = c("-3", "-2", "-1", "0")),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom")

del_ddg <- del_dplddt_ddg %>%
  select(pos, ddg)
del_ddg_compl <- t(del_ddg)
del_ddg_heatmap <- matrix(del_ddg_compl[-1,], ncol = 187)
colnames(del_ddg_heatmap) <- del_ddg_compl[1,]
draw(Heatmap(del_ddg_heatmap,
             cluster_rows = FALSE, 
             cluster_columns = FALSE,
             colorRamp2(c(-10, 0, 10), c("darkgreen", "white", "darkorange")),
             heatmap_legend_param = list(direction = "horizontal",
                                         title = expression(Delta*Delta*"G"),
                                         title_position = "topcenter", 
                                         legend_width = unit(6, "cm"),
                                         labels_gp = gpar(fontsize = 12),
                                         title_gp = gpar(fontsize = 14, fontface = "bold"),
                                         border = "black", 
                                         at = c(-10, -5, 0, 5, 10),
                                         labels = c("-10", "-5", "0", "5", "10")),
             row_names_gp = gpar(fontsize = 14),
             column_names_gp = gpar(fontsize = 5),
             row_names_side = "left",
             border_gp = gpar(col = "black", lwd = 1),
             rect_gp = gpar(col = "black", lwd = 1),
             na_col = "dimgrey"),
     heatmap_legend_side = "bottom")


#Fig. S13 (sequence alignment)

#Number of DHFR homologs with deletions per position
seq_align <- ggplot(del_sequence_alignment, aes(x = pos, y = del_freq)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 187, 20)) +
  labs(x = "", y = "# of sequences with deletions") +
  theme_bw()
cpop_score_line <- ggplot() +
  geom_line(del_sequence_alignment, mapping = aes(x = pos, y = score_30, color = "30 °C")) +
  geom_line(del_sequence_alignment, mapping = aes(x = pos, y = score_30_MTX, color = "30 °C + MTX")) +
  geom_line(del_sequence_alignment, mapping = aes(x = pos, y = score_37, color = "37 °C")) +
  geom_line(del_sequence_alignment, mapping = aes(x = pos, y = score_37_MTX, color = "37 °C + MTX")) +
  scale_x_continuous(breaks = seq(0, 187, 20)) +
  scale_y_continuous(limits = c(-0.3, 1.1), breaks = c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  scale_color_manual(values = c("orange", "limegreen", "magenta3", "dodgerblue2")) + 
  labs(x = "Position", y = "CPOP score", colour = "") +
  theme_bw()
print(seq_align/cpop_score_line)

#Fig. S14 (ESM1b ROC curves)

#ESM1b ROC curves for insertions
ins_thres <- cpop_data_ROC_ins %>%
  filter(mut_type == "ins") %>%
  merge(ins_esm1b, by = "pos") %>%
  mutate(cutoff_30 = ifelse(score_30 > 0.5, 1, 0),
         cutoff_30_MTX = ifelse(score_30_MTX > 0.5, 1, 0),
         cutoff_37 = ifelse(score_37 > 0.5, 1, 0),
         cutoff_37_MTX = ifelse(score_37_MTX > 0.5, 1, 0))
ggroc(list(esm1b1b_30 = roc(ins_thres$cutoff_30, ins_thres$esm1b_score),
           esm1b_30_MTX = roc(ins_thres$cutoff_30_MTX, ins_thres$esm1b_score),
           esm1b_37 = roc(ins_thres$cutoff_37, ins_thres$esm1b_score),
           esm1b_37_MTX = roc(ins_thres$cutoff_37_MTX, ins_thres$esm1b_score)),
      legacy.axes = TRUE, linewidth = 2, aes = c("color", "linetype")) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", labels = c(bquote(esm1b1b: "30 °C " (AUC == .(round(auc(ins_thres$cutoff_30, ins_thres$esm1b_score), 2)))),
                                       bquote(esm1b1b: "30 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_30_MTX, ins_thres$esm1b_score), 2)))),
                                       bquote(esm1b1b: "37 °C " (AUC == .(round(auc(ins_thres$cutoff_37, ins_thres$esm1b_score), 2)))),
                                       bquote(esm1b1b: "37 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_37_MTX, ins_thres$esm1b_score), 2))))),
                        values = c(1, 1, 1, 1)) +
  scale_color_manual("", labels = c(c(bquote(esm1b1b: "30 °C " (AUC == .(round(auc(ins_thres$cutoff_30, ins_thres$esm1b_score), 2)))),
                                      bquote(esm1b1b: "30 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_30_MTX, ins_thres$esm1b_score), 2)))),
                                      bquote(esm1b1b: "37 °C " (AUC == .(round(auc(ins_thres$cutoff_37, ins_thres$esm1b_score), 2)))),
                                      bquote(esm1b1b: "37 °C + MTX" (AUC == .(round(auc(ins_thres$cutoff_37_MTX, ins_thres$esm1b_score), 2)))))),
                     values = c("orange", "limegreen", "magenta3", "dodgerblue2")) +
  geom_abline(linetype = "dashed", linewidth = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0, aspect.ratio = 1)

#ESM1b ROC curves for deletions
del_thres <- cpop_data_ROC_del %>%
  filter(mut_type == "del") %>%
  merge(del_esm1b, by = "pos") %>%
  mutate(cutoff_30 = ifelse(score_30 > 0.5, 1, 0),
         cutoff_30_MTX = ifelse(score_30_MTX > 0.5, 1, 0),
         cutoff_37 = ifelse(score_37 > 0.5, 1, 0),
         cutoff_37_MTX = ifelse(score_37_MTX > 0.5, 1, 0))
ggroc(list(esm1b1b_30 = roc(del_thres$cutoff_30, del_thres$esm1b_score),
           esm1b_30_MTX = roc(del_thres$cutoff_30_MTX, del_thres$esm1b_score),
           esm1b_37 = roc(del_thres$cutoff_37, del_thres$esm1b_score),
           esm1b_37_MTX = roc(del_thres$cutoff_37_MTX, del_thres$esm1b_score)),
      legacy.axes = TRUE, linewidth = 2, aes = c("color", "linetype")) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  scale_linetype_manual("", labels = c(bquote(ESM1b: "30 °C " (AUC == .(round(auc(del_thres$cutoff_30, del_thres$esm1b_score), 2)))),
                                       bquote(ESM1b: "30 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_30_MTX, del_thres$esm1b_score), 2)))),
                                       bquote(ESM1b: "37 °C " (AUC == .(round(auc(del_thres$cutoff_37, del_thres$esm1b_score), 2)))),
                                       bquote(ESM1b: "37 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_37_MTX, del_thres$esm1b_score), 2))))),
                        values = c(1, 1, 1, 1)) +
  scale_color_manual("", labels = c(c(bquote(ESM1b: "30 °C " (AUC == .(round(auc(del_thres$cutoff_30, del_thres$esm1b_score), 2)))),
                                      bquote(ESM1b: "30 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_30_MTX, del_thres$esm1b_score), 2)))),
                                      bquote(ESM1b: "37 °C " (AUC == .(round(auc(del_thres$cutoff_37, del_thres$esm1b_score), 2)))),
                                      bquote(ESM1b: "37 °C + MTX" (AUC == .(round(auc(del_thres$cutoff_37_MTX, del_thres$esm1b_score), 2)))))),
                     values = c("orange", "limegreen", "magenta3", "dodgerblue2")) +
  geom_abline(linetype = "dashed", linewidth = 2) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text = element_text(size = 14), legend.text.align = 0, aspect.ratio = 1)
