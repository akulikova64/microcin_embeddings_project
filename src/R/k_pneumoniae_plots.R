# K pneumoniae data:

library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)
library(ggridges)
library(stringr)
library(ggrepel)
library(ggpmisc)
library(gridExtra)
library(grid)

#in K pneumoniae:
# E492: [257-566]
# G492: [4964-5231]

# loading data:
microcin_dist <- read.csv(file = "../../microcin_files/distance_bw_microcins.csv", header=TRUE, sep=",")
E492 <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_E492_sp_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
G492 <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_G492_tr_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
H47 <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_H47_sp_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
I47 <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_I47_tr_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
L <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_L_tr_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
M <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_M_tr_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
N <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_N_tr_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
PDI <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_PDI_tr_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
S <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_S_tr_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")
V <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_V_sp_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")


# making a heat map between the microcins:
microcin_dist2 <- microcin_dist %>%
  mutate(microcin1 = str_sub(microcin1, 0, -4)) %>%
  mutate(microcin2 = str_sub(microcin2, 0, -4)) %>%
  mutate(distance_label = format(round(distance, 2), nsmall = 2))

microcin_dist_plot <- microcin_dist2 %>%
  ggplot(aes(x = microcin1, y = microcin2, fill = distance)) +
  scale_fill_gradient(low = "white", #old low: #e0dede
                      high = "#1b1c3b") + #old high: #31336b
  geom_tile() +
  labs(x = "",
       y = "") +
  ggtitle("Distances between known microcins") +
  geom_text(aes(
    label = distance_label),
    color = "white") +
  scale_x_discrete(
    expand = c(0,0)
  ) +
  scale_y_discrete(
    expand = c(0,0)
  ) +
  theme_bw() +
  theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white"),
      axis.text = element_text(color = "black", size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_line(color = "grey92", size=0.5)
      ) 

microcin_dist_plot

ggsave(filename = "../../analysis/figures/distance_bw_microcins.png", plot = microcin_dist_plot, width = 6, height = 5)


# lets now make a table of microcin distances.

microcin_dist_table <- microcin_dist2 %>%
  filter(microcin1 == "G492"| microcin1 == "E492") %>%
  filter(distance != 0.000) %>%
  distinct(distance, .keep_all = TRUE) %>%
  select(-distance_label)

grid.table(microcin_dist_table)

#-----------------------------------------------------------------------------------
# k pneumomiae plots
#-----------------------------------------------------------------------------------

microcin_dist <- microcin_dist %>%
  filter(distance != 0.0000)

# # point plot of distances
# distance_E492 <- E492 %>%
# ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
#   geom_point()
# 
# distance_E492
# 
# density_E492 <- E492 %>%
#   ggplot(aes(x = distance)) +
#   geom_density(fill = "lightblue") +
#   scale_x_continuous(
#     name = "distance",
#     limits = c(0.0, 11.0),
#     breaks = seq(0.0, 11.0, by = 1),
#     expand = c(0, 0)) + 
#   theme_cowplot()
# 
# density_E492 
# 
# mean_E492 <- mean(E492$distance)  
# mean_E492



E492_new <- E492 %>%
  mutate(microcin = "E492")
G492_new <- G492 %>%
  mutate(microcin = "G492")
H47_new <- H47 %>%
  mutate(microcin = "H47")
I47_new <- I47 %>%
  mutate(microcin = "I47")
L_new <- L %>%
  mutate(microcin = "L")
M_new <- M %>%
  mutate(microcin = "M")
N_new <- N %>%
  mutate(microcin = "N")
PDI_new <- PDI %>%
  mutate(microcin = "PDI")
S_new <- S %>%
  mutate(microcin = "S")
V_new <- V %>%
  mutate(microcin = "V")

all_data <- rbind(E492_new, G492_new, H47_new, I47_new, L_new, M_new, N_new, PDI_new, S_new, V_new)

# make ridgeline plot of distance density for all microcins 

all_density <- all_data %>%
  mutate(microcin = fct_reorder(microcin, distance, mean)) %>%
  ggplot(aes(x = distance, y = microcin)) +
  geom_density_ridges(fill = "cornsilk3", 
                      color = "grey25",
                      alpha = 0.6,
                      rel_min_height=.001,
                      scale = 1) +
  scale_x_continuous(
    name = "distance",
    limits = c(3.0, 10.5),
    breaks = seq(3.0, 10.5, by = 1),
    expand = c(0, 0)) + 
  theme_cowplot() +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

all_density

ggsave(filename = "../../analysis/figures/k_pneumoniae/all_dist_density.png", plot = all_density, width = 6, height = 4)

#---------------------------------------------------------------------------------------------------------
# Comparing distribution of distances between ORF-microcins and microcins-microcins
#---------------------------------------------------------------------------------------------------------

all_data_labeled <- all_data %>%
  select(c(distance, microcin)) %>%
  mutate(group = "ORF-microcin")

microcin_labeled <- microcin_dist %>%
  select(distance) %>%
  mutate(
    group = "microcin-microcin",
    microcin = NA)

box_data <- rbind(all_data_labeled, microcin_labeled)

boxplot <- box_data %>%
  ggplot(aes(y = distance, x = group, fill = group)) +
  geom_boxplot(alpha = 0.5) +
  xlab("") +
  scale_y_continuous(
    name = "distance",
    limits = c(0.0, 12),
    breaks = seq(0.0, 10, by = 2.5),
    expand = c(0, 0)) +
  theme_cowplot(14)  +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14)) +
  scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.minor.y = element_line(color = "grey92", size=0.5))

boxplot

ggsave(filename = "../../analysis/figures/k_pneumoniae/micro_ORF_boxplot.png", plot = boxplot, width = 5, height = 5)

mean_microcins <- mean(microcin_dist$distance)
mean_all_data <- mean(all_data$distance)

#-----------------------------------------------------------------------------
### density plot

two_density_crop <- box_data %>%
  ggplot(aes(x = distance, fill = group)) +
  geom_density(
    alpha = 0.5,
    aes(y = after_stat(count))) +
  #facet_wrap(vars(microcin)) +
  scale_x_continuous(
    name = "distance",
    limits = c(0.0, 5.5),
    breaks = seq(0.0, 5.5, by = 1),
    expand = c(0, 0)) + 
  scale_y_continuous(
    name = "count",
    limits = c(0.0, 250), 
    breaks = seq(0, 240, by = 20),
    expand = expansion(add = c(0,0)))+
  theme_cowplot(14)  +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "black", size = 14)) +
  scale_fill_brewer(palette="Dark2") +
  labs(fill = "") +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_line(color = "grey92", size=0.5),
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

two_density_crop

#ggsave(filename = "../../analysis/figures/k_pneumoniae/two_density_crop.png", plot = two_density_crop, width = 9, height = 7)

two_density <- box_data %>%
  ggplot(aes(x = distance, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(
    name = "distance",
    limits = c(0.0, 10.5),
    breaks = seq(0.0, 10.5, by = 1),
    expand = c(0, 0)) + 
  scale_y_continuous(
    expand = c(0, 0)) + 
  theme_cowplot(14)  +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "black", size = 14)) +
  scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

two_density

#ggsave(filename = "../../analysis/figures/k_pneumoniae/two_density.png", plot = two_density, width = 8, height = 8)

#-------------------------------------------------------------------------------------
# now looking at point plot of distances

# genome annotations:
k_pneum <- read.csv(file = "../../genome_annotation_csvs/k_pneumoniae_cluster_annotation.csv", header=TRUE, sep=",")

k_pneum2 <- k_pneum %>%
  select(c(orf_location, is_microcin)) %>%
  mutate(label = ifelse(is_microcin == "non-microcin", "black", "red"))

E492_2 <- inner_join(k_pneum2, E492)


ordered_E492 <- E492_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_E492 <- ordered_E492 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_E492, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("E492 vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_E492

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_E492.png", plot = distance_E492, width = 6, height = 6)


#------repeating for G492
G492_2 <- inner_join(k_pneum2, G492)


ordered_G492 <- G492_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_G492 <- ordered_G492 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0, NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_G492, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("G492 vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_G492

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_G492.png", plot = distance_G492, width = 6, height = 6)

#------repeating for H47
H47_2 <- inner_join(k_pneum2, H47)


ordered_H47 <- H47_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_H47 <- ordered_H47 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_H47, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("H47 vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_H47

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_H47.png", plot = distance_H47, width = 6, height = 6)

#------repeating for I47 
I47_2 <- inner_join(k_pneum2, I47)


ordered_I47 <- I47_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_I47 <- ordered_I47 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_I47, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("I47 vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_I47

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_I47.png", plot = distance_I47, width = 6, height = 6)

#------repeating for L 
L_2 <- inner_join(k_pneum2, L)


ordered_L <- L_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_L <- ordered_L %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_L, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("L vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_L

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_L.png", plot = distance_L, width = 6, height = 6)

#------repeating for M
M_2 <- inner_join(k_pneum2, M)


ordered_M <- M_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_M <- ordered_M %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_M, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("M vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_M

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_M.png", plot = distance_M, width = 6, height = 6)

#------repeating for N
N_2 <- inner_join(k_pneum2, N)


ordered_N <- N_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_N <- ordered_N %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_N, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("N vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_N

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_N.png", plot = distance_N, width = 6, height = 6)

#------repeating for PDI 
PDI_2 <- inner_join(k_pneum2, PDI)


ordered_PDI <- PDI_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_PDI <- ordered_PDI %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_PDI, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("PDI vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_PDI

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_PDI.png", plot = distance_PDI, width = 6, height = 6)

#------repeating for S
S_2 <- inner_join(k_pneum2, S)


ordered_S <- S_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_S <- ordered_S %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_S, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("S vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_S

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_S.png", plot = distance_S, width = 6, height = 6)

#------repeating for V
V_2 <- inner_join(k_pneum2, V)


ordered_V <- V_2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_V <- ordered_V %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_V, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("V vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_V

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_V.png", plot = distance_V, width = 6, height = 6)

#----------------------------------
#-------using average embedding----
#----------------------------------
# genome annotations:
k_pneum <- read.csv(file = "../../genome_annotation_csvs/k_pneumoniae_cluster_annotation.csv", header=TRUE, sep=",")
average <- read.csv(file = "../../distance_csvs/kpneumoniae_cluster/distances_average_vs_kpneumoniae_cluster.csv", header=TRUE, sep=",")

k_pneum2 <- k_pneum %>%
  select(c(orf_location, is_microcin)) %>%
  mutate(label = ifelse(is_microcin == "non-microcin", "black", "red"))

average2 <- inner_join(k_pneum2, average)

ordered_ave <- average2 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_ave <- ordered_ave %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance, color = label)) +
  geom_point() +
  scale_color_manual(values = c("black", "red")) +
  scale_y_discrete(
    name = "location of ORF in genome") +
  scale_x_continuous(
    limits = c(0,NA),
    breaks = seq(0, 10, by = 1),
    expand = c(0, 0)) + 
  geom_text_repel(
    data = subset(ordered_ave, is_microcin != "non-microcin"),
    aes(
      label = is_microcin),
    max.overlaps = Inf,
    box.padding = 0.7,
    color = "red") +
  theme_classic() +
  ggtitle("Average emb. vs K.pneumoniae cluster") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "grey92", size = 0.5),
    panel.grid.minor.x = element_line(color = "grey92", size = 0.5))

distance_ave

ggsave(filename = "../../analysis/figures/k_pneumoniae/dist_point_average.png", plot = distance_ave, width = 6, height = 6)



