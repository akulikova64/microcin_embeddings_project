library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)
library(ggridges)

# loading data:
microcin_dist <- read.csv(file = "../../microcin_files/distance_bw_microcins.csv", header=TRUE, sep=",")
E492 <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_E492_sp_vs_ecoli_s88.csv", header=TRUE, sep=",")
G492 <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_G492_tr_vs_ecoli_s88.csv", header=TRUE, sep=",")
H47 <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_H47_sp_vs_ecoli_s88.csv", header=TRUE, sep=",")
I47 <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_I47_tr_vs_ecoli_s88.csv", header=TRUE, sep=",")
L <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_L_tr_vs_ecoli_s88.csv", header=TRUE, sep=",")
M <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_M_tr_vs_ecoli_s88.csv", header=TRUE, sep=",")
N <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_N_tr_vs_ecoli_s88.csv", header=TRUE, sep=",")
PDI <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_PDI_tr_vs_ecoli_s88.csv", header=TRUE, sep=",")
S <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_S_tr_vs_ecoli_s88.csv", header=TRUE, sep=",")
V <- read.csv(file = "../../distance_csvs/ecoli_s88/distances_V_sp_vs_ecoli_s88.csv", header=TRUE, sep=",")

microcin_dist <- microcin_dist %>%
  filter(distance != 0.0000)

# point plot of distances
distance_E492 <- E492 %>%
  mutate()
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point()

distance_E492
  
density_E492 <- E492 %>%
  ggplot(aes(x = distance)) +
  geom_density(fill = "lightblue") +
  scale_x_continuous(
    name = "distance",
    limits = c(0.0, 11.0),
    breaks = seq(0.0, 11.0, by = 1),
    expand = c(0, 0)) + 
  theme_cowplot()
  
density_E492 

mean_E492 <- mean(E492$distance)  
mean_E492


  
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

#ggsave(filename = "../../analysis/figures/all_dist_density.png", plot = all_density, width = 8, height = 6)

#---------------------------------------------------------------------------------------------------------
# Comparing distribution of distances between ORF-microcins and microcins-microcins
#---------------------------------------------------------------------------------------------------------

all_data_labeled <- all_data %>%
  select(distance) %>%
  mutate(group = "ORF-microcin")

microcin_labeled <- microcin_dist %>%
  select(distance) %>%
  mutate(group = "microcin-microcin")

box_data <- rbind(all_data_labeled, microcin_labeled)

boxplot <- box_data %>%
  ggplot(aes(y = distance, x = group, fill = group)) +
  geom_boxplot(alpha = 0.5) +
  xlab("") +
  ylab("distance") +
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

#ggsave(filename = "../../analysis/figures/micro_ORF_boxplot.png", plot = boxplot, width = 8, height = 8)

mean_microcins <- mean(microcin_dist$distance)
mean_all_data <- mean(all_data$distance)

#-----------------------------------------------------------------------------
### density plot

two_density_crop <- box_data %>%
  ggplot(aes(x = distance, fill = group)) +
  geom_density(
    alpha = 0.5,
    aes(y = after_stat(count))) +
  scale_x_continuous(
    name = "distance",
    limits = c(0.0, 5.5),
    breaks = seq(0.0, 5.5, by = 1),
    expand = c(0, 0)) + 
  scale_y_continuous(
    name = "count",
    limits = c(0.0, 1000), 
    breaks = seq(0, 1000, by = 250),
    expand = c(0, 0))+
  theme_cowplot(14)  +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "black", size = 14)) +
  scale_fill_brewer(palette="Dark2") +
  theme(
    panel.grid.major.x = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

two_density_crop

#ggsave(filename = "../../analysis/figures/two_density_crop.png", plot = two_density_crop, width = 8, height = 8)

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

#ggsave(filename = "../../analysis/figures/two_density.png", plot = two_density, width = 8, height = 8)

#-------------------------------------------------------------------------------------
# now looking at point plot of distances

ordered_E492 <- E492 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_E492 <- ordered_E492 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("E492 vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_E492

#ggsave(filename = "../../analysis/figures/dist_point_E492.png", plot = distance_E492, width = 8, height = 8)


#------repeating for G492
ordered_G492 <- G492 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_G492 <- ordered_G492 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() + 
  ggtitle("G492 vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_G492

#ggsave(filename = "../../analysis/figures/dist_point_G492.png", plot = distance_G492, width = 8, height = 8)

#------repeating for H47
ordered_H47 <- H47 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_H47 <- ordered_H47 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("H47 vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_H47

ggsave(filename = "../../analysis/figures/dist_point_H47.png", plot = distance_H47, width = 8, height = 8)

#------repeating for I47 
ordered_I47 <- I47 %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_I47 <- ordered_I47 %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("I47 vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_I47

ggsave(filename = "../../analysis/figures/dist_point_I47.png", plot = distance_I47, width = 8, height = 8)

#------repeating for L 
ordered_L <- L %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_L <- ordered_L %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("L vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_L

ggsave(filename = "../../analysis/figures/dist_point_L.png", plot = distance_L, width = 8, height = 8)

#------repeating for M
ordered_M <- M %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_M <- ordered_M %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("M vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_M

ggsave(filename = "../../analysis/figures/dist_point_M.png", plot = distance_M, width = 8, height = 8)

#------repeating for N
ordered_N <- N %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_N <- ordered_N %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("N vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_N

ggsave(filename = "../../analysis/figures/dist_point_N.png", plot = distance_N, width = 8, height = 8)

#------repeating for PDI 
ordered_PDI <- PDI %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_PDI <- ordered_PDI %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("PDI vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_PDI

ggsave(filename = "../../analysis/figures/dist_point_PDI.png", plot = distance_PDI, width = 8, height = 8)

#------repeating for S
ordered_S <- S %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_S <- ordered_PDI %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("S vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_S

ggsave(filename = "../../analysis/figures/dist_point_S.png", plot = distance_S, width = 8, height = 8)

#------repeating for V
ordered_V <- V %>%
  slice_min(order_by = distance, n = 50)

# point plot of distances
distance_V <- ordered_V %>%
  ggplot(aes(y = fct_reorder(orf_location, distance), x = distance)) +
  geom_point() +
  scale_y_discrete(
    name = "location of ORF in genome") +
  theme_classic() +
  ggtitle("V vs Ecoli S88") +
  theme(
    panel.grid.major = element_line(color = "grey92", size=0.5),
    panel.grid.minor.x = element_line(color = "grey92", size=0.5))

distance_V

ggsave(filename = "../../analysis/figures/dist_point_V.png", plot = distance_V, width = 8, height = 8)




#-------
# genome annotations:

e_coli_s88 <- read.csv(file = "../../genome_annotation_csvs/ecoli_s88_annotation.csv", header=TRUE, sep=",")
e_coli_nissle <- read.csv(file = "../../genome_annotation_csvs/ecoli_nissle_annotation.csv", header=TRUE, sep=",")
k_pneum <- read.csv(file = "../../genome_annotation_csvs/k_pneumoniae_cluster_annotation.csv", header=TRUE, sep=",")










