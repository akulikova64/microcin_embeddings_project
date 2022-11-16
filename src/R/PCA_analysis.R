library(tidyverse)
library(cowplot)
options(scipen = 999)
library(colorspace)
library(broom)  # for augment(), tidy()

# PCA between 10 microcins and ORFs

# loading data:
embeddings <- read.csv(file = "../../for_PCA/kpneumoniae_PCA_data.csv", header=TRUE, sep=",")

embeddings_filtered <- embeddings %>%
  filter(name != 378,
         name != 354,
         name != 176)

pca_fit <- embeddings_filtered %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  scale() %>%                   # scale to zero mean and unit variance
  prcomp()  

point_plot1 <- pca_fit %>%
  # add PCs to the original dataset
  augment(embeddings_filtered) %>%
  ggplot(aes(.fittedPC1, .fittedPC2)) +
  geom_point(aes(color = group)) +
  labs(x = "PC1", y = "PC2") +
  theme_bw() +
  scale_color_brewer(palette = "Set1")

point_plot1

point_plot2 <- pca_fit %>%
  # add PCs to the original dataset
  augment(embeddings_filtered) %>%
  ggplot(aes(.fittedPC2, .fittedPC3)) +
  geom_point(aes(color = group)) +
  labs(x = "PC2", y = "PC3") +
  theme_bw() +
  scale_color_brewer(palette = "Set1")

point_plot2

# rotation matrix
arrow_style <- arrow(
  angle = 20, length = grid::unit(8, "pt"),
  ends = "first", type = "closed"
)

rotation_matrix <- pca_fit %>%
  # extract rotation matrix
  tidy(matrix = "rotation") %>%
  pivot_wider(
    names_from = "PC", values_from = "value",
    names_prefix = "PC"
  ) %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(
    xend = 0, yend = 0,
    arrow = arrow_style
  ) +
  geom_text(aes(label = column), hjust = 1) +
  coord_fixed()

rotation_matrix

#------------------------------------------------------------------------------------------
# make variance explained plot
variance <- pca_fit %>%
  # extract eigenvalues
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) + 
  geom_col(fill = "darkslategrey", alpha = 0.8) + 
  scale_x_continuous(
    # create one axis tick per PC
    limits = c(0, 21),
    breaks = 1:20,
    expand = expansion(add = c(0, 0))
  ) +
  scale_y_continuous(
    name = "variance explained",
    # format y axis ticks as percent values
    label = scales::label_percent(accuracy = 1),
    expand = expansion(add = c(0, 0.005))
  ) +
  theme_classic()
variance
