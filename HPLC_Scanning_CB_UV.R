rm(list = ls())
library(readxl)
setwd("C:/Users/mafata/Desktop/WORK/HPLC scanning/Charts/CBUV")
library(FactoMineR) # For producing  PCA
library(tidyverse) # This has ggplot2 and tidy packages as dependencies
library(tidyselect) # To select inside a dplyr pipe
library(factoextra) # Additional visualization commands using fviz
library(glue) # Add python-style f-strings
library(plotly) # To create interactive charts
library(htmlwidgets) # To save interactive charts
library(rlist) # Add the list type
library(ggsci) #colour pallet for scientific plots
library(sjPlot) #for saving the images

#################import data sets##########################
CBUVabsorbance <-
  read_excel(
    "C:/Users/mafata/Desktop/TOTALLY A WORK FILE/HPLC scanning/CB_HPLC.xlsx",
    col_names = TRUE,
    sheet = 'CBUVabsorbance'
  )
CBUVabsorbance <- as.data.frame(CBUVabsorbance)
rownames(CBUVabsorbance) <- CBUVabsorbance$PrimaryID

################### Make a list of the three regions######
mylist <- list("AVN", "CDB", "FRV", "DTK", "KZC", "PDB")

for (farm in mylist) {
  ################filter by farm###########################
  farm_table <- CBUVabsorbance %>%
    filter(Winery == farm)
  
  ##################Generate PCA############################
  farm_pca <- PCA(
    farm_table,
    quali.sup = c(1:3),
    ind.sup = 3,
    ncp = 25,
    graph = FALSE
  )
  ################### scree plots ###########################
  farm_scree <-
    fviz_screeplot(farm_pca, ncp = 25, addlabels = TRUE)
  save_plot(
    glue("{farm}_scree.jpg"),
    fig = farm_scree,
    width = 29.7,
    height = 21
  )
  
  gg_farm_scree <- ggplotly(farm_scree)
  saveWidget(gg_farm_scree, file = glue("gg_{farm}_scree.html"))
  browseURL(glue("gg_{farm}_scree.html"))
  
  ################### Generate the PCA plots ################
  
  ####Graphical plots for individuals, variables, and biplot####
  farm_pca_ind <- fviz_pca_ind(
    farm_pca,
    axes = c(1, 2),
    geom.ind = c("point"),
    habillage = 2,
    repel = TRUE,
    addEllipses = TRUE,
    ellipse.level = 0.95,
    palette = "lancet",#colour palette for the pca
    label = "ind"
  ) + theme(
    text = element_text(size = 5),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.title = element_text(size = 5),
    axis.title.x = element_text(size = 5, hjust = 0.75),
    axis.title.y = element_text(size = 5, hjust = 0.6)
  )
  save_plot(
    glue("{farm}_pca_ind.jpg"),
    fig = farm_pca_ind,
    width = 29.7,
    height = 21
  )
  
  gg_farm_pca_ind <-
    ggplotly(farm_pca_ind, tooltip = c("text"))
  saveWidget(gg_farm_pca_ind, file = glue("{farm}_pca_ind.html"))
  browseURL(glue("{farm}_pca_ind.html"))
  
  ###########clustering of individuals ##########################
  farm_clustering <-
    HCPC(
      farm_pca,
      nb.clust = -1,
      consol = FALSE,
      iter.max = 10,
      min = 3,
      max = 4,
      metric = "euclidean",
      method = "ward",
      order = TRUE,
      graph.scale = "inertia",
      nb.par = 5,
      graph = FALSE,
      description = TRUE
    )
  farm_cluster_plot <-
    plot.HCPC(
      farm_clustering,
      axes = c(1, 2),
      choice = "tree",
      cex = 0.5,
      rect = FALSE,
      draw.tree = TRUE,
      ind.names = TRUE,
      title = NULL,
      tree.barplot = TRUE,
      centers.plot = FALSE
    )
  save_plot(
    glue("{farm}_pca_cluster_ind.jpg"),
    fig = farm_cluster_plot,
    width = 29.7,
    height = 21
  )
  ##########generate the PCA variables plots ###########
  farm_pca_var <- fviz_pca_var(
    farm_pca,
    axes = c(1, 2),
    repel = TRUE,
    geom.ind = c("point"),
    label = "none",
    geom = "point",
    col.var = "contrib"
  ) + theme(
    text = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.position = "bottom",
    plot.title = element_text(size = 8),
    axis.title.x = element_text(size = 8, hjust = 0.75),
    axis.title.y = element_text(size = 8, hjust = 0.6)
  )
  save_plot(
    glue("{farm}_pca_ind.jpg"),
    fig = farm_pca_ind,
    width = 29.7,
    height = 21
  )
  gg_farm_pca_var <-
    ggplotly(farm_pca_var, tooltip = c("x", "y", "text"))
  saveWidget(gg_farm_pca_var, file = glue("{farm}_pca_var.html"))
  browseURL(glue("{farm}_pca_var.html"))
}
