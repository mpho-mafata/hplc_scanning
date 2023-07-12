# HPLC scanning data
{
  rm(list = ls())#clear the environment to save on driver space and run faster
  setwd("C:/Users/mafata/Desktop/WORK/Collaborative Work/Cody")
  library("ncdf4")
  library("tidyverse")
  library("hrbrthemes")
  
  # DAD 280nm
  signal_1_file <- 'C:/Users/mafata/Desktop/WORK/Collaborative Work/HPLC scanning/CDF files/AVN_CB_T0_R1.aia/SIGNAL01.cdf'
  signal_1 <- nc_open(signal_1_file, write=FALSE, readunlim=FALSE, verbose=FALSE, 
                      auto_GMT=TRUE, suppress_dimvals=FALSE, return_on_error=FALSE ) 
  peak_retention_time <- ncvar_get( signal_1, signal_1$var$peak_retention_time)
  peak_retention_time <- as.data.frame(peak_retention_time)
  peak_area <- ncvar_get( signal_1, signal_1$var$peak_area)
  peak_area <- as.data.frame(peak_area)
  signal_1_df <- as.data.frame(c(peak_retention_time, peak_area))
  signal_1_df$peak_retention_time <- (signal_1_df$peak_retention_time)/60 
  
  # DAD 320nm
  signal_2_file <- 'C:/Users/mafata/Desktop/WORK/Collaborative Work/HPLC scanning/CDF files/AVN_CB_T0_R1.aia/SIGNAL02.cdf'
  signal_2 <- nc_open(signal_2_file, write=FALSE, readunlim=FALSE, verbose=FALSE, 
                      auto_GMT=TRUE, suppress_dimvals=FALSE, return_on_error=FALSE ) 
  peak_retention_time <- ncvar_get( signal_2, signal_2$var$peak_retention_time)
  peak_retention_time <- as.data.frame(peak_retention_time)
  peak_area <- ncvar_get( signal_2, signal_2$var$peak_area)
  peak_area <- as.data.frame(peak_area)
  signal_2_df <- as.data.frame(c(peak_retention_time, peak_area))
  signal_2_df$peak_retention_time <- (signal_2_df$peak_retention_time)/60 
  
  # DAD 360nm
  signal_3_file <- 'C:/Users/mafata/Desktop/WORK/Collaborative Work/HPLC scanning/CDF files/AVN_CB_T0_R1.aia/SIGNAL03.cdf'
  signal_3 <- nc_open(signal_3_file, write=FALSE, readunlim=FALSE, verbose=FALSE, 
                      auto_GMT=TRUE, suppress_dimvals=FALSE, return_on_error=FALSE ) 
  peak_retention_time <- ncvar_get( signal_3, signal_2$var$peak_retention_time)
  peak_retention_time <- as.data.frame(peak_retention_time)
  peak_area <- ncvar_get( signal_3, signal_3$var$peak_area)
  peak_area <- as.data.frame(peak_area)
  signal_3_df <- as.data.frame(c(peak_retention_time, peak_area))
  signal_3_df$peak_retention_time <- (signal_3_df$peak_retention_time)/60 
  
  # DAD 420nm
  signal_4_file <- 'C:/Users/mafata/Desktop/WORK/Collaborative Work/HPLC scanning/CDF files/AVN_CB_T0_R1.aia/SIGNAL04.cdf'
  signal_4 <- nc_open(signal_4_file, write=FALSE, readunlim=FALSE, verbose=FALSE, 
                      auto_GMT=TRUE, suppress_dimvals=FALSE, return_on_error=FALSE ) 
  peak_retention_time <- ncvar_get( signal_4, signal_2$var$peak_retention_time)
  peak_retention_time <- as.data.frame(peak_retention_time)
  peak_area <- ncvar_get( signal_4, signal_4$var$peak_area)
  peak_area <- as.data.frame(peak_area)
  signal_4_df <- as.data.frame(c(peak_retention_time, peak_area))
  signal_4_df$peak_retention_time <- (signal_4_df$peak_retention_time)/60 
  
  # FLD , Ex=280, Em=320
  signal_5_file <- 'C:/Users/mafata/Desktop/WORK/Collaborative Work/HPLC scanning/CDF files/AVN_CB_T0_R1.aia/SIGNAL05.cdf'
  signal_5 <- nc_open(signal_5_file, write=FALSE, readunlim=FALSE, verbose=FALSE, 
                      auto_GMT=TRUE, suppress_dimvals=FALSE, return_on_error=FALSE ) 
  peak_retention_time <- ncvar_get( signal_4, signal_2$var$peak_retention_time)
  peak_retention_time <- as.data.frame(peak_retention_time)
  peak_area <- ncvar_get( signal_4, signal_4$var$peak_area)
  peak_area <- as.data.frame(peak_area)
  signal_4_df <- as.data.frame(c(peak_retention_time, peak_area))
  signal_4_df$peak_retention_time <- (signal_4_df$peak_retention_time)/60 
  
  {
    signal_01_linegraph <- ggplot(signal_1_df,
                                  aes(x = peak_retention_time, y = peak_area))+
      geom_line(linewidth = 0.2)+ 
      geom_point(size = 0.2, shape = 21)+
      theme_ipsum(base_family = "Arial",axis_text_size = 9) +
      theme(axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5),
            axis.title.x=element_text(hjust=0.5,vjust=0.5),
            axis.title.y=element_text(hjust=0.5,vjust=0.5),
            axis.ticks.length.x = unit(0, "cm"),
            axis.ticks.length.y = unit(0, "cm"),
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(0.3, "cm"),
            legend.direction = "vertical",
            legend.position = c(0.18, 0.55),
            legend.box.background = element_rect(fill = 'white', colour = 'black', linewidth = 0.1),
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 10)
      )
    # scale_x_continuous(breaks = seq(from = min(get(patent_type)$year), to = max(get(patent_type)$year), by = 2),
    # minor_breaks = 2)+
    # scale_y_continuous(breaks = my_breaks, minor_breaks = my_minor_breaks)+
    # labs(x = "RT(min)", y = glue("intensity"), colour = "organizations", fill = "organizations")
    ggsave(glue("signal_01.jpg"),
           plot = signal_01_linegraph,
           width = 30,
           height = 15,
           units = 'cm')
    browseURL(glue("signal_01.jpg"))
  }
}

# FactomineR PCA plots
{
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

}