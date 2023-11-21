rm(list = ls())#clear the environment to save on driver space and run faster
setwd("C:/Users/mafata/Documents/GitHub/hplc_scanning/assets")
library("ncdf4") # for accessing information in a cdf file
library("tidyverse") # to wrangle data frames
library("glue") 
library("ggplot2") # to plot the spectra
library("scales") # for formating numerical data in plots
# Create a data frame with all the .aia file folder names and paths
hplc_wines <-
  data.frame(
    filename = list.files(
      "C:/Users/mafata/Desktop/WORK/Collaborative Work/HPLC scanning/CDF files"
    )
  )
hplc_wines <- hplc_wines  %>%
  mutate(
    filepath = paste0(
      "C:/Users/mafata/Desktop/WORK/Collaborative Work/HPLC scanning/CDF files/",
      filename
    )
  )
# 
for (i in 1:length(hplc_wines$filepath)){
  hplc_wines <- hplc_wines  %>%
    mutate(
      samples = str_sub (hplc_wines$filename, end = -8),
      repeats = str_sub (hplc_wines$filename, end = -5),
      uv_280 = paste0(filepath,"/SIGNAL01.cdf"),
      uv_320 = paste0(filepath,"/SIGNAL02.cdf"),
      uv_360 = paste0(filepath,"/SIGNAL03.cdf"),
      uv_420 = paste0(filepath,"/SIGNAL04.cdf"),
      fluo = paste0(filepath,"/SIGNAL05.cdf") # FLD , Ex=280, Em=320
    )
}
# Fetch the data from each signal file
datasets = c("uv_280","uv_320","uv_360","uv_420","fluo")
for (dataset in datasets){
dataset_list = list()
for (i in 1:length(hplc_wines[[dataset]])){
  dataset_file <- 
    nc_open(
      hplc_wines[[dataset]][i],
      write = FALSE,
      readunlim = FALSE,
      verbose = FALSE,
      auto_GMT = TRUE,
      suppress_dimvals = FALSE,
      return_on_error = FALSE
    )
  dataset_list = append(dataset_list, list(dataset_file))
}
names(dataset_list) = hplc_wines$filename
if (dataset == "uv_280"){uv280_list = dataset_list}
else if (dataset == "uv_320"){uv320_list = dataset_list}
else if (dataset == "uv_360"){uv360_list = dataset_list}
else if (dataset == "uv_420"){uv420_list = dataset_list}
else {fluo_list = dataset_list}
}

# fetch the spectra for each signal into separate data frames
uv_280_spectra = list()
for (i in 1:length(uv280_list)){
  # Add verbosity to the script
  print(glue(". . . generating sample number {i} uv 280 nm spectra"))
  peak_retention_time <-as.data.frame(ncvar_get(
  uv280_list[[i]],
  uv280_list[[i]]$var$peak_retention_time))
  colnames(peak_retention_time) = "peak_retention_time"
  peak_retention_time <-format(round(peak_retention_time, 0), nsmall = 0)
  peak_area <-as.data.frame(ncvar_get(
  uv280_list[[i]],
  uv280_list[[i]]$var$peak_area))
  colnames(peak_area) = hplc_wines$repeats[i]
  uv_280_spectra <- append(uv_280_spectra,
                              list(name = c(peak_retention_time, peak_area)),
                              after = length(uv_280_spectra))
}
names(uv_280_spectra) = hplc_wines$filename

# retention time deviations. Correction through merging.
# create a merged data frame block for MFA analysis
merged_uv280_spectra <- full_join(x=as.data.frame(uv_280_spectra[[1]]), 
                          y=as.data.frame(uv_280_spectra[[2]]), 
                          by="peak_retention_time")
for (i in 3:length(uv_280_spectra)){
merged_uv280_spectra <- full_join(x=as.data.frame(merged_uv280_spectra), 
                                  y=as.data.frame(uv_280_spectra[[i]]), 
                                  by="peak_retention_time")
}

# print the spectra for the 280 spectra
plt <- merged_uv280_spectra %>% pivot_longer(!peak_retention_time, 
                                             names_to = "samples",
                                             values_to = "peak_area",
                                             values_drop_na = TRUE)
plt$peak_retention_time = as.numeric(plt$peak_retention_time)/60
uv_280 <- ggplot(
  plt,
  aes(x = peak_retention_time, y = peak_area, colour = samples)
) + geom_line(linewidth = 0.1) + theme_bw(base_family = "Arial") + theme(
  axis.text.x = element_text(
    angle = 0,
    hjust = 0.5,
    vjust = 0.5
  ),
  axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
  axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
  axis.ticks.length.x = unit(0, "cm"),
  axis.ticks.length.y = unit(0, "cm"),
  legend.key.height = unit(0.5, "cm"),
  legend.key.width = unit(1, "cm"),
  legend.direction = "horizontal",
  legend.position = "none",
  legend.box.background = element_rect(
    fill = 'white',
    colour = 'black',
    linewidth = 0.1
  ),
  legend.text = element_text(size = 4),
  legend.title = element_text(size = 4)
) +
  scale_x_continuous(breaks = seq(
  from = 0,
  to = 30,
  by = 1
)) +
  scale_y_continuous(
  breaks = seq(
    from = 0,
    to = 650,
    length = 10
  ),
  labels = scales::label_comma()
)

ggsave("uv_280_spectral_overlay.svg",
       plot = uv_280,
       width = 30,
       height = 12,
       units = "cm")


# fetch the 320 nm spectra for each signal into separate data frames
uv_320_spectra = list()
for (i in 1:length(uv320_list)){
  # Add verbosity to the script
  print(glue(". . . generating sample number {i} uv 320 nm spectra"))
  peak_retention_time <-as.data.frame(ncvar_get(
    uv320_list[[i]],
    uv320_list[[i]]$var$peak_retention_time))
  colnames(peak_retention_time) = "peak_retention_time"
  peak_retention_time <-format(round(peak_retention_time, 0), nsmall = 0)
  peak_area <-as.data.frame(ncvar_get(
    uv320_list[[i]],
    uv320_list[[i]]$var$peak_area))
  colnames(peak_area) = hplc_wines$repeats[i]
  uv_320_spectra <- append(uv_320_spectra,
                           list(name = c(peak_retention_time, peak_area)),
                           after = length(uv_320_spectra))
}
names(uv_320_spectra) = hplc_wines$filename

# retention time deviations. Correction through merging.
# create a merged data frame block for MFA analysis
merged_uv320_spectra <- full_join(x=as.data.frame(uv_320_spectra[[1]]), 
                                  y=as.data.frame(uv_320_spectra[[2]]), 
                                  by="peak_retention_time")
for (i in 3:length(uv_320_spectra)){
  merged_uv320_spectra <- full_join(x=as.data.frame(merged_uv320_spectra), 
                                    y=as.data.frame(uv_320_spectra[[i]]), 
                                    by="peak_retention_time")
}

# print the spectra for the 320 spectra

plt2 <- merged_uv320_spectra %>% pivot_longer(!peak_retention_time, 
                                             names_to = "samples",
                                             values_to = "peak_area",
                                             values_drop_na = TRUE)
plt2$peak_retention_time = as.numeric(plt2$peak_retention_time)/60

uv_320 <- ggplot(
  plt2,
  aes(x = peak_retention_time, y = peak_area, colour = samples)
) + geom_line(linewidth = 0.1) + theme_bw(base_family = "Arial") + theme(
  axis.text.x = element_text(
    angle = 0,
    hjust = 0.5,
    vjust = 0.5
  ),
  axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
  axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
  axis.ticks.length.x = unit(0, "cm"),
  axis.ticks.length.y = unit(0, "cm"),
  legend.key.height = unit(0.5, "cm"),
  legend.key.width = unit(1, "cm"),
  legend.direction = "horizontal",
  legend.position = "none",
  legend.box.background = element_rect(
    fill = 'white',
    colour = 'black',
    linewidth = 0.1
  ),
  legend.text = element_text(size = 4),
  legend.title = element_text(size = 4)
) +
  scale_x_continuous(breaks = seq(
    from = 0,
    to = 30,
    by = 1
  )) +
  scale_y_continuous(
    breaks = seq(
      from = 0,
      to = ceiling(max(plt2$peak_area)),
      length = 10
    ),
    labels = scales::label_comma()
  )
ggsave("uv_320_spectral_overlay.svg",
       plot = uv_320,
       width = 30,
       height = 12,
       units = "cm")





# run the MFA
library("FactoMineR")
library("factoextra")
df1 <- merged_uv280_spectra
row.names(df1) <- df1$peak_retention_time
df1 <- t(df1)
df1 <- as.data.frame(df1[c(-1),])
df1$samples <- row.names(df1)
  
df2 <- merged_uv320_spectra
row.names(df2) <- df2$peak_retention_time
df2 <- t(df2)
df2 <- as.data.frame(df2[c(-1),])
df2$samples <- row.names(df2)

merged_df <- full_join(x=df1,
                       y=df2,
                       by="samples")
merged_df <- merged_df[,c(-1234)]
merged_df[] <- lapply(df, as.numeric)
  
mfa_plot <- MFA(
  merged_df,
  group = c(770,463),
  type = c(rep("s", 2)),
  ncp = 40,
  name.group = c("CB",
                 "SB"),
  graph = TRUE)
fviz_mfa_ind(mfa_plot,
             geom = c("point", "text"),
             axes = c(1,2)
)
