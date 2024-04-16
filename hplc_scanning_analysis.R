#clear the environment to save on driver space and run faster
rm(list = ls())
setwd("/Users/mphomafata/Documents/GitHub/hplc_scanning")

# import the necessary libraries
library("tidyverse") # to wrangle data frames
library("glue") 
library("ggplot2") # to plot the spectra
library("scales") # for formating numerical data in plots
library("FactoMineR")
library("factoextra")

# Read-in the data tables
uv_280_sb <- readxl::read_excel(path = "/Users/mphomafata/Documents/Work_file/Collaborative Work/HPLC scanning/CSV data/HPLC scouting.xlsx",
                                 sheet = "280 nm SB")

uv_280_cb <- readxl::read_excel(path = "/Users/mphomafata/Documents/Work_file/Collaborative Work/HPLC scanning/CSV data/HPLC scouting.xlsx",
                                 sheet = "280 nm CB")
fl_sb <- readxl::read_excel(path = "/Users/mphomafata/Documents/Work_file/Collaborative Work/HPLC scanning/CSV data/HPLC scouting.xlsx",
                                sheet = "FLD SB")
fl_cb <- readxl::read_excel(path = "/Users/mphomafata/Documents/Work_file/Collaborative Work/HPLC scanning/CSV data/HPLC scouting.xlsx",
                                sheet = "FLD CB")

files_list <- list(uv_280_sb, uv_280_cb, fl_sb, fl_cb)

# split the files into their respective samples
uv_280_sb_list <- list()
for (i in 1:(ncol(uv_280_sb))) {
  first <- as.numeric((2 * i) - 1)
  second <- as.numeric(2 * i)
  sample_i <- uv_280_sb[c(first, second)]
  colnames(sample_i)[1] <- c("rt")
  sample_i <- list(sample_i)
  uv_280_sb_list <- append(uv_280_sb_list, sample_i)
}
name_list <- row.names(as.data.frame(t(uv_280_sb)) %>% filter(row_number() %% 2 == 0))
names(uv_280_sb_list) <- name_list

# uv 280 cb
uv_280_cb_list <- list()
for (i in 1:(ncol(uv_280_cb))) {
  first <- as.numeric((2 * i) - 1)
  second <- as.numeric(2 * i)
  sample_i <- uv_280_cb[c(first, second)]
  colnames(sample_i)[1] <- c("rt")
  sample_i <- list(sample_i)
  uv_280_cb_list <- append(uv_280_cb_list, sample_i)
}
name_list <- row.names(as.data.frame(t(uv_280_cb)) %>% filter(row_number() %% 2 == 0))
names(uv_280_cb_list) <- name_list

# fluorescence sb
fl_sb_list <- list()
for (i in 1:(ncol(fl_sb))) {
  first <- as.numeric((2 * i) - 1)
  second <- as.numeric(2 * i)
  sample_i <- fl_sb[c(first, second)]
  colnames(sample_i)[1] <- c("rt")
  sample_i <- list(sample_i)
  fl_sb_list <- append(fl_sb_list, sample_i)
}
name_list <- row.names(as.data.frame(t(fl_sb)) %>% filter(row_number() %% 2 == 0))
names(fl_sb_list) <- name_list

# fluorescence cb
fl_cb_list <- list()
for (i in 1:(ncol(fl_cb))) {
  first <- as.numeric((2 * i) - 1)
  second <- as.numeric(2 * i)
  sample_i <- fl_cb[c(first, second)]
  colnames(sample_i)[1] <- c("rt")
  sample_i <- list(sample_i)
  fl_cb_list <- append(fl_cb_list, sample_i)
}
name_list <- row.names(as.data.frame(t(fl_cb)) %>% filter(row_number() %% 2 == 0))
names(fl_cb_list) <- name_list


# create a merged data frame from the lists
# uv 280 sb
merged_uv280_sb_spectra <- full_join(x = as.data.frame(uv_280_sb_list[[1]]),
                                  y = as.data.frame(uv_280_sb_list[[2]]),
                                  by = "rt")
for (i in 3:length(uv_280_sb_list)) {
  merged_uv280_sb_spectra <-
    full_join(
      x = as.data.frame(merged_uv280_sb_spectra),
      y = as.data.frame(uv_280_sb_list[[i]]),
      by = "rt"
    )
}

# uv 280 cb
merged_uv280_cb_spectra <- full_join(x = as.data.frame(uv_280_cb_list[[1]]),
                                     y = as.data.frame(uv_280_cb_list[[2]]),
                                     by = "rt")
for (i in 3:length(uv_280_cb_list)) {
  merged_uv280_cb_spectra <-
    full_join(
      x = as.data.frame(merged_uv280_cb_spectra),
      y = as.data.frame(uv_280_cb_list[[i]]),
      by = "rt"
    )
}

# fluorescence sb
merged_fl_sb_spectra <- full_join(x = as.data.frame(fl_sb_list[[1]]),
                                  y = as.data.frame(fl_sb_list[[2]]),
                                  by = "rt",
                                  relationship = "many-to-many")
for (i in 3:length(fl_sb_list)) {
  merged_fl_sb_spectra <-
    full_join(
      x = as.data.frame(merged_fl_sb_spectra),
      y = as.data.frame(fl_sb_list[[i]]),
      by = "rt",
      relationship = "many-to-many"
    )
}

# fluorescence cb
merged_fl_cb_spectra <- full_join(x = as.data.frame(fl_cb_list[[1]]),
                                  y = as.data.frame(fl_cb_list[[2]]),
                                  by = "rt")
for (i in 3:length(fl_cb_list)) {
  merged_fl_cb_spectra <-
    full_join(
      x = as.data.frame(merged_fl_cb_spectra),
      y = as.data.frame(fl_cb_list[[i]]),
      by = "rt"
    )
}

# run the MFAs
# uv 280 sb
mfa_plot_uv28_sb <- MFA(
  merged_uv280_sb_spectra[-1],
  group = c(26, 25, 25, 25, 25, 25),
  type = c(rep("s", 6)),
  ncp = 4,
  name.group = c("AVN", "CDB", "DTK", "FRV", "KZN", "PDB"),
  graph = TRUE
)
mfa_groups_uv280_sb <- plot.MFA(mfa_plot_uv28_sb, choix = "group")
ggsave(
  "mfa_groups_uv28_sb.jpg",
  plot = mfa_groups_uv280_sb,
  width = 30,
  height = 15,
  units = 'cm',
  dpi = 300
)
browseURL("mfa_groups_uv28_sb.jpg")

# uv 280 cb
mfa_plot_uv28_cb <- MFA(
  merged_uv280_cb_spectra[-1],
  group = c(26, 25, 25, 25, 25, 25),
  type = c(rep("s", 6)),
  ncp = 4,
  name.group = c("AVN", "CDB", "DTK", "FRV", "KZN", "PDB"),
  graph = TRUE
)
mfa_groups_uv280_cb <- plot.MFA(mfa_plot_uv28_cb, choix = "group")
ggsave("mfa_groups_uv28_cb.jpg",
       plot = mfa_groups_uv280_cb,
       width = 30,
       height = 15,
       units = 'cm',
       dpi = 300)
browseURL("mfa_groups_uv28_cb.jpg")

# Fluorescence cb
mfa_plot_fl_cb <- MFA(
  merged_fl_cb_spectra[-1],
  group = c(26, 25, 25, 25, 25, 25),
  type = c(rep("s", 6)),
  ncp = 4,
  name.group = c("AVN", "CDB", "DTK", "FRV", "KZN", "PDB"),
  graph = TRUE
)
mfa_groups_fl_cb <- plot.MFA(mfa_plot_fl_cb, choix = "group")
ggsave("mfa_groups_fl_cb.jpg",
       plot = mfa_groups_fl_cb,
       width = 30,
       height = 15,
       units = 'cm',
       dpi = 300)
browseURL("mfa_groups_fl_cb.jpg")

# Fluorescence cb
mfa_plot_fl_sb <- MFA(
  merged_fl_sb_spectra[-1],
  group = c(26, 25, 25, 25, 25, 25),
  type = c(rep("s", 6)),
  ncp = 4,
  name.group = c("AVN", "CDB", "DTK", "FRV", "KZN", "PDB"),
  graph = TRUE
)
mfa_groups_fl_sb <- plot.MFA(mfa_plot_fl_sb, choix = "group", graph.type = "ggplot")
ggsave("mfa_groups_fl_sb.jpg",
       plot = mfa_groups_fl_sb,
       width = 30,
       height = 15,
       units = 'cm',
       dpi = 300)
browseURL("mfa_groups_fl_sb.jpg")

# Visualize the spectra
# uv 280 sb
merged_uv280_sb_spectra <- merged_uv280_sb_spectra %>%
  pivot_longer(
    !rt,
    names_to = "samples",
    values_to = "peak_area",
    values_drop_na = TRUE
  )
spectra_280_sb <- ggplot(
  merged_uv280_sb_spectra,
  aes(x = rt, y = peak_area, colour = samples)
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
  ))

ggsave("uv_280_sb_new_spectral_overlay.jpg",
       plot = spectra_280_sb,
       width = 30,
       height = 12,
       units = "cm")
browseURL("uv_280_sb_new_spectral_overlay.jpg")

# uv 280 cb
merged_uv280_cb_spectra <- merged_uv280_cb_spectra %>%
  pivot_longer(
    !rt,
    names_to = "samples",
    values_to = "peak_area",
    values_drop_na = TRUE
  )
spectra_280_cb <- ggplot(
  merged_uv280_cb_spectra,
  aes(x = rt, y = peak_area, colour = samples)
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
  ))

ggsave("uv_280_cb_new_spectral_overlay.jpg",
       plot = spectra_280_cb,
       width = 30,
       height = 12,
       units = "cm")
browseURL("uv_280_cb_new_spectral_overlay.jpg")

# fluorescence sb
merged_fl_sb_spectra <- merged_fl_sb_spectra %>%
  pivot_longer(
    !rt,
    names_to = "samples",
    values_to = "peak_area",
    values_drop_na = TRUE
  )
spectra_fl_sb <- ggplot(
  merged_fl_sb_spectra,
  aes(x = rt, y = peak_area, colour = samples)
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
  ))

ggsave("fl_sb_new_spectral_overlay.jpg",
       plot = spectra_fl_sb,
       width = 30,
       height = 12,
       units = "cm")
browseURL("fl_sb_new_spectral_overlay.jpg")

# fluorescence cb
merged_fl_cb_spectra <- merged_fl_cb_spectra %>%
  pivot_longer(
    !rt,
    names_to = "samples",
    values_to = "peak_area",
    values_drop_na = TRUE
  )
spectra_fl_cb <- ggplot(
  merged_fl_cb_spectra,
  aes(x = rt, y = peak_area, colour = samples)
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
  ))

ggsave("fl_cb_new_spectral_overlay.jpg",
       plot = spectra_fl_cb,
       width = 30,
       height = 12,
       units = "cm")
browseURL("fl_cb_new_spectral_overlay.jpg")
