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

# Read-in the data
hplc_wines <- readxl::read_excel(path = "/Users/mphomafata/Documents/Work_file/Collaborative Work/HPLC scanning/CSV data/HPLC scouting.xlsx",
                                 sheet = "280 nm SB") # "280 nm CB", "FLD SB", "FLD CB"

# split the file into respective samples

uv_280 <- list()
for (i in 1:(ncol(hplc_wines))) {

first <- as.numeric((2*i)-1)
second <- as.numeric(2*i)
sample_i <- hplc_wines[c(first,second)]
colnames(sample_i)[1] <- c("rt")
sample_i <- list(sample_i)
uv_280 <- append(uv_280, sample_i)
}
name_list <- row.names(as.data.frame(t(hplc_wines)) %>% filter(row_number() %% 2 == 0))
names(uv_280) <- name_list

# create seperate dataframes
# list2env(uv_280,globalenv()) 
# my_data <- c(name_list)

# create a merged data frame block for MFA analysis
merged_uv280_spectra <- full_join(x=as.data.frame(uv_280[[1]]), 
                                  y=as.data.frame(uv_280[[2]]), 
                                  by="rt")
for (i in 3:length(uv_280)){
  merged_uv280_spectra <- full_join(x=as.data.frame(merged_uv280_spectra), 
                                    y=as.data.frame(uv_280[[i]]), 
                                    by="rt"
                                    )
}

# run the MFA
mfa_plot <- MFA(
  merged_uv280_spectra[-1],
  group = c(26,25,25,25,25,25),
  type = c(rep("s", 6)),
  ncp = 4,
  name.group = c("AVN",
                 "CDB",
                 "DTK",
                 "FRV",
                 "KZN",
                 "PDB"),
  graph = TRUE)
mfa_scores <- fviz_mfa_ind(mfa_plot,
             geom = c("point", "text"),
             axes = c(1,2)
)
ggsave(glue("mfa_scores.jpg"),
       plot = mfa_scores,
       width = 30,
       height = 15,
       units = 'cm',
       dpi = 300)
browseURL(glue("mfa_scores.jpg"))


