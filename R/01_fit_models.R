
# ============================================================================== -
# Process all RGB, Spectral, and scoring data
# Author: Jonas Anderegg; jonas.anderegg@usys.ethz.ch
# Last edited: 2023-11-27
# ============================================================================== -

# ============================================================================== -
# Prepare work space ----
# ============================================================================== -

rm(list = ls())


list.of.packages <- c("tidyverse", "mvoutlier", "data.table", 
                      "ggsci", "ggpubr", "doParallel", "foreach", "nls.multstart",
                      "scam")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE, repos = 'https://stat.ethz.ch/CRAN/')

library(tidyverse)
library(data.table)
library(foreach)
library(doParallel)
library(nls.multstart)
library(scam)
library(progress)
library(R.utils)
library(ggpubr)
library(ggsci)

git_dir <- "C:/dynamic_senescence_modeling"

# set working directory
path_to_data <- file.path(git_dir,"data")
path_to_funcs <- file.path(git_dir,"R")

# load functions
source(file.path(path_to_funcs, "A_Utils.R"))
source(file.path(path_to_funcs, "B_Helpers.R"))

setwd(path_to_data)

# ============================================================================== -

# list all files to process ----

# RGB
files_rgb_example <- list.files("RGB_example", pattern = ".csv", full.names = T)

# SPC
files_spc_example <- list.files("SPC_example", pattern = "_LM.csv", full.names = T)

# gather all
files <- c(files_rgb_example, files_spc_example)


# ============================================================================== -
# Run and save ----
# ============================================================================== -

# Writing to network drives not working - write to scratch
out_dir <- file.path(git_dir,"Temp")
out_pars <- file.path(git_dir,"Temp","pars")
out_fits <- file.path(git_dir,"Temp","fits")

# set up cluster and iterate over the output files
numCors <- round(detectCores()*0.8)
cl <- parallel::makeCluster(numCors)
doParallel::registerDoParallel(cl)
foreach(i = 1:length(files),
        .packages = c("tidyverse", "nls.multstart", "scam",  "progress", "R.utils"),
        .export = c("Gompertz_constrained", "Gompertz_flex"),
        .verbose = TRUE) %dopar% {

          # extract ids
          bn <- basename(files[i])
          print(bn)
          bn <- gsub("HSV_", "HSV-", bn)
          exp <- extract_experiment_id(name = bn)
          id <- extract_full_id(name = bn)
          
          # Read data
          dat <- data.table::fread(files[i]) %>%
            dplyr::select(plot.UID, date, index_name, median_value) %>%
            mutate(index_name = paste0("SI_", index_name)) %>% 
            mutate(date = as.Date(date, format="%Y-%m-%d"))

          # Filter relevant time span
          dat <- filter_dates(dat, id)
          
          # reshape for processing
          d_wide <- dat %>% pivot_wider(names_from = index_name, values_from = median_value) %>%
            rename("plot_UID" = plot.UID) %>%
            mutate(date = as.Date(date)) %>%
            group_by(plot_UID, date) %>% nest() %>%
            rename("SVI" = data,
                   "meas_date" = date)

          # get to required format
          exp_plots <- unique(d_wide$plot_UID)
          SVI <- d_wide %>%
            dplyr::filter(plot_UID %in% exp_plots) %>%
            dplyr::select(plot_UID, meas_date, SVI) %>%
            mutate(SVI = purrr::map(SVI, data.table::as.data.table)) %>%
            mutate(Treatment = "None")
          
          # indices to process
          svi <- gsub("SI_", "", names(SVI$SVI[[1]]))
          
          # at least 5 measurements must be available for each plot
          # remove otherwise
          SVI <- check_series(data = SVI, min_length = 5)
          
          # scale SVI data
          data <- scale_SVI(data = SVI,
                            plotid = "plot_UID",
                            treatid = "Treatment",
                            plot = F)

          # # subset data - ONLY FOR TESTING
          # plots <- unique(data)$Plot_ID[1:2]
          # data <- data %>% filter(Plot_ID %in% plots)

          # run and save
          parameters <- get_svi_dynamics(data = data,
                                         svi = svi,
                                         type = "scaled",
                                         method = c("linear", "cgom", "fgom", "pspl"),
                                         timevar = "dafm",
                                         pred_freq = 0.25,
                                         plot_dynamics = TRUE,
                                         output_dir = out_dir)

          # save only parameters
          saveRDS(parameters, paste0(out_fits, "/", id, ".rds"))

          # save only parameters
          pars <- parameters$parameters
          saveRDS(pars, paste0(out_pars, "/", id, ".rds"))
          
        }

parallel::stopCluster(cl)
registerDoSEQ()

# ============================================================================== -


