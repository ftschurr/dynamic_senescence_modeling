
# ============================================================================== -
# Helper function required for the joint processing of time-series of 
# spectral, rgb, and scoring data.
# Author: Jonas Anderegg; jonas.anderegg@usys.ethz.ch
# Last edited: 2023-11-27
# ============================================================================== -

extract_experiment_id <- function(name){
  parts <- unlist(strsplit(name, "_"))
  exp_id <- parts[grepl("FPWW|ESWW", parts)]
}

extract_full_id <- function(name){
  
  if(grepl("^ETHZ", name)){
    id <- paste(strsplit(bn, "_")[[1]][3:5], collapse ="_")
  } else if (grepl("^FPWW", name)){
    id <- paste(strsplit(bn, "_")[[1]][1:3], collapse ="_")
  } else if (grepl("^ESWW", name)){
    if(grepl("lot", name)){
      id <- paste(strsplit(bn, "_")[[1]][1:3], collapse ="_")
    } else {
      id <- paste(strsplit(bn, "_")[[1]][1:2], collapse ="_")
    }
  }
}

filter_dates <- function(data, id){
  if(grepl("FPWW022", id)){
    data <- data %>% 
      dplyr:: filter(date >= as.Date("2018-05-28") & date <= as.Date("2018-07-13"))
    } else if(grepl("FPWW032", id)){
      data <- data %>% 
        dplyr::filter(date >= as.Date("2021-05-27") & date <= as.Date("2021-07-19"))
    } else if(grepl("FPWW033", id)){
      data <- data %>% 
        dplyr::filter(date >= as.Date("2022-05-27") & date <= as.Date("2022-07-18"))
    } else if(grepl("ESWW006", id)){
      data <- data %>% 
        dplyr::filter(date >= as.Date("2022-05-27") & date <= as.Date("2022-07-14"))
    }
}

check_series <- function(data, min_length){
  
  # identify incomplete series
  DD <- data %>% group_by(plot_UID) %>% nest() %>% 
    mutate(n = purrr::map_dbl(data, nrow))
  incomplete_plot_series <- DD %>% filter(n < 5) %>% pull(plot_UID)
  
  # remove from data set
  out <- data %>% 
    filter(!plot_UID %in% incomplete_plot_series)
  
  return(out)
  
}
