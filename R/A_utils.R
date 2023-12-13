
# ============================================================================== -
# Functions for the processing of time-series of Spectral, rgb, and scoring data.
# Author: Jonas Anderegg; jonas.anderegg@usys.ethz.ch
# Last edited: 2023-11-27
# ============================================================================== -

# main functions ----

#' Scales spectral indices to range from 0 to 10 
#' representing the minimum and maximum value observed in a time series of measurements
#' @param data A tibble with SVI values in a list column of data.tables
#' @param plotid The variable name of the plot identifier
#' @return A tibble with list columns
#' @details This function reverts the scale for SVI with an increase during the measurement period
scale_SVI <- function(data, plotid, treatid, invert_only = F, 
                      plot = T, topdf = F) {
  
  # fix identifiers
  colid <- which(names(data)==plotid)
  names(data)[colid] <- "Plot_ID"
  colid <- which(names(data)==treatid)
  names(data)[colid] <- "Treatment"
  data$meas_date <- as.Date(data$meas_date, "%Y%m%d")
  
  # keep only plots with measurements covering the entire measurement period
  # get start and end date
  min_date <- min(data$meas_date)
  max_date <- max(data$meas_date)
  # transform measurement date in days after first measurement
  # this should be replaced by a measure of chronological or thermal time after heading
  data$dafm <- as.numeric(data$meas_date - min_date)
  
  # extract SVI from list column and add required metadata
  meta <- data[c("Plot_ID", "meas_date", "dafm", "Treatment")]
  SVI_dat <- data.table::rbindlist(data[["SVI"]])
  d <- cbind(meta, SVI_dat)
  
  # filter dataset
  d_ <- d %>% group_by(Plot_ID) %>% nest() %>%
    mutate(max = purrr::map_chr(data,  ~as.character(max(.$meas_date))),
           min = purrr::map_chr(data,  ~as.character(min(.$meas_date)))) %>%
    filter(max == max_date & min == min_date)
  # remove helper columns
  d_ <- d_[!(names(d_) %in% c("max", "min"))]
  
  # scale SVI
  if(!invert_only){
    d_mutated <- d_ %>% 
      transmute(SVI_sc = purrr::map(data, col_scaling)) %>% 
      unnest(cols = c(SVI_sc))
  } else {
    d_mutated <- d_ %>% 
      unnest(cols = c(data))
  }
  
  # revert scale where required
  ids <- d_mutated[c("Plot_ID", "meas_date", "dafm", "Treatment")] %>% ungroup()
  SVI_sc_dat <- d_mutated[grepl("^SI_", names(d_mutated))]
  name_maintainer <- names(SVI_sc_dat)
  
  SVI_sc_dat <- SVI_sc_dat %>% 
    mutate_all(list(r = revert)) %>% 
    #select original or reversed values
    dplyr::select_if(function(col) col[1] > 5) %>% 
    data.table::as.data.table()
  
  if(names(SVI_sc_dat) == "r"){
    names(SVI_sc_dat) <- paste0(name_maintainer, "_r")
  }
  
  if(plot){
    
    # transform measurement date in days after first measurement
    # this should be replaced by a measure of chronological or thermal time after heading
    pd <- cbind(ids, SVI_sc_dat)
    pd$dafm <- as.numeric(as.Date(pd$meas_date, "%Y%m%d") - as.Date(min_date, "%Y%m%d"))
    # reshape data for plotting
    pd <- pd %>% dplyr::select(1:Treatment, dafm, starts_with("SI_")) %>% 
      tidyr::gather(SVI, value, starts_with("SI_"))
  
    p <- ggplot(pd) +
      geom_line(aes(x = dafm, y = value, group = Plot_ID, col = Treatment), alpha  = 0.3) +
      facet_wrap(~SVI) +
      ggsci::scale_color_npg() +
      theme_bw(base_size = 7)
    
    # check if Output directory exists
    if(!file.exists(paste0(path_to_data, "Output"))){
      dir.create(paste0(path_to_data, "Output"))
    }
    
    if(topdf){
      # save plot to pdf
      pdf(paste0(path_to_data, "Output/SVI_ts.pdf"), width = 35, height = 35)
      plot(p)
      dev.off()
    } else {
      plot(p)
    }

  }
  
  # create data output (list column)
  spc_pre_list <- map(purrr::transpose(SVI_sc_dat), data.table::as.data.table)
  if(!invert_only){
    ids[, "SVI_sc"] <- list(spc_pre_list)
    # join with input data
    data_out <- full_join(data, ids, by = c("Plot_ID", "meas_date", "dafm", "Treatment"))
  } else {
    ids[, "SVI"] <- list(spc_pre_list)
    # join with input data - remove old SVI data
    data <- data %>% dplyr::select(-SVI)
    data_out <- full_join(data, ids, by = c("Plot_ID", "meas_date", "dafm", "Treatment"))
  }
  
} 


#' Gets dynamics parameters for each index and plot
#' @param data A tibble with scaled SVI values in a list column of data.tables, as returned by "scale_SVI"
#' @param svi A vector of index names to process
#' @param type "scaled" or "raw". Required for meaningful plotting
#' @param timevar A character string specifying the name of the column containing the time variable
#' @param method A character string specifying the method to be used for modelling of the temporal index dynamics. 
#' Any combination of "linear", "cgom", "fgom", "pspl"
#' @param plot_dynamics Boolean, indicating whether or not to create plots of the fits
#' @param output_dir The output directory
#' @return A tibble containing the dynamics parameters for each Plot and SVI. 
get_svi_dynamics <- function(data, 
                             svi, 
                             type,
                             timevar, 
                             pred_freq,
                             method, 
                             plot_dynamics,
                             output_dir){
  
  # function to fit spline
  p_spline <- function(data, x_name, y_name) {
    pb$tick()
    y <- data[[y_name]]
    x <- data[[x_name]]
    tryCatch({
      spl <- withTimeout({
        scam(y ~ s(x, k = round(length(x) * 0.8), bs = "mpd"),
             optimizer = "bfgs")
      }, timeout = 2)
    }, TimeoutException = function(ex) {
      message("Timeout. Skipping.")
    })
    return(spl)
  }
  
  # function to fit constrained Gompertz
  nls_mult_const <- function(data){
    pb$tick()
    fit <- nls_multstart(value ~ Gompertz_constrained(b, M, tempsum = timevar),
                         data = data, 
                         iter = 300, 
                         start_lower = c(b = -0.2, M = 15),
                         start_upper = c(b = 0.1, M = 25),
                         convergence_count = 150, 
                         supp_errors = "Y")
    return(fit)
    }
   
  # function to fit flexible Gompertz
   nls_mult_flex <- function(data){
     pb$tick()
     fit <- nls_multstart(value ~ Gompertz_flex(A, C, b, M, tempsum = timevar),
                          data = data, 
                          iter = 500, 
                          start_lower = c(A = -1000, C = 5, b = -0.2, M = 15),
                          start_upper = c(A = 1000, C = 10000, b = 0.1, M = 25),
                          convergence_count = 150, 
                          supp_errors = "Y")
     return(fit)
   }
  
  print(paste("processing", length(svi), "indices ..."))
  
  # ============================================================================ -
  # prepare data 
  # ============================================================================ -
  
  # for variable selection
  pattern <- paste0(svi, collapse = "|")
  
  # un-nest SVI data
  if(type == "scaled"){
    dat_svi <- data.table::rbindlist(data[["SVI_sc"]])
  } else {
    dat_svi <- data.table::rbindlist(data[["SVI"]])
  }
  namevec <- grep(pattern = pattern, names(dat_svi), value = T)
  dat_svi <- dat_svi[, ..namevec]
  
  # fix variable names
  names(data)[which(names(data)==timevar)] <- "timevar"
  
  # reshape
  meta <- data[, names(data) %in% c("Plot_ID", "Treatment", "timevar")]
  dat <- cbind(meta, dat_svi)
  dat_long <- reshape2::melt(
    dat, 
    id.vars = c("Plot_ID", "Treatment", "timevar"),
    measure.vars = c(grep("^SI_", names(dat), value = T)))
    # measure.vars = c(grep("_cover", names(dat), value = T)))
  
  # ============================================================================ -
  # fit models
  # ============================================================================ -

  # group data
  fits <- dat_long %>%
    dplyr::group_by(Plot_ID, variable, Treatment) %>%
    tidyr::nest()

  # perform linear interpolations
  if("linear" %in% method){

    print(" -- Interpolating ...")

    # get linear interpolations
    fits <- fits %>%
      mutate(preds_lin = purrr::map(.x = data, pred_freq = pred_freq,
                                    .f = lin_approx))
  }

  # df for predictions
  new_preds <- dat %>% ungroup() %>%
    do(., data.frame(
      timevar = seq(min(.$timevar), max(.$timevar), by = pred_freq),
      stringsAsFactors = FALSE)
    )

  # Fit Gompertz models
  if ("cgom" %in% method){

    print(" -- Fitting constrained Gompertz ...")

    num_ticks <- n_groups(fits)
    pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
                           total = num_ticks)

    # fit Gompertz models
    fits <- fits %>%
      mutate(fit_cgom = purrr::map(.x = data,
                                   .f = nls_mult_const))

    # max and min for each curve
    max_min <- group_by(dat, Plot_ID) %>%
      summarise(., min_gGDDAH = min(timevar), max_gGDDAH = max(timevar)) %>%
      ungroup()
    # get predictions
    fits <- fits %>%
      mutate(preds_cgom = purrr::map(fit_cgom, broom::augment, newdata = new_preds)) %>%
      mutate(preds_cgom = purrr::map(preds_cgom, tidy_gompertz_preds, method = "cgom"))
  }

  # Fit Gompertz models
  if ("fgom" %in% method){

    print(" -- Fitting flexible Gompertz ...")

    num_ticks <- n_groups(fits)
    pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
                           total = num_ticks)

    # fit Gompertz models
    fits <- fits %>%
      mutate(fit_fgom = purrr::map(.x = data,
                                   .f = nls_mult_flex))

    # max and min for each curve
    max_min <- group_by(dat, Plot_ID) %>%
      summarise(., min_gGDDAH = min(timevar), max_gGDDAH = max(timevar)) %>%
      ungroup()
    # get predictions
    fits <- fits %>%
      mutate(preds_fgom = purrr::map(fit_fgom, broom::augment, newdata = new_preds)) %>%
      mutate(preds_fgom = purrr::map(preds_fgom, tidy_gompertz_preds, method = "fgom"))
  }

  # Fit p-splines
  if ("pspl" %in% method){

    print(" -- Fitting p-splines ...")

    num_ticks <- n_groups(fits)
    pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) elapsed :elapsed eta :eta",
                           total = num_ticks)

    # fit splines
    fits <-  fits %>%
      mutate(fit_pspl = purrr::map(.x = data, x_name="timevar", y_name="value",
                                   # .f = p_spline))
                                   .f = possibly(p_spline, otherwise = NA_real_)))

    # make predictions
    fits <- fits %>%
      mutate(preds_pspl = purrr::map(.x = fit_pspl, newdata = new_preds,
                                     .f = possibly(predict_p_spline, otherwise = NA_real_)))

  }

  if(!any(method %in% c("linear", "pspl", "cgom", "fgom"))) {
    stop(paste("interpolation method not supported"))
  }

  # ============================================================================ -
  # Plot fits
  # ============================================================================ -

  if (plot_dynamics){

    print(" -- Creating plots ...")

    # remove plots with missing predictions
    if("preds_pspl" %in% names(fits)){
      fits <- fits %>%
        mutate(preds_ = purrr::map_lgl(preds_pspl, is.data.frame)) %>%
        filter(preds_ == TRUE) %>%
        dplyr::select(-preds_)
    }

    # extract data
    preds <- fits[c("Plot_ID", "Treatment", "variable",
                  grep("preds_", colnames(fits), value = T))] %>%
      unnest(cols = starts_with("preds_"))
    preds <- cbind(new_preds, preds) %>% as_tibble %>%
      tidyr::pivot_longer(., cols = 6:ncol(.),
                          names_to = "method",
                          values_to = "prediction") %>%
      mutate(method = gsub("preds_", "", method))

    if (type == "any"){
      axis_limits = c(-0.1, 1.1)
    } else {
      axis_limits = c(-.5, 10.5)}

    # basic plot layout
    p0 <- ggplot() +
      scale_colour_manual(values = c("green4", "yellow", "black", "red")) +
      scale_fill_manual(values = c("green4", "yellow", "black", "red")) +
      scale_y_continuous(breaks = seq(0,10,2), limits = axis_limits) +
      geom_abline(slope = 0, intercept = 8.5) +
      geom_abline(slope = 0, intercept = 5.0) +
      geom_abline(slope = 0, intercept = 1.5) +
      theme_bw( base_family = "Helvetica") +
      theme(panel.grid = element_blank()) +
      ylab("Index value") +
      xlab(timevar)

    # treatment-wise plots ====================================================== -

    # empty lists for plots
    mets_list = list()

    # loop over interpolation methods
    for (met in unique(preds$method)){

      plot_id <- paste(met)

      pd <- preds[preds$method == met,]
      # pd <- pd[pd$Treatment %in% c("F0I", "FI"),]

      mets_list[[met]] <- p0 +
        geom_line(aes(timevar, prediction, group = interaction(Treatment, Plot_ID), colour = Treatment),
                  alpha = 0.75, pd) +
        # scale_colour_manual(values = c("green4", "black")) +
        facet_wrap(~variable) +
        ggtitle(plot_id)
    }

    # plot-wise plots including raw data points ================================ -

    subset_plots <- sample(unique(fits$Plot_ID), min(100, length(unique(fits$Plot_ID))))

    for (var in unique(fits$variable)){

      # extract data

      pd <- fits[fits$variable == var, ]
      pd <- pd[pd$Plot_ID %in% subset_plots, ]
      obs <- pd[c("Plot_ID", "Treatment", "variable", "data")] %>%
        unnest(c(data))
      preds <- pd[c("Plot_ID", "Treatment", "variable",
                    grep("preds_", colnames(pd), value = T))] %>%
        unnest(cols = grep("^preds_", names(.), value = T))
      preds <- cbind(new_preds, preds) %>% as_tibble %>%
        tidyr::pivot_longer(., cols = 5:ncol(.),
                            names_to = "method",
                            values_to = "prediction")
      # create plot
      p <- p0 +
        geom_point(aes(timevar, value), shape = 1, obs) +
        geom_line(aes(timevar, prediction, group = method, colour = method),
                  alpha = 0.75, preds) +
        facet_wrap(~Plot_ID, labeller = labeller(.multi_line = FALSE)) +
        ggtitle(var)

      # save to pdf
      pdf(paste0(output_dir, "/", id, ".pdf"), width = 35, height = 35)
      plot(p)
      dev.off()
      }
  }

  # ============================================================================ -
  # Extract model parameters
  # ============================================================================ -

  print(" -- Extracting parameters ...")

  if ("linear" %in% method){
    mod_pars <- fits %>%
      mutate(pars_lin = purrr::map(.x = preds_lin,
                                   .f = extract_pars,
                                   new_data = new_preds) %>%
               purrr::map(cbind.data.frame))
  }
  if ("cgom" %in% method){
    # mod_pars <- mod_pars
    mod_pars <- mod_pars %>%
      mutate(pars_cgom = purrr::map(fit_cgom, broom::tidy)) %>%
      mutate(pars_cgom = purrr::map(pars_cgom, tidy_gompertz_output)) %>%
      mutate(pars2_cgom = purrr::map2(.x = preds_cgom,
                                      .y = fit_cgom,
                                      .f = extract_pars,
                                      new_data = new_preds) %>%
               purrr::map(cbind.data.frame))
  }
  if ("fgom" %in% method){
    mod_pars <- mod_pars %>%
      mutate(pars_fgom = purrr::map(fit_fgom, broom::tidy)) %>%
      mutate(pars_fgom = purrr::map(pars_fgom, tidy_gompertz_output)) %>%
      mutate(pars2_fgom = purrr::map2(.x = preds_fgom,
                                      .y = fit_fgom,
                                      .f = extract_pars,
                                      new_data = new_preds) %>%
               purrr::map(cbind.data.frame))
  }
  if ("pspl" %in% method){
    mod_pars <- mod_pars %>%
      mutate(pars_pspl = purrr::map(.x = preds_pspl, new_data = new_preds,
                                    .f = extract_pars) %>%
               purrr::map(cbind.data.frame))
  }

  # ============================================================================ -
  # create a tidy output tibble holding parameter values
  # ============================================================================ -

  parameters <- mod_pars[c("Plot_ID", "Treatment", "variable",
                           grep("pars", colnames(mod_pars), value = T))] %>%
    unnest(grep("pars_|pars2_", names(.)), names_sep = "_")

  names(parameters) <- gsub("pars_|pars2_|pars_timepoints", "", names(parameters))

  # ============================================================================ -

  out <- list("fits" = fits, "parameters" = parameters, "timepoints" = new_preds)
  
  return(out)
  # return(fits)

}

# ============================================================================================================= -

# helper functions ---- 

revert <- function(x){10 - x}

col_scaling <- function(d) {
  ids <- d[!grepl("^SI_", names(d))]
  d <- d[grep("^SI_", names(d))]
  d <- purrr::map_df(d, function(X) (X - min(X, na.rm = TRUE))/diff(range(X, na.rm = TRUE))*10)
  out <- bind_cols(ids, d)
  return(out)
}

lin_approx <- function(data, pred_freq){
  data <- as.data.frame(data)
  # linearly interpolate between measurement time points
  preds <- approx(data[, "timevar"], data[,"value"],
                xout = seq(min(data$timevar), max(data$timevar), 
                           by = pred_freq))
  preds <- tibble(preds_lin = preds$y)
  return(preds)
}

# Constrained Gompertz equation
Gompertz_constrained <- function(b, M, tempsum) {
  grenness_decay <- 10*exp(-exp(-b*(tempsum-M)))
  return(grenness_decay)
}

# Flexible Gompertz equation
Gompertz_flex <- function(A, C, b, M, tempsum) {
  grenness_decay <- A + C*exp(-exp(-b*(tempsum-M)))
  return(grenness_decay)
}

p_spline <- function(data, x_name, y_name) {
  pb$tick()
  y <- data[[y_name]]
  x <- data[[x_name]]
  spl <- scam(y ~ s(x, k = round(length(x) * 0.8), bs = "mpd"),
              optimizer = "bfgs")
  return(spl)
  
}

predict_p_spline <- function(spl, newdata = new_preds) {
  names(newdata) <- "x"
  preds <- predict.scam(spl, newdata = newdata)
  preds <- tibble(preds_pspl = unname(preds))
  return(preds)
}

extract_pars <- function(data, fit, new_data){
  
  # check if model could be fitted
  if(is_tibble(data) && !any(is.na(data[[1]]))){

    # to avoid "extremely bad integrand behavior" error
    data[data < 0] <- 0
    # find method predictions column name
    colid <- grep("^preds_", names(data), value = T)
    # add the time frame of predictions
    data <- bind_cols(data, new_data)
    
    # extract parameters
    t80 <- data[which(data[[colid]] < 8)[1], "timevar"]$timevar
    t50 <- data[which(data[[colid]] < 5)[1], "timevar"]$timevar
    t20 <- data[which(data[[colid]] < 2)[1], "timevar"]$timevar
    dur1 <- t50 - t80
    dur2 <- t20 - t80
    
    # extract interval under curve/interpolation
    # helper functions
    f1 <- approxfun(data$timevar, data[[colid]])   
    f2 <- function(x) abs(f1(x))
    # extract the difference between the curves
    diff <- approx(data$timevar, data[[colid]], n = 100)$y
    # get the integral of the difference curve
    Integral <- integrate(f2, min(data$timevar), max(data$timevar), 
                          subdivisions = 2000)$value
    
    # get the minimum and maximum of the second derivatives as onset and end
    # ONLY possible for Gompertz models
    if(grepl("gom", colid)){
      if(grepl("cgom", colid)){
        expr <- expression(10*exp(-exp(-b*(x-M))))
        coeff <- coefficients(fit)
        fun <- do.call(substitute, list(as.list(expr)[[1]], env=list(b=as.numeric(coeff["b"]), M=as.numeric(coeff["M"]))))
      }
      if(grepl("fgom", colid)){
        expr <- expression(A + C*exp(-exp(-b*(x-M))))
        coeff <- coefficients(fit)
        fun <- do.call(substitute, list(as.list(expr)[[1]], env=list(A=as.numeric(coeff["A"]),C=as.numeric(coeff["C"]),b=as.numeric(coeff["b"]), M=as.numeric(coeff["M"]))))
      }
      
      asymptote1deriv <- DD(fun,"x", 1)
      data_1deriv <- eval(asymptote1deriv,list(x=c(data$timevar)))
      
      asymptote2deriv <- DD(fun,"x",2)
      data_2deriv <- eval(asymptote2deriv,list(x=c(data$timevar)))
      
      onsen_idx <- which.min(data_2deriv)
      onsen <- data$timevar[onsen_idx]
      endsen_idx <- which.max(data_2deriv)
      endsen <- data$timevar[endsen_idx]
      
    } else {
      onsen <- endsen <- NA
    }
    
  # if no model object exists, make parameters NA
  } else {
    onsen <- endsen <- t80 <- t50 <- t20 <- dur1 <- dur2 <- Integral <- NA
  }
  
  pars <- cbind(onsen, endsen, t80, t50, t20, dur1, dur2, Integral)
  
  return(pars)
  
}

DD <- function(expr, name, order = 1) {
  # calculate higher order derivative
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}


tidy_gompertz_preds <- function(data, method){
  out <- tibble(preds = data$.fitted)
  names(out)[1] <- paste0("preds_", method)
  return(out)
}

tidy_gompertz_output <- function(data){
  out <- data[,c("term", "estimate")] %>% 
    tidyr::spread(term, estimate)
  return(out)
}

