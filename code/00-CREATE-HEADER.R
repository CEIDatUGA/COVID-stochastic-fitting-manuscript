# This script creates header information for all states that the
# cluster array script can read in for fitting the model.



# Clean workspace ---------------------------------------------------------

rm(list = ls(all.names = TRUE))



# Set global parameters ---------------------------------------------------

TEST <- TRUE
DATA_SOURCE <- "JHU"
END_DATE <- as.Date("2020-12-31")
STATES <- ifelse(TEST, "Georgia", c(state.name, "District of Columbia"))



# Load libraries ----------------------------------------------------------

library('lubridate') #needs to loaded first so it doesn't mess with 'here'
library('readr')
library('doParallel')
library('foreach')
library('purrr')
library('ggplot2')
library('pomp')
library('dplyr')
library('tidyr')
library('here')
library('vctrs')
library('dplyr')



# Source functions --------------------------------------------------------

# setting all parameters, specifying those that are  fitted
source(here::here("code/model-setup/setparsvars.R")) 

# loads data from JHU
source(here::here("code/data-processing/loadcleandata.R"))

# function that processes and retrieves mobility covariate
source(here::here("code/data-processing/loadcleanucmobility.R")) 




# Generate time stamp -----------------------------------------------------

tm <- .POSIXct(Sys.time(), "US/Eastern")  # time stamp with date, hours, minutes
timestamp <- paste(lubridate::date(tm),
               stringr::str_pad(as.character(lubridate::hour(tm)), 
                                width = 2, side = "left", pad = "0"),
               stringr::str_pad(as.character(lubridate::minute(tm)), 
                                width = 2, side = "left", pad = "0"),
               sep='-')




# Specify initial RO conditions for fitting -------------------------------

statedf <- readRDS(here::here("data/us_popsize.rds")) %>% 
  dplyr::filter(state_full %in% STATES) %>% 
  dplyr::mutate(init = "fresh") %>% 
  # R0 at beginning of epidemic for each state, choosing max possible
  # since mobility and latent trend can only reduce
  dplyr::mutate(initR0 = 7) 



# Get data ----------------------------------------------------------------

all_states_pomp_data <- loadcleandata(
  datasource = DATA_SOURCE,
  locations = STATES,
  vars = c("cases", "hosps", "deaths"),
  smooth = FALSE,  # use raw data no smoothing
  # trim leading zeros (trim to first reported case or death for each state)
  trim = TRUE
) %>%  
  dplyr::filter(date <= END_DATE)

# add in initial NA data at t = 0 for all states to 
# make initial estimation easier (pomp trick)
first_date <- all_states_pomp_data %>%
  group_by(location) %>% 
  filter(date == min(date)) %>%
  ungroup() %>%
  mutate(date = date - 1,
         time = 0,
         cases = NA, hosps = NA, deaths = NA)

all_states_pomp_data <- bind_rows(first_date, all_states_pomp_data)



# Get the mobility metrics ------------------------------------------------

all_states_pomp_covar <- loadcleanucmobility(
  location = STATES,
  pomp_data = all_states_pomp_data
) %>%
  dplyr::filter(date <= END_DATE)




# Set up parameter bounds for latin hypercube sample ----------------------

est_these_pars <- c(
  "log_sigma_dw",
  "min_frac_dead",
  "max_frac_dead",
  "log_half_dead",
  "log_theta_cases",
  "log_theta_deaths"
)

est_these_inivals <- c(
  "E1_0", 
  "Ia1_0",
  "Isu1_0", 
  "Isd1_0"
)

par_var_bounds <- list(
  lowers = c(log_sigma_dw = -5,  # log scale
             min_frac_dead = -6,  # logit scale
             max_frac_dead = -6,  # logit scale
             log_half_dead = -5,  # log scale
             log_theta_cases = -5,  # log scale
             log_theta_deaths = -5,  # log scale
             E1_0 = 0,  # log scale
             Ia1_0 = 0,  # log scale
             Isu1_0 = 0,  # log scale
             Isd1_0 = 0),  # log scale
  
  uppers = c(log_sigma_dw = 5,  # log scale
             min_frac_dead = 6,  # logit scale
             max_frac_dead = 6,  # logit scale
             log_half_dead = 5,  # log scale
             log_theta_cases = 5,  # log scale
             log_theta_deaths = 5,  # log scale
             E1_0 = 10,  # log scale
             Ia1_0 = 10,  # log scale
             Isu1_0 = 10,  # log scale
             Isd1_0 = 10)  # log scale
)



# Loop over states and store model fitting information --------------------

# initialize large list that will hold pomp model and other info for each state
pomp_list <- vector("list", length(STATES))

for(i in 1:length(STATES)) {
  dolocation <- STATES[i]
  print(sprintf('starting state %s', dolocation))
  
  # This will be appended to each saved file 
  filename_label <- dolocation
  
  # Get pomp data for location
  pomp_data <- all_states_pomp_data %>%
    filter(location == dolocation)
  
  # calculate number of knots for location
  n_knots <- round(nrow(pomp_data) / 21)
  knot_coefs <-  paste0("b", 1:n_knots)
  
  # Set the parameter values and initial conditions
  par_var_list <- setparsvars(
    est_these_pars = c(est_these_pars, knot_coefs),
    est_these_inivals = est_these_inivals,
    population = statedf %>%
      filter(state_full == dolocation) %>% pull(total_pop),
    n_knots = n_knots,
    # set R0 at beginning of epidemic
    rnaught = statedf %>%
      filter(state_full == dolocation) %>% pull(initR0)
  )  
  
  # add beta inits to par_var_bounds
  beta_names <- paste0("b", 1:n_knots)
  beta_lowers <- rep(-10, n_knots)
  names(beta_lowers) <- beta_names
  beta_uppers <- rep(10, n_knots)
  names(beta_uppers) <- beta_names
  par_var_bounds_all <- list()
  par_var_bounds_all$lowers <- c(par_var_bounds$lowers, beta_lowers)
  par_var_bounds_all$uppers <- c(par_var_bounds$uppers, beta_uppers)
  par_var_list$par_var_bounds <- par_var_bounds_all  # for latin hypercube sample
  
  
  # Get covariate 
  tmp_covar <- all_states_pomp_covar %>%
    filter(location == dolocation)
  
  covar <- covariate_table(
    t = pomp_data$time,
    seas = bspline.basis(
      x=t,
      nbasis=n_knots,
      degree=3
    ),
    rel_beta_change = as.matrix(tmp_covar$rel_beta_change),
    trend_sim = as.matrix(rep(10, times = nrow(tmp_covar))),  # this is a placeholder only needed for simulation
    fit = 1,  # 1 = fitting; 0 = simulating
    times="t",
    order = "constant"
  )
  
  # Save all pieces for each state in a list
  # pomp_list[[i]]$pomp_model <- pomp_model 
  pomp_list[[i]]$filename_label <- filename_label
  pomp_list[[i]]$pomp_data <- pomp_data
  pomp_list[[i]]$pomp_covar <- covar
  pomp_list[[i]]$location <- dolocation
  pomp_list[[i]]$par_var_list <- par_var_list
  
} # done serial loop over all states that creates info for fitting



# Save header outputs for array job ---------------------------------------


saveRDS(pomp_list, file = here::here("header/pomp_list.rds"))
saveRDS(timestamp, file = here::here("header/timestamp.rds"))
saveRDS(statedf, file = here::here("header/statedf.rds"))

