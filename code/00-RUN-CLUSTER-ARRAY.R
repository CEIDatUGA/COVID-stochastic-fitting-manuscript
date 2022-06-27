# This script is designed to run the entire workflow for estimating
# parameters for the stochastic COVID-19 model and for producing
# estimates of observed states. This script is designed to be run
# on a computing cluster.



# Clean the workspace -----------------------------------------------------

rm(list = ls(all.names = TRUE))



# Set global options ------------------------------------------------------

TEST <- TRUE  # FALSE for production run
LOCAL <- FALSE  # FALSE for cluster run

# set state id if running on local machine
if(LOCAL) a = 1

# run in parallel (TRUE) or not (FALSE)
PARALLEL_FLAG <- ifelse(TEST, TRUE, TRUE)

# number of cores to parallelize across
NUM_CORES <- ifelse(TEST, 5, 32)  

# number of particles for each mif2 iterations
MIF_PARTICLES <- ifelse(TEST, 200, 3500)  

# cooling fractions for each of 4 mif2 rounds
MIF_COOLING <- c(1, 0.9,  # rounds 1 & 2 in batch 1
                 1, 0.5)  # rounds 1 & 2 in batch 2

# number of iterations per each of 4 mif2 rounds
if(TEST) {
  MIF_ITERATIONS <- c(20, 20,  # iterations for rounds 1 & 2 in batch 1
                      10, 10)  # iterations for rounds 1 & 2 in batch 2
} else {
  MIF_ITERATIONS <- c(200, 100,  # iterations for rounds 1 & 2 in batch 1
                      50, 50)  # iterations for rounds 1 & 2 in batch 2
}

# number of independent mif2 runs, run in parallel with different 
# parameter starting conditions
MIF_RUNS <- ifelse(TEST, 10, 100)  

# number of particles for pfilter loglik estimation
PF_LL_PARTICLES <- ifelse(TEST, 10, 5000)

# number of replicate pfilters for loglik estimation
PF_LL_REPS <- ifelse(TEST, 2, 20)

# number of particles per pfitler for smoothed posterior states
PF_STATES_PARTICLES <- ifelse(TEST, 100, 2500)

# number of replicate particle filters for smoothed posterior states
PF_STATES_REPS <- ifelse(TEST, 10, 1000)



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
library('tibble')



# Read in id argument from the bash script --------------------------------
# this indicates which state to run
args <- (commandArgs(TRUE))
## args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  print("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
myargument <- as.numeric(a)  # numeric id for indexing list of pomp objects



# Source functions --------------------------------------------------------

source("../code/model-setup/makepompmodel.R") # function that generates the pomp model
source("../code/model-fitting/runmif_allstates_array.R") # runs mif fitting
source("../code/result-exploration/exploremifresults.R") # explore/summarise mif results



# Load the pomp information for this state --------------------------------

pomp_listr <- readRDS("../header/pomp_list.rds")  # pomp info for all states
this_pomp <- pomp_listr[[myargument]]  # pomp info for this state

# set the number of spline knots to one every 21 days
# this defines the knots for the latent trend spline fit
n_knots <- round(nrow(this_pomp$pomp_data) / 21 ) 

# time stamp for current run
timestamp <- readRDS("../header/timestamp.rds")



# Set parallel options ----------------------------------------------------

parallel_info <- list()  # an empty list to hold all arguments
parallel_info$parallel_run <- PARALLEL_FLAG  # run in parallel
parallel_info$num_cores <- NUM_CORES  # number of cores to parallelize across



# Set arguments for pomp::mif2 --------------------------------------------

mif_settings <- list()  # an empty list for mif2 settings
mif_settings$mif_num_particles  <- MIF_PARTICLES
mif_settings$mif_num_iterations <- MIF_ITERATIONS[1:2]
mif_settings$mif_cooling_fracs <- MIF_COOLING[1:2]
mif_settings$replicates <- MIF_RUNS
mif_settings$pf_num_particles <- PF_LL_PARTICLES
mif_settings$pf_reps <- PF_LL_REPS

  


# Generate the pomp model object ------------------------------------------

pomp_model <- makepompmodel(
  par_var_list = this_pomp$par_var_list,
  pomp_data = this_pomp$pomp_data,
  pomp_covar = this_pomp$pomp_covar,
  n_knots = n_knots
)

# save in the `this_pomp` list
this_pomp$pomp_model <- pomp_model



# Run first batch of 2 mif rounds -----------------------------------------

mif_res1 <- runmif_allstates(
  parallel_info = parallel_info,
  mif_settings = mif_settings,
  pomp_list = this_pomp,
  par_var_list = this_pomp$par_var_list
)

pomp_res1 <- this_pomp  # current state of affairs
pomp_res1$mif_res <- mif_res1  # add the first mif result

# generate traceplots and summarize fits
mif_explore1 <- exploremifresults(
  pomp_res = pomp_res1,
  par_var_list = pomp_res1$par_var_list,
  n_knots = n_knots
) 

pomp_res1$mif_explore <- mif_explore1

# intermediate save of fits
outdir <- "../output/"
state <- this_pomp$location

# batch 1 mif results
saveRDS(pomp_res1, file = paste0(outdir, state, "-mif_res_1.RDS"))



# Generate new starting parameter sets for final mif batch ----------------

# sample parameter sets with loglikelihood weights
nsample <- mif_settings$replicates
sets <- mif_explore1$est_partable %>%
  filter(!is.nan(LogLik)) %>%
  arrange(-LogLik) %>%
  filter(LogLik >= (max(LogLik)-10))  # all sets within 10 loglik
lls <- sets$LogLik
weights <- exp(lls-mean(lls))
weights <- weights / sum(weights)
rids <- sample(1:nrow(sets), size = nsample, replace = TRUE, prob = weights)
newparams <- tibble()
for(i in rids) {
  tmp <- sets[i, ] %>%
    dplyr::select(-MIF_ID, -LogLik, -LogLik_SE)
  for(j in 1:ncol(tmp)) {
    tmp[1,j] <- rnorm(1, mean = tmp[1,j], sd = abs(tmp[1,j]/4))
  }
  # tack on the fixed parameters
  fixed <- this_pomp$par_var_list$allparvals
  fixed <- fixed[!names(fixed) %in% c(this_pomp$par_var_list$params_to_estimate,
                                     this_pomp$par_var_list$inivals_to_estimate)]
  fixed <- enframe(fixed) %>%
    pivot_wider(names_from = "name", values_from = "value")
  tmp <- bind_cols(fixed, tmp)
  newparams <- bind_rows(newparams, tmp)
}

# overwrite the initial par_var_bound
bounds <- this_pomp$par_var_list$par_var_bounds  # add back in later
this_pomp$par_var_list$par_var_bounds <- NULL  # empty out the ranges
this_pomp$par_var_list$allparvals_update <- newparams


# Run final mif batch (2 rounds per replicate) ----------------------------

# redefine mif2 settings
mif_settings <- list()  # an empty list for mif2 settings
mif_settings$mif_num_particles  <- MIF_PARTICLES
mif_settings$mif_num_iterations <- MIF_ITERATIONS[3:4]
mif_settings$mif_cooling_fracs <- MIF_COOLING[3:4]
mif_settings$replicates <- MIF_RUNS
mif_settings$pf_num_particles <- PF_LL_PARTICLES
mif_settings$pf_reps <- PF_LL_REPS

# run batch 2 of mif
mif_res2 <- runmif_allstates(
  parallel_info = parallel_info,
  mif_settings = mif_settings,
  pomp_list = this_pomp,
  par_var_list = this_pomp$par_var_list
)

pomp_res2 <- this_pomp  # current state of affairs
pomp_res2$mif_res <- mif_res2  # store mif results

# generate traceplots and summarize fits
mif_explore2 <- exploremifresults(
  pomp_res = pomp_res2,
  par_var_list = pomp_res2$par_var_list,
  n_knots = n_knots
)


# add results computed to the pomp_res object
pomp_res2$traceplot <- mif_explore2$traceplot
pomp_res2$all_partable <- mif_explore2$all_partable
pomp_res2$est_partable <- mif_explore2$est_partable
pomp_res2$partable_natural <- mif_explore2$partable_natural

# batch 2 mif results
saveRDS(pomp_res2, file = paste0(outdir, state, "-mif_res_2.RDS"))


# Estimate smoothed posteriors of states for analysis ---------------------
# see Fox et al. PNAS, 2022 and, 
# https://github.com/kingaa/pomp/issues/74

# extract MLEs
params <- pomp_res2$all_partable %>%
  arrange(-LogLik) %>%
  slice(1) %>%
  dplyr::select(-MIF_ID, -LogLik, -LogLik_SE) %>%
  gather() %>%
  tibble::deframe()

foreach(i = 1:PF_STATES_REPS, 
        .packages = c("pomp")) %dopar% {
  pfilter(
    this_pomp$pomp_model,
    params = params,
    Np = PF_STATES_PARTICLES,
    filter.traj = TRUE,
    save.states = TRUE,
    max.fail = Inf
  )
} -> pf


# Combine the smoothed posterior state trajectories 
# for smoothed posterior distribution
filter_states <- tibble()
for(i in 1:length(pf)) {
  tmp <- pf[[i]]
  tmparr <- tmp@filter.traj[,1,]
  tmpout <- t(tmparr)
  colnames(tmpout) <- rownames(tmp@saved.states[[1]])
  tmpout <- as.data.frame(tmpout) %>%
    mutate(time = 1:n(),
           pf_rep = i,
           loglik = tmp@loglik) %>%
    dplyr::select(time, pf_rep, loglik, everything()) %>%
    as_tibble()
  filter_states <- bind_rows(filter_states, tmpout)
}

# save final states
# filtered state distribution
saveRDS(filter_states, file = paste0(outdir, state, "-filtered_states.RDS"))

# save pomp object list
# empty pomp model object
this_pomp$par_var_list$par_var_bounds <- bounds  # add back in
saveRDS(this_pomp, file = paste0(outdir, state, "-pomp_model.RDS"))

