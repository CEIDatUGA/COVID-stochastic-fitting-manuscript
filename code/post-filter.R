# code to conduct post-fitting particle filtering to yield smoothed posterior
# states for runs that failed on the GACRC.

library(tidyverse)
library(tidyr)
library(tibble)
library(pomp)  # must be v3.1
# devtools::install_version("pomp", version = "3.1", repos = "http://cran.us.r-project.org")
library(foreach)
library(doParallel)
library(parallel)

PF_REPS <- 500
PF_PARTICLES <- 2500
NUM_CORES <- 17


states_to_filter <- state.name[1:17]  # states with fits but no filters

mif_path <- "L:/1186 - UGA/MIFS/"
all_mif_files <- list.files(mif_path)[grep("res_2", list.files(mif_path))]

registerDoParallel(NUM_CORES)
foreach(i = 1:length(states_to_filter), 
        .packages = c("pomp", "dplyr", "tibble", "tidyr")) %dopar% {
          # load fit
          fname <- all_mif_files[grep(states_to_filter[i], all_mif_files)]
          pomp_res2 <- readRDS(paste0(mif_path, fname))
          
          # extract MLEs
          params <- pomp_res2$all_partable %>%
            arrange(-LogLik) %>%
            slice(1) %>%
            dplyr::select(-MIF_ID, -LogLik, -LogLik_SE) %>%
            tidyr::gather() %>%
            tibble::deframe()
          
          out <- tibble::tibble()
          for(j in 1:PF_REPS){
            # run filter
            pf <- pfilter(
              pomp_res2$pomp_model,
              params = params,
              Np = PF_PARTICLES,
              filter.traj = TRUE,
              save.states = TRUE,
              max.fail = Inf
            )
            tmparr <- pf@filter.traj[,1,]
            tmpout <- t(tmparr)
            colnames(tmpout) <- rownames(pf@saved.states[[1]])
            tmpfilter <- as.data.frame(tmpout) %>%
              dplyr::mutate(time = 1:n(),
                            pf_rep = j,
                            loglik = pf@loglik) %>%
              dplyr::select(time, pf_rep, loglik, dplyr::everything()) %>%
              tibble::as_tibble()
            out <- dplyr::bind_rows(out, tmpfilter)
          }
          saveRDS(out, paste0("../output/",state.name[i],"-filtered_states.RDS"))
        }
stopImplicitCluster()
