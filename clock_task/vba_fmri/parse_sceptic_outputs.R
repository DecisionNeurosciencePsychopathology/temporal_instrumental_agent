parse_sceptic_outputs <- function(outdir, subjects_dir) {
  require(readr)
  require(dplyr)
  require(entropy)
  require(LaplacesDemon)
  sceptic_files <- list.files(path=outdir, pattern=".*sceptic.*\\.csv", full.names=TRUE)
  isexplore<-any(grepl("explore",sceptic_files))
  basis_file <- grep("basis.csv", sceptic_files, fixed=TRUE, value=TRUE)
  global_file <- grep("global_statistics.csv", sceptic_files, fixed=TRUE, value=TRUE)
  trial_file <- grep("trial_outputs_by_timestep.csv", sceptic_files, fixed=TRUE, value=TRUE)
  stopifnot(all(lengths(list(basis_file, global_file, trial_file)) == 1)) #only support one match per folder per now

  basis <- readr::read_csv(basis_file)
  global <- readr::read_csv(global_file)
  trial_df <- readr::read_csv(trial_file) %>% dplyr::rename(rt_next=u_1, score_next=u_2, run_boundary=u_3) %>% mutate(id=as.character(id))

  has_u <- any(grepl("^U_\\d+", names(trial_df), perl=TRUE)) #are u outputs present?
  has_d <- any(grepl("^D_\\d+", names(trial_df), perl=TRUE)) #are d outputs present?
  
  basis <- as.matrix(basis) #24 x 40 (nbasis x ntimesteps)

  #response vector (y1-y40)
  y_mat <- trial_df %>% dplyr::select(matches("y_\\d+")) %>% as.matrix()

  #identify the chosen time bin
  y_chosen <- apply(y_mat, 1, function(r) { which(r==1) })

  #####
  # Expected value statistics

  v_mat <- trial_df %>% dplyr::select(matches("V_\\d+")) %>% as.matrix() # (nsubjs * ntrials) x nbasis

  #put value back onto time grid
  v_func <- v_mat %*% basis

  #RT of the value max
  rt_vmax <- apply(v_func, 1, function(r) { ifelse(sd(r) < .1, NA, which.max(r)) })

  #Max value
  v_max <- apply(v_func, 1, function(r) { ifelse(sd(r) < .1, NA, max(r)) })

  #AUC of value function
  v_auc <- apply(v_func, 1, function(r) { ifelse(sd(r) < .1, NA, sum(r)) })

  #SD of value function
  v_sd <- apply(v_func, 1, function(r) {
    sdr <- sd(r)
    ifelse(sdr < .1, NA, sdr)
  })

  #value of the chosen time bin (action)
  v_chosen <- unname(sapply(1:nrow(v_func), function(r) { v_func[r, y_chosen[r]] }))

  #calculate entropy of the basis weights (as in the Cognition paper)
  v_entropy <- apply(v_mat, 1, function(basis_weights) {
    w_norm <- basis_weights/sum(basis_weights)
    nz <- w_norm[w_norm > 0]
    entropy <- -sum(nz * log10(nz))
    return(entropy)
  })

  #calculate entropy on the full value function, discretizing into 20 bins
  v_entropy_func <- apply(v_func, 1, function(r) {
    if (all(is.na(r)) || sd(r) < .01) { #essentially no variability, so hard to say there's a max
      return(NA)
    } else {
      #normalize entropy by discretizing into 20 bins
      dd <- discretize(r, numBins=20)
      entropy.empirical(dd)/log(20) #ML-based Shannon entropy, normalized to 0--1
    }
  })

  #loss functions based on K-L distance
  v_mat_lag <- trial_df %>% select(id, asc_trial) %>% bind_cols(as.data.frame(v_mat)) %>%
    mutate(asc_trial=as.integer(asc_trial), run_trial=case_when(
      asc_trial >= 1 & asc_trial <= 50 ~ asc_trial,
      asc_trial >= 51 & asc_trial <= 100 ~ asc_trial - 50L, #dplyr/rlang has gotten awfully picky about data types!!
      asc_trial >= 101 & asc_trial <= 150 ~ asc_trial - 100L,
      asc_trial >= 151 & asc_trial <= 200 ~ asc_trial - 150L,
      asc_trial >= 201 & asc_trial <= 250 ~ asc_trial - 200L,
      asc_trial >= 251 & asc_trial <= 300 ~ asc_trial - 250L,
      asc_trial >= 301 & asc_trial <= 350 ~ asc_trial - 300L,
      asc_trial >= 351 & asc_trial <= 400 ~ asc_trial - 350L,
      TRUE ~ NA_integer_)) %>%
    group_by(id) %>%
    mutate_at(vars(starts_with("V_")), funs(lag=lag(., 1, order_by=asc_trial))) %>% ungroup()

  kld_est <- sapply(1:nrow(v_mat_lag), function(r) {
    v_cur <- v_mat_lag %>% select(matches("V_\\d+$")) %>% slice(r) %>% unlist()
    v_lag <- v_mat_lag %>% select(matches("V_\\d+_lag$")) %>% slice(r) %>% unlist()
    if (sum(v_cur) < .001 || sum(v_lag) < .001) {
      kld <- c(NA, NA, NA, NA)
    } else {
      v_cur_norm <- v_cur/sum(v_cur)
      v_lag_norm <- v_lag/sum(v_lag)
      kl_out <- KLD(v_cur_norm,v_lag_norm)
      #KLD.px.py is distance from py (lag) to px (cur) here. New learning
      #KLD.py.px is distance from px (cur) to py (lag) here. Forgetting
      kld <- c(kl_out$mean.sum.KLD, kl_out$intrinsic.discrepancy,
               kl_out$sum.KLD.px.py, kl_out$sum.KLD.py.px)
    }
    return(kld)
  })

  kld_est <- t(kld_est) %>% as.data.frame() %>%
    setNames(c("mean_kld", "intrinsic_discrepancy", "kld_newlearn", "kld_forget"))

  #####
  # Prediction error statistics

  pe_mat <- trial_df %>% dplyr::select(matches("PE_\\d+")) %>% as.matrix()
  pe_func <- pe_mat %*% basis

  pe_max <- apply(pe_func, 1, function(r) {
    #note that this returns NA if there was no V in the time bin and an omission occurs
    #this should probably be a zero PE, not NA
    #if (sd(r) < .1) {
    #  return(NA)
    #} else if (min(r) >= 0) {

    if (min(r) >= 0) {
      return(max(r)) #positive PE
    } else {
      return(min(r)) #negative PE
    }
  })

  #Note that because PE is calculated as a function of the value curve, it's possible that the max/min PE is
  #not the PE of the chosen timestep. They tend to be highly correlated, but here's the PE of the chosen timestep.
  #Because the PE columns are right shifted by 1 compared to y_chosen, we need to align them here or we will not be grabbing
  #the right element of pe_func that corresponds to y_chosen
  y_chosen_prev <- trial_df %>% select(id, asc_trial) %>% bind_cols(y_chosen=y_chosen) %>% group_by(id) %>%
    mutate(y_chosen_prev=lag(y_chosen, order_by=asc_trial)) %>% pull(y_chosen_prev)

  pe_chosen <- sapply(1:nrow(pe_func), function(r) { pe_func[r, y_chosen_prev[r]] })

  #####
  # Decay statistics
  
  if (has_d) {
    
    d_mat <- trial_df %>% dplyr::select(matches("D_\\d+")) %>% as.matrix()
    if (ncol(d_mat) > 0) { #will be absent for non-decay models
      d_func <- d_mat %*% basis
      d_auc <- apply(d_mat, 1, function(r) {
        if (sd(r) < .001) {
          return(NA)
        } else {
          return(sum(r))
        }
      })
    } else {
      d_auc <- NULL #so that bind_cols proceeds
    }
  }

  #####
  # Uncertainty statistics
  if (has_u) {

    u_mat <- trial_df %>% dplyr::select(matches("U_\\d+")) %>% as.matrix() # (nsubjs * ntrials) x nbasis
    
    #put value back onto time grid
    u_func <- u_mat %*% basis

    #uncertainty of the chosen time bin (action)
    u_chosen <- unname(sapply(1:nrow(u_func), function(r) { u_func[r, y_chosen[r]] }))
  }

  
  #For PE- AND D-based statistics, these have been right-shifted by 1 position (in u) to line up the inputs and outputs in SCEPTIC VBA fitting.
  #This is necessary to get the t versus t-1 in the evolution of learning in VBA.

  #Consequently, the PE of trial 1 is actually in position 2, etc. And, for now, the PE and D of the last trial (50) is not collected/estimated

  #compile trial statistics
  trial_stats <- trial_df %>% select(id, dataset, model, asc_trial, rt_next, score_next) %>%
    bind_cols(y_chosen=y_chosen, v_chosen=v_chosen, rt_vmax=rt_vmax, v_max=v_max, v_auc=v_auc,
              v_sd=v_sd, v_entropy=v_entropy, v_entropy_func=v_entropy_func, pe_max=pe_max, pe_chosen=pe_chosen) %>%
    group_by(id) %>%
    mutate(
      rt=lead(rt_next, 1, order_by=asc_trial), #shift rt back onto original time grid
      score=lead(score_next, 1, order_by=asc_trial), #same for score
      pe_max=lead(pe_max, 1, order_by=asc_trial), #same for PE_max
      pe_chosen=lead(pe_chosen, 1, order_by=asc_trial) #same for PE_chosen
    ) %>% ungroup() %>% dplyr::rename(rt_vba=rt, score_vba=score) %>% cbind(kld_est)

  if (has_u) {
    trial_stats <- trial_stats %>% bind_cols(u_chosen=u_chosen) %>% group_by(id) %>%
      mutate(
        u_chosen_lag = lag(u_chosen, 1, order_by=asc_trial),
        u_chosen_change = u_chosen - u_chosen_lag
      ) %>% ungroup()
  }

  if (has_d) {
    trial_stats <- trial_stats %>% bind_cols(d_auc=d_auc) %>% group_by(id) %>%
      mutate(d_auc=lead(d_auc, 1, order_by=asc_trial)) %>% ungroup() #as with RTs and PE, need to left-shift DAUC
  }
  
  #remove readr detritus
  attr(trial_stats, "spec") <- NULL

  #merge with original data
  if(isexplore){
    subj_files <- list.files(path=subjects_dir, pattern=".*\\.csv", full.names=TRUE, recursive=TRUE)
  } else {
    subj_files <- list.files(path=subjects_dir, pattern=".*\\.csv", full.names=TRUE, recursive=FALSE)
  }

  subj_data <- bind_rows(lapply(subj_files, function(ff) {
    df <- read_csv(ff) %>% dplyr::rename(rt_csv=rt, score_csv=score) %>% mutate(asc_trial=1:n()) #asc_trial used to match with vba outputs
    if(isexplore){
      a<-strsplit(ff,.Platform$file.sep)[[1]]
      df$id<-a[length(a)-1]
    } else {
      df$id <- sub(".*(?<=MEG_|fMRIEmoClock_)([\\d_]+)(?=_tc|_concat).*", "\\1", ff, perl=TRUE)
    }
    df <- df %>% select(id, run, asc_trial, rewFunc, emotion, everything())
    return(df)
  }))

  trial_stats <- trial_stats %>% left_join(subj_data, by=c("id", "asc_trial")) %>%
    select(dataset, model, id, run, trial, asc_trial, rewFunc, emotion, rt_csv, score_csv, magnitude, probability, ev, everything()) %>%
    arrange(dataset, model, id, trial)

  return(trial_stats)
}
