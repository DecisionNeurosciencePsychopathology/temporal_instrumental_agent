parse_sceptic_outputs <- function(outdir, subjects_dir) {
  library(readr)
  library(dplyr)
  sceptic_files <- list.files(path=outdir, pattern=".*sceptic.*\\.csv", full.names=TRUE)
  basis_file <- grep("basis.csv", sceptic_files, fixed=TRUE, value=TRUE)
  global_file <- grep("global_statistics.csv", sceptic_files, fixed=TRUE, value=TRUE)
  trial_file <- grep("trial_outputs_by_timestep.csv", sceptic_files, fixed=TRUE, value=TRUE)
  stopifnot(all(lengths(list(basis_file, global_file, trial_file)) == 1)) #only support one match per folder per now

  basis <- readr::read_csv(basis_file)
  global <- readr::read_csv(global_file)
  trial <- readr::read_csv(trial_file) %>% dplyr::rename(rt_next=u_1, score_next=u_2) %>% mutate(id=as.character(id))

  basis <- as.matrix(basis) #24 x 40 (nbasis x ntimesteps)

  #response vector (y1-y40)
  y_mat <- trial %>% dplyr::select(matches("y_\\d+")) %>% as.matrix()

  #identify the chosen time bin
  y_chosen <- apply(y_mat, 1, function(r) { which(r==1) })

  #####
  # Expected value statistics
  
  v_mat <- trial %>% dplyr::select(matches("V_\\d+")) %>% as.matrix()

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
  v_chosen <- unname(sapply(1:nrow(v_func), function(r) {
    v_func[r, y_chosen[r]]
  }))

  #####
  # Prediction error statistics

  pe_mat <- trial %>% dplyr::select(matches("PE_\\d+")) %>% as.matrix()
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

  #####
  # Decay statistics
  d_mat <- trial %>% dplyr::select(matches("D_\\d+")) %>% as.matrix()
  if (ncol(d_mat) > 0) { #will be absent for non-decay models
    d_func <- d_mat %*% basis
    d_auc <- apply(d_mat, 1, function(r) {
      if (sd(r) < .1) {
        return(NA)
      } else {
        return(sum(r))
      }
    })    
  } else {
    d_auc <- NULL #so that bind_cols proceeds
  }

  #For PE- AND D-based statistics, these have been right-shifted by 1 position (in u) to line up the inputs and outputs in SCEPTIC VBA fitting.
  #This is necessary to get the t versus t-1 in the evolution of learning in VBA.
  
  #Consequently, the PE of trial 1 is actually in position 2, etc. And, for now, the PE and D of the last trial (50) is not collected/estimated
  #Thus, flip the first trial (value of 0) to the end (currently unspecified/not estimated)

  #compile trial statistics
  trial_stats <- trial %>% select(id, dataset, model, trial, rt_next, score_next) %>%
    bind_cols(y_chosen=y_chosen, v_chosen=v_chosen, rt_vmax=rt_vmax, v_max=v_max, v_auc=v_auc,
      v_sd=v_sd, pe_max=pe_max, d_auc=d_auc) %>%
    group_by(id) %>%
    mutate(
      rt=lead(rt_next, 1, order_by=trial), #shift rt back onto original time grid
      score=lead(score_next, 1, order_by=trial), #same for score
      pe_max=lead(pe_max, 1, order_by=trial), #same for PE
      d_auc=lead(d_auc, 1, order_by=trial) #same for DAUC
    ) %>% ungroup() %>% dplyr::rename(rt_vba=rt, score_vba=score)

  #remove readr detritus
  attr(trial_stats, "spec") <- NULL

  #merge with original data
  subj_files <- list.files(path=subjects_dir, pattern=".*\\.csv", full.names=TRUE, recursive=FALSE)
  subj_data <- bind_rows(lapply(subj_files, function(ff) {
    df <- read_csv(ff) %>% dplyr::rename(rt_csv=rt, score_csv=score) %>% mutate(asc_trial=1:n()) #asc_trial used to match with vba outputs
    df$id <- sub(".*(?<=MEG_|fMRIEmoClock_)([\\d_]+)(?=_tc|_concat).*", "\\1", ff, perl=TRUE)
    df <- df %>% select(id, run, trial, rewFunc, emotion, everything())
    return(df)
  }))

  trial_stats <- trial_stats %>% left_join(subj_data, by=c("id", "trial")) %>%
    select(dataset, model, id, run, trial, asc_trial, rewFunc, emotion, rt_csv, score_csv, magnitude, probability, ev, everything()) %>%
    arrange(dataset, model, id, trial)
  
  return(trial_stats)
}
