library(tidyverse)
library(R.matlab)
l <- readMat("~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/model_simulations/confusion_out.mat")
results_l <- l$oresults

models <- unique(unlist(l$fitmodelnames))
model_list = list()
n_contingencies = 40
n_trials = 50
n_sets = 3
# mdf = as_tibble(NULL)
for (model_index in 1:length(models)) {
  model_name <- unlist(names(results_l[model_index,,]))
  # get costs, rts, rewards
  m <- results_l[model_index,,][[unlist(names(results_l[model_index,,]))]]
  # costs <- m[1,,][[1]]
  rts <- m[3,,][[1]][1,,] # 40 contingencyjects by 3 (param sets), each containing a list of 50 trials
  rewards <- m[4,,][[1]][1,,] # 40 contingencyjects by 3 (param sets), each containing a list of 50 trials
  set_list <- list()
  contingency_list <- list()
  for (contingency in 1:n_contingencies) {
    for (set in 1:n_sets) {
      df <- as_tibble(unlist(rts[contingency,set][[1]])) %>% rename(rt = value)
      df$rt_swing <- c(NA, abs(diff(df$rt)))
      df$trial <- 1:n_trials
      df$contingency <- contingency
      df$set <- set
      df$reward <- unlist(rewards[contingency,set][[1]])
      df$model <- model_name
      set_list[[set]] <- df
    }
    contingency_list[[contingency]] <- data.table::rbindlist(set_list)
  }
  mdf <- data.table::rbindlist(contingency_list)
  model_list[[model_index]] <- mdf
}
bdf <- data.table::rbindlist(model_list)

# we know that betas are the same across models

ggplot(bdf, aes(trial, rt, lty = model)) + geom_smooth()
ggplot(bdf %>% filter(str_detect(model, "fixed")), aes(trial, rt_swing, color = model)) + geom_smooth(se = F) + geom_line(alpha = .1) + facet_grid(~set)
ggplot(bdf, aes(trial, reward, color = model)) + geom_smooth() + facet_wrap(~set)

# get the sinusoidals
optmat <- l$optmat[,1,]
#get the EV of the models' choices (not just the score)
# idea: loop over all contingencies
models <- unique(bdf$model)
for (contingency_index in 1:n_contingencies) {
  bdf$sinusoidal[bdf$contingency==contingency_index] <- as.character(optmat[,contingency_index]$name)
  for (model_index in 1:length(models)) {
    for (set in 1:n_sets) {
      for (trial in 1:n_trials) {
        bdf$ev_chosen[bdf$contingency==contingency_index & bdf$model==models[model_index] & bdf$set==set & bdf$trial==trial] <- 
          optmat[,contingency]$ev[bdf$rt[bdf$contingency==contingency_index & bdf$model==models[model_index] & bdf$set==set & bdf$trial==trial]]
      }
    }
  }
}
bdf <- bdf %>% group_by(model, contingency, set) %>% arrange(trial, by_group = T) %>% mutate(ev_lag = lag(ev_chosen),
                                                                                             rt_lag = lag(rt)) %>% ungroup()
setwd("~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/model_simulations/")
pdf("ev_chosen_by_model_trial_set.pdf", height = 30, width = 10)
ggplot(bdf, aes(trial,ev_chosen, color = model)) + geom_smooth() + facet_wrap(contingency~set)
dev.off()

# sanity check: SOMETHING WRONG
ggplot(bdf, aes(rt,ev_chosen, color = model)) + geom_smooth(method = "loess")# + facet_wrap(contingency~set)

pdf("ev_chosen_vs_rtswing_by_model.pdf", height = 30, width = 10)
ggplot(bdf, aes(rt_swing,ev_chosen, color = model)) + geom_smooth(method = "loess")# + facet_wrap(contingency~set)
dev.off()

ggplot(bdf %>% filter(trial>1), aes(rt_swing,ev_chosen, color = model)) + geom_smooth(method = "loess") + facet_wrap(~ev_lag > 22)# + facet_wrap(contingency~set)
