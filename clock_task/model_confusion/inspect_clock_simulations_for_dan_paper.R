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
