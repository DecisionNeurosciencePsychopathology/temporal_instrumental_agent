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

ev_list <- list()
for (contingency_index in 1:n_contingencies) {
  cont_df <- as_tibble(as.vector(optmat[,contingency_index]$ev)) %>% rename(ev_chosen = value)
  cont_df$contingency <- contingency_index
  cont_df$contingency_name <- as.character(optmat[1,][contingency_index])
  cont_df$rt <- as.numeric(1:500)
  ev_list[[contingency_index]] <- cont_df
}
ev_df <- data.table::rbindlist(ev_list)

ggplot(bdf, aes(rt, reward)) + geom_smooth() + facet_wrap(~contingency)


bdf <- bdf %>% inner_join(ev_df, by = c("rt", "contingency")) %>% 
  # remove regular KF, dead on arrival
  filter(!str_detect(pattern = "kalman",model))%>% 
  group_by(model, contingency, set) %>% arrange(trial, by_group = T) %>% 
  mutate(ev_lag = lag(ev_chosen),
         rt_lag = lag(rt),
         reward_lag = lag(reward)) %>% ungroup() 
ggplot(ev_df, aes(rt, ev_chosen)) + geom_line() + facet_wrap(~contingency)
ggplot(bdf, aes(rt, ev_chosen)) + geom_line() + facet_wrap(~contingency)

ggplot(ev_df, aes(rt, ev_chosen)) + geom_smooth()

setwd("~/OneDrive - University of Pittsburgh/Documents/SCEPTIC_fMRI/model_simulations/")
pdf("ev_chosen_by_model_trial_set.pdf", height = 30, width = 10)
ggplot(bdf, aes(trial,ev_chosen, color = model)) + geom_smooth() + facet_wrap(contingency~set)
dev.off()

pdf("ev_chosen_by_model.pdf", height = 4, width = 5)
ggplot(bdf, aes(trial,ev_chosen, color = model)) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs")) + facet_wrap(~set)
dev.off()

pdf("rt_swing_by_model.pdf", height = 4, width = 5)
ggplot(bdf, aes(trial,rt_swing, color = model)) + 
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5, bs = "cs")) + facet_wrap(~set)
dev.off()

# sanity check
# ggplot(bdf, aes(rt,ev_chosen, color = model)) + geom_smooth(method = "loess") + facet_wrap(contingency~set)

pdf("ev_chosen_vs_rtswing_by_model.pdf", height = 4, width = 5)
ggplot(bdf, aes(rt_swing,ev_chosen, color = model)) + geom_smooth(method = "loess") + facet_wrap(~set)
dev.off()

pdf("rt_swing_ev_by_model_and_reward.pdf", height = 4, width = 8)
ggplot(bdf %>% filter(trial>1), aes(rt_swing,ev_chosen, color = model)) + geom_smooth(method = "loess") + facet_wrap(~reward_lag>0)# + facet_wrap(contingency~set)
dev.off()
