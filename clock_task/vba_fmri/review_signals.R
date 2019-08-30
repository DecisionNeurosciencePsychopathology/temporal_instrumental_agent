setwd("~/Box Sync/skinner/projects_analyses/SCEPTIC/fMRI_paper/signals_review/")
library(readr)
library(lme4)
# library(lmerTest)
library(ggplot2)
library(tidyverse)
library(readr)
library(multcompView)
library(stargazer)


df <- read_csv("compiled_outputs/mmclock_fmri_decay_factorize_selective_psequate_mfx_trial_statistics.csv.gz")
df <- as.tibble(df)
p1 <- ggplot(df,aes(trial,ev,color = rewFunc)) + geom_line() + facet_wrap(~id)
ggsave("ev_fse_mfx.pdf", p1, width = 8, height = 6)

p2 <- ggplot(df,aes(trial,v_entropy,color = rewFunc)) + geom_line() + facet_wrap(~id)
ggsave("ventropy_fse_mfx.pdf", p2, width = 8, height = 6)
# some people just have entropy shelves after the first few trials -- must be low decay

p3 <- ggplot(df,aes(trial,d_auc,color = rewFunc)) + geom_line() + facet_wrap(~id)
ggsave("dauc_fse_mfx.pdf", p3, width = 8, height = 6)

p3a <- ggplot(df,aes(trial,sqrt(-d_auc),color = rewFunc)) + geom_line() + facet_wrap(~id)
ggsave("dauc_sqrt_fse_mfx.pdf", p3a, width = 8, height = 6)


# exactly -- the same subjects have flat decay curves
p4 <- ggplot(df,aes(trial,v_chosen,color = rewFunc)) + geom_line() + facet_wrap(~id)
ggsave("vchosen_fse_mfx.pdf", p4, width = 8, height = 6)
