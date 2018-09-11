setwd("~/Box Sync/skinner/projects_analyses/SCEPTIC/compiled_sceptic_vba_datasets/")
library(readr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(corrplot)
library(ggpubr)

# meaning fMRI decay subject
fds <- read_csv("mmclock_fmri_decay_mfx_sceptic_global_statistics.csv")
mds <- read_csv("mmclock_meg_decay_mfx_sceptic_global_statistics.csv")

# "fMRI decay trial"
fdt <- read_csv("mmclock_fmri_decay_mfx_trial_stats.csv.gz")
mdt <- read_csv("mmclock_meg_decay_mfx_trial_stats.csv.gz")

# "meg basis"
mb <- read_csv("mmclock_meg_decay_mfx_sceptic_basis.csv")

md <- merge(mds,mdt)
fd <- merge(fds,fdt)

# figure out missing RTs
plot(md$rt_vba,(md$isi_onset-md$clock_onset))
md$pred_rt <- md$isi_onset-md$clock_onset
boxplot(md$pred_rt - md$rt_csv~md$id)
boxplot(md$rt_vba*100 - md$rt_csv~md$id)
boxplot(md$rt_vba*100 - md$pred_rt~md$id)

unique(md$id[abs(md$rt_vba*100 - md$rt_csv)>100])

boxplot(md$rt_csv~md$id)

plot(md$rt_csv,md$rt_vba)

describe(md$pred_rt - md$rt_csv)
plot(md$trial, md$rt_csv - md$isi_onset + md$clock_onset)

unique(md$id[md$rt_csv!=(md$isi_onset-md$clock_onset)])
unique(md$id[md$rt_csv==(md$isi_onset-md$clock_onset)])



z <- as.tibble(md[md$R2>.2,])
# parameter scatterplots
scatter.smooth(fds$alpha,fds$gamma)
scatter.smooth(mds$alpha,mds$gamma)

scatter.smooth(fds$R2,fds$beta)
scatter.smooth(mds$R2,mds$beta)

setwd("plots")

mds_params <-  mds[,c("alpha", "gamma", "beta", "R2")]
#val_rois <- val_rois[,-grep("ACC",names(val_rois))]
mcormat <- (cor(mds_params))
p1 <- corrplot(mcormat,  type = "upper")

fds_params <-  fds[,c("alpha", "gamma", "beta", "R2")]
#val_rois <- val_rois[,-grep("ACC",names(val_rois))]
fcormat <- (cor(fds_params))
p2 <- corrplot(fcormat,  type = "upper")

# what is wrong with the cluster of low-temp, R2~1 subjects?
md = mdt %>% group_by(id) %>% arrange(id, run, trial) %>% dplyr::summarize(
  totRew = mean(na.omit(score_vba))) %>% ungroup()
mds <- merge(mds,md)
scatter.smooth(mds$R2,mds$totRew)
mdt <- merge(mdt,mds)
mdt$perfectR2 <- mdt$R2>.3
p <- ggplot(mdt, aes(x = trial, y = rt_csv)) + geom_point() + facet_wrap(~perfectR2)
ggsave("perfect_R2_diag.pdf",p)
# they don't have any RTs

# Time-courses of decay, entropy, PE, v_max

fdt = fdt %>% group_by(id) %>% arrange(id, trial) %>% dplyr::mutate(
  meanRew = mean(na.omit(score_vba))) %>% ungroup()

fdt$aboveMedian <- fdt$meanRew > median(fdt$meanRew)

p <- ggplot(fdt, aes(x = trial, y = d_auc, color = rewFunc)) + geom_point() + 
  facet_wrap(~meanRew>median(fdt$meanRew), labeller = label_both)
ggsave("fmri_decay_timecourse.pdf",p)


p <- ggplot(fdt, aes(x = trial, y = pe_max, color = rewFunc)) + geom_point() + 
  facet_wrap(~meanRew>median(fdt$meanRew), labeller = label_both)
ggsave("fmri_PEmax_timecourse.pdf",p)

p1 <- ggplot(fdt[fdt$aboveMedian,], aes(x = trial, y = pe_max, color = rewFunc)) + geom_point() + 
  facet_wrap(~id)
p2 <- ggplot(fdt[!fdt$aboveMedian,], aes(x = trial, y = pe_max, color = rewFunc)) + geom_point() + 
  facet_wrap(~id)
pall <- ggarrange(p1,p2, labels = c("good", "bad"))
ggsave("fmri_PEmax_timecourse_by_subject.pdf",pall, width = 20, height = 20)


