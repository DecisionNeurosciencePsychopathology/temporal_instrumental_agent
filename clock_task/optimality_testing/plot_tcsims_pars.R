library(R.matlab)
library(ggplot2)
library(tidyr)
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "optimality_testing"))

tcsimfiles <- list.files(pattern="tcsims_frank09.*\\.mat")

#sims <- c("good", "nu0p1", "alphan1p5", "lambda0p5")
#simnames <- c("Baseline", "nu = 0.1", "alphaN = 1.5", "lambda = 0.5")

sims <- c("good", "rho10k", "alphan1p0", "lambda0p5")
simnames <- c("Baseline", "rho = 10,000", "alphaN = 1.0", "lambda = 0.5")


simmat <- c()
for (s in 1:length(sims)) {
  pars <- readMat(paste0("tcsims_frank09_", sims[s], ".mat"))

  rtpred <- data.frame(t(pars$allRTsmoothavg))
  names(rtpred) <- c('DEV',    'IEV',    'CEV',    'CEVR')
  rtpred$trial <- 1:nrow(rtpred)
  rtpred$sim <- simnames[s]
  simmat <- rbind(simmat, rtpred)
}

#good <- readMat("tcsims_good.mat")
##good <- readMat("tcsims_franksugg.mat")
#
#rtpred <- data.frame(t(good$allRTsmoothavg))
#names(rtpred) <- c('DEV',    'IEV',    'CEV',    'CEVR')
#rtpred$trial <- 1:nrow(rtpred)
#rtpred$sim <- 'good'

#bad alphas
#good <- readMat("tcsims_good.mat")
#
#rtpred <- data.frame(t(good$allRTsmoothavg))
#names(rtpred) <- c('DEV',    'IEV',    'CEV',    'CEVR')
#rtpred$trial <- 1:nrow(rtpred)
#rtpred$sim <- 'good'


m <- simmat %>% gather('contingency', 'rt', -trial, -sim)
m$sim <- ordered(m$sim, levels=c(simnames))


pdf("tcsims_50trials.pdf", width=10, height=10)
ggplot(subset(m, trial <= 50), aes(x=trial, y=rt, color=contingency)) + geom_line(size=1.5) + facet_wrap(~sim, ncol=2) + theme_bw(base_size=20) +
    xlab("Trial") + ylab("Average response time (ms)") + scale_color_brewer("Contingency", palette="Set1") +
    theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)))
dev.off()

pdf("tcsims_500trials.pdf", width=10, height=10)
ggplot(m, aes(x=trial, y=rt, color=contingency)) + geom_line(size=1.5) + facet_wrap(~sim, ncol=2) + theme_bw(base_size=20) +
    xlab("Trial") + ylab("Average response time (ms)") + scale_color_brewer("Contingency", palette="Set1") +
    theme(axis.title.x=element_text(margin = margin(t = 20)), axis.title.y=element_text(margin = margin(r = 20)))
dev.off()