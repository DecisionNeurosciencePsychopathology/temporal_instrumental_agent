library(R.matlab)
library(dplyr)
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "tc_identifiability_vetting"))
#tc <- readMat("tc_identifiability.mat")
tc <- readMat("tc_identifiability_focal_4000reps.mat")

origpars <- data.frame(tc$parsets)
names(origpars) <- paste0(c("lambda", "epsilon", "alphaG", "alphaN", "K", "nu", "rho"), "_orig")

fittedpars <- data.frame(tc$fittedpars)
names(fittedpars) <- paste0(c("lambda", "epsilon", "alphaG", "alphaN", "K", "nu", "rho"), "_fitted")

both <- cbind(origpars, fittedpars)

both$replication <- 1:nrow(both) #need this to uniquely identify observations

library(tidyr)
df <- gather(both, key, value, -replication) %>% 
  separate(key, c("parameter", "type"), sep="_") %>% spread(key=type, value=value)

library(ggplot2)
#pdf("tc identifiability.pdf", width=15, height=10)
#ggplot(df, aes(x=orig, y=fitted)) + geom_point(alpha=0.5) + 
#  facet_wrap(~parameter, nrow=2, scales="free") + theme_bw(base_size=15) +
#  stat_smooth(se=FALSE, color="blue", size=2)
#dev.off()

df$origfactor <- factor(df$orig)

pdf("tc identifiability focal.pdf", width=15, height=10)
ggplot(df, aes(x=origfactor, y=fitted)) + geom_jitter(alpha=0.4) + 
    facet_wrap(~parameter, nrow=2, scales="free") + theme_bw(base_size=15)
dev.off()

pdf("tc identifiability focal boxplots.pdf", width=15, height=10)
ggplot(df, aes(x=origfactor, y=fitted)) + geom_boxplot() + 
    facet_wrap(~parameter, nrow=2, scales="free") + theme_bw(base_size=15)
dev.off()



df %>% group_by(parameter) %>% summarize(cor(orig, fitted))

df %>% group_by(parameter, origfactor) %>% summarize(mean(fitted), sd(fitted))