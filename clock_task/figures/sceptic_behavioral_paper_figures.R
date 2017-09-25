##just a quick plot of the basis
library(ggplot2)
library(dplyr)
ntimesteps=13
maxt=4000
centers <- seq(0,4000, length=ntimesteps)
gaussmat <- sapply(centers, function(v) {
    dnorm(x=0:maxt, mean=v, sd=300)
})

#matplot(gaussmat, type="l", lty=1, lwd=5, col=colors(distinct=TRUE))
pdf("gauss_basis.pdf", width=6, height=5)
matplot(gaussmat, type="l", lty=1, lwd=5, col=colorRampPalette(c("blue", "red"))(13), ann=FALSE, xaxt='n', yaxt='n', bty='n')
dev.off()

##test out weights and centers to show example
weights<- 100*c(0.0776, 7.4801, 2.0792, 1.4008, 2.0000, 1.0000, 3.0000, 9.0000, 25.0000, 22.0024, 30.2286, 12.0326, 30.2332)
v <- 5* apply( sapply(1:nrow(gaussmat), function(r) {
    gaussmat[r,]*weights
}), 2, sum)



pdf("learned_ev.pdf", width=6, height=5)
df <- data.frame(time=seq(0,4,length=4001), ev=v)
ggplot(df, aes(x=time, y=ev)) + geom_line(size=1.5) + ylab("Expected value (learned)\n") + xlab("\nTime (s)") + theme_bw(base_size=24)
dev.off()


##uncertainty weights
uweights<- -10*c(0.3776, 1.4801, 0.9792, 1.4008, 2.0000, 1.0000, 3.0000, 5.0000, 7.0000, 8.0024, 10.2286, 3.0326, 1.2332)
u <- 5* apply( sapply(1:nrow(gaussmat), function(r) {
    gaussmat[r,]*uweights
}), 2, sum)


library(ggplot2)
pdf("learned_uncertainty.pdf", width=6, height=5)
df <- data.frame(time=seq(0,4,length=4001), u=u)
ggplot(df, aes(x=time, y=u)) + geom_line(size=1.5) + ylab("Uncertainty (experienced)\n") + xlab("\nTime (s)") + theme_bw(base_size=24)
dev.off()

    
pdf("weights.pdf", width=6, height=3.5)
par(mar=c(1, 5, 2, 1) + 0.1)

plot(c(0, 4000), c(0, max(weights) + 50), type = "n", xlab = "", ylab = "Basis weight (AU)", yaxs="i",
     cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, bty='n', xaxt='n')
rect(centers-50, 0, centers+50, weights, col="gray")
#axis(1, at=c(0,4), labels=c("0s", "4s"))#, pos=, lty=, col=, las=, tck=, ...)
dev.off()


##true underlying IEV contingency.
setwd("~/Data_Analysis/clock_analysis/td_model")
source("~/Data_Analysis/clock_analysis/td_model/getrew.R")
fm <- getMagFreq(0:4000, "IEV")
f <- fm[4002:8002]
m <- fm[1:4001]
ev <- f*m

library(ggplot2)
pdf("true_ev.pdf", width=6, height=5)
df <- data.frame(time=seq(0,4,length=4001), ev=ev)
ggplot(df, aes(x=time, y=ev)) + geom_line(size=1.5) + ylab("Expected value\n") + xlab("\nTime (s)") + theme_bw(base_size=24)
dev.off()

##real curve
pdf("iev_func.pdf", width=5, height=4)
ggplot(df, aes(x=time, y=ev)) + geom_line(size=1.5) + ylab("Expected value\n") + xlab("\nTime (ms)") + theme_bw(base_size=24)
dev.off()

df <- c()
#take CEVR out since CEV and CEVR are identical wrt EV (and we are not showing prob + freq)
for (cont in c("IEV", "DEV", "CEV", "CEVR")) { #, "CEVR"
  fm <- getMagFreq(0:4000, cont)
  #if (cont=="CEV") { cont="CEV/CEVR" } #for plot name
  df <- rbind(df, data.frame(contingency=cont, time=0:4000, mag=fm$Mag, freq=fm$Freq, ev=fm$Mag*fm$Freq))
}

#Figure: plot of EV in clock task
pdf("Clock contingencies.pdf", width=5, height=3.4)
ggplot(df, aes(x=time/1000, y=ev, color=contingency)) + geom_line(size=2) + ylab("Expected value (points)") + xlab("Time (seconds)") + scale_color_brewer("Contingency", palette="Dark2") +
    theme_bw(base_size=18) + theme(axis.title.x=element_text(margin = margin(t = 10)), axis.title.y=element_text(margin = margin(r = 10)), legend.margin = unit(0.15, "cm"), 
        plot.margin=margin(r=3, l=3, t=10, b=10))
dev.off()

#freq, prob, and ev
library(cowplot)
gcommon <- list(geom_line(size=2), xlab("Time (seconds)"), scale_color_brewer("Contingency", palette="Dark2"), theme_bw(base_size=18), 
  theme(axis.title.x=element_text(margin = margin(t = 10)), axis.title.y=element_text(margin = margin(r = 8)), 
        legend.margin = unit(0.15, "cm"), 
        plot.margin=margin(r=10, l=10, t=10, b=5)))


g1 <- ggplot(df, aes(x=time/1000, y=mag, color=contingency)) + gcommon + ylab("Reward magnitude (points)") + 
  theme(legend.position="none", plot.margin=margin(r=10, l=0, t=10, b=5))

g2 <- ggplot(df, aes(x=time/1000, y=freq, color=contingency)) + gcommon + ylab("Reward probability") + theme(legend.position="none")

df$ev[df$contingency=="CEVR"] <- df$ev[df$contingency=="CEVR"] + 0.5 #offset for plotting
df$ev[df$contingency=="CEV"] <- df$ev[df$contingency=="CEV"] - 0.5 #offset for plotting
g3 <- ggplot(df, aes(x=time/1000, y=ev, color=contingency)) + gcommon + ylab("Expected value (points)") + theme(legend.position="none")

pdf("Clock contingencies with freq mag.pdf", width=8.5, height=4)
pg <- plot_grid(g1, g2, g3, nrow=1)
plot(pg)
dev.off()

pdf("Clock contingencies legend.pdf", width=2, height=2)
legend_b <- get_legend(g1 + theme(legend.position="right") + theme_bw(base_size=25) + guides(color = guide_legend(keywidth=2, keyheight=2)))

p <- plot_grid(legend_b, nrow=1)
plot(p)
dev.off()


# add the legend underneath the row we made earlier. Give it 10% of the height
# of one plot (via rel_heights).
# p <- plot_grid(pg, legend_b, nrow=1, rel_widths = c(.9, .15))
# plot(p)


setwd("~/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri")
#load(file="dataframe_for_entropy_analysis_Oct2016.RData")
#this contains data with 24 basis functions and post-Niv learning rule
#load(file="dataframe_for_entropy_analysis_Nov2016.RData")
load(file="dataframe_for_entropy_analysis_Mar2017.RData") #has the random priors entropy

bdf = bdf %>% rename(subject=rowID) %>% group_by(subject) %>% arrange(subject, run, trial) %>% mutate(totreward=sum(score), cumreward=cumsum(score)) %>% ungroup() %>%
  mutate(medreward=median(totreward), #between subjects
         msplit=factor(as.numeric(totreward > medreward), levels=c(0,1), labels=c("< Median", ">= Median")))

bdf <- bdf %>% mutate(msplit=recode(msplit, "< Median" = "Total~earnings<median", ">= Median"="Total~earnings>=median"))
bdf$rewFunc <- factor(bdf$rewFunc, levels=c("IEV", "DEV", "CEV", "CEVR")) #to match contingency plot

# Figure 1c
pdf("Fig_1c.pdf", width = 8.5, height = 3.75)
ggplot(bdf, aes(x=trial, y=rt/1000, color = rewFunc )) + stat_smooth(method="loess", size = 2) +
  scale_color_brewer("Contingency", palette="Dark2") +
  theme_bw(base_size=18) + facet_wrap(~msplit, labeller=label_parsed) + xlab("Trial") +
  ylab("Response time (seconds)") + scale_y_continuous(breaks=c(1.25, 1.5, 1.75, 2, 2.25)) +
  theme(axis.title.x=element_text(margin = margin(t = 12)), axis.title.y=element_text(margin = margin(r = 12)), 
        legend.margin = margin(t=0, r=2, b=0, l=5), plot.margin=margin(r=10, l=10, t=10, b=5)) + theme(legend.position="none")
  
dev.off()

# 1d
# pdf("Fig_1d.pdf", width = 10, height = 4)
# ggplot(subset(bdf), aes(x=trial, y=abstschange*100, color = rewFunc)) + stat_smooth(method="loess", size = 2) + theme_bw(base_size=25) + facet_wrap(~msplit) + ylab("RT swings, ms") + labs(colour = "Contingency") #facet_wrap(~msplit) #geom_jitter(alpha=0.2) +
# dev.off()

pdf("Fig_1d.pdf", width = 8.5, height = 3.75)
ggplot(bdf, aes(x=trial, y=abstschange/10, color = rewFunc )) + stat_smooth(method="loess", size = 2) +
  scale_color_brewer("Contingency", palette="Dark2") +
  theme_bw(base_size=18) + facet_wrap(~msplit, labeller=label_parsed) + xlab("Trial") +
  ylab("Change in RT (seconds)") +
  theme(axis.title.x=element_text(margin = margin(t = 12)), axis.title.y=element_text(margin = margin(r = 12)), 
        legend.margin = margin(t=0, r=2, b=0, l=5), plot.margin=margin(r=10, l=10, t=10, b=5)) + theme(legend.position="none")
dev.off()




#Supplementary figure: Sinusoidal contingency
ntimesteps=500

ev = 10*sin(2*pi*(1:ntimesteps)*1/ntimesteps) + 2.5*sin(2*pi*(1:ntimesteps)*2/ntimesteps) + 2.0*cos(2*pi*(1:ntimesteps)*4/ntimesteps)
ev = ev + abs(min(ev)) + 10;
prb = 25*cos(2*pi*(1:ntimesteps)*1/ntimesteps) + 10*cos(2*pi*(1:ntimesteps)*3/ntimesteps) + 6*sin(2*pi*(1:ntimesteps)*5/ntimesteps)
prb_max=0.7
prb_min=0.3
prb = (prb - min(prb))*(prb_max-prb_min)/(max(prb)-min(prb)) + prb_min

allshift = array(NA_real_, dim=c(ntimesteps, ntimesteps, 3))

for (i in 1:ntimesteps) {
  if (i > 1) {
    shift = c(i:ntimesteps, 1:(i-1))
  } else { shift <- 1:ntimesteps }
  evi = ev[shift]
  prbi = prb[shift]
  
  allshift[i,,1] = evi
  allshift[i,,2] = prbi
  allshift[i,,3] = evi/prbi
  
}

shift3 <- rbind(data.frame(time=1:ntimesteps/100, EV=allshift[,1,1], Probability=allshift[,1,2], Magnitude=allshift[,1,3], name="shift = 0"),
    data.frame(time=1:ntimesteps/100, EV=allshift[,100,1], Probability=allshift[,100,2], Magnitude=allshift[,100,3], name="shift = 100"),
    data.frame(time=1:ntimesteps/100, EV=allshift[,200,1], Probability=allshift[,200,2], Magnitude=allshift[,200,3], name="shift = 200"))

library(reshape2); library(ggplot2)
m3 <- melt(shift3, id.vars=c("time", "name"))

pdf("Sinusoid contingency.pdf", width=10, height=8)
ggplot(m3, aes(x=time, y=value)) + geom_line(size=1.5) + facet_grid(variable ~ name, scales="free_y") + theme_bw(base_size=24) + xlab("Time (seconds)") + ylab("") +
    theme(panel.margin = unit(20, "pt"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "grey90", colour = "grey50", size = 0.2))
dev.off()

##BMC Figure
library(R.matlab)
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "figures"))

scepticbmc <- readMat("finalicissimo_BMC_for_eLife_fig.mat")

#out.Ef contains estimated frequencies
#out.Vf contains the variance-covariance matrix of frequencies
#inside plotUncertainTimeSeries, which is called from VBA_groupBMC, it appears the SEs are derived by the sqrt of the diagonal of out.Vf

#somehow got mangled -- went into MATLAB and just saved these in a simpler .mat 
Ef <- scepticbmc$out[,,1]$Ef
scepticbmc$out[,,1]$Vf

#bmcef <- readMat("elife_bmc_frequencies.mat")
bmcef <- readMat("~/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/figures/ploscompbio_bmc_frequencies.mat")

df <- data.frame(model=unlist(bmcef$modelnames), freq=bmcef$Ef, se=sqrt(diag(bmcef$Vf)))
#df$m_ordered <- ordered(df$model, levels=c("fixed", "fixed_uv", "fixed_decay", "kalman_softmax",
#        "kalman_uv_sum", "kalman_logistic", "kalman_processnoise", "kalman_sigmavolatility", "Qstep"),
#    labels=c("Fixed LR V", "Fixed LR U + V", "Fixed LR V Decay", "KF V", "KF U + V", 
#        "KF U -> V", "KF Process Noise", "KF Volatility", "TD"))

df$m_ordered <- ordered(df$model, levels=c("fixed", "fixed_uv", "fixed_decay", "kalman_softmax",
        "kalman_uv_sum", "Qstep"),
    labels=c("Fixed LR V", "Fixed LR U + V", "Fixed LR V Sel. Maint.", "KF V", "KF U + V", "TD"))


#for ggplot with coord_flip, need to reverse
df$m_ordered <- factor(df$m_ordered, levels=rev(levels(df$m_ordered))) 

library(ggplot2)

#updated version with smaller model set for PLoS Comp Bio
pdf("SCEPTIC Main BMC v5 May2017.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.3) +
    annotate("text", x=4, y=0.65, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
        data=data.frame(x=0.65, y=0.48, label=as.character(expression(paste("BOR < ",10^{-51})))),
        hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
    ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank(), plot.margin=margin(t=5, r=10, b=5, l=0)) +
    scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()



pdf("SCEPTIC Main BMC.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity") + geom_errorbar(width=0.5) +
    annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    annotate("text", x=1.2, y=0.30, label=as.character(expression(paste("BOR = ",8.03," x ",10^{-49}))), hjust=0, vjust=0, parse=TRUE, size=6) +
    ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)))
dev.off()


library(ggplot2)
pdf("SCEPTIC Main BMC v2.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_bar(stat="identity", fill="grey92", color="black") + geom_errorbar(width=0.5) +
    annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    annotate("text", x=1.2, y=0.30, label=as.character(expression(paste("BOR = ",8.03," x ",10^{-49}))), hjust=0, vjust=0, parse=TRUE, size=6) +
    ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)))
dev.off()


library(ggplot2)
pdf("SCEPTIC Main BMC v3.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_bar(stat="identity", fill="grey92", color="black") + geom_errorbar(width=0.5) +
    annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
        data=data.frame(x=0.8, y=0.47, label=as.character(expression(paste("BOR < ",10^{-49})))),
        hjust=0, vjust=0, parse=TRUE, size=6) +
    ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)))
dev.off()

pdf("SCEPTIC Main BMC v4.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_bar(stat="identity", fill="grey92", color="black") + geom_errorbar(width=0.5) +
    annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
        data=data.frame(x=0.8, y=0.47, label=as.character(expression(paste("BOR < ",10^{-49})))),
        hjust=0, vjust=0, parse=TRUE, size=6) +
    ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())
dev.off()

pdf("SCEPTIC Main BMC v5.pdf", width=6, height=4)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
    annotate("text", x=7, y=0.54, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
        data=data.frame(x=0.8, y=0.47, label=as.character(expression(paste("BOR < ",10^{-43})))),
        hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
    ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
    scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()

#ar1 and schoenberg results
arfreqs <- readMat("ar_modelfreqs_Sep2016.mat")


#manual entry from Jon email 19Sep2016
mnames <- c("Fixed LR V",	"Fixed LR U + V", "Fixed LR V Sel. Maint.", "KF V", "KF Process Noise", "KF U + V", "KF Volatility")
df <- data.frame(model=ordered(mnames), freq=arfreqs$ar1Ef, se=sqrt(diag(arfreqs$ar1Vf)))

#for ggplot with coord_flip, need to reverse
df$m_ordered <- factor(df$model, levels=rev(levels(df$model))) 

pdf("Ar1 Main BMC v5.pdf", width=6, height=4.3)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
    annotate("text", x=5, y=0.45, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
        data=data.frame(x=0.8, y=0.40, label=as.character(expression(paste("BOR < ",10^{-32})))),
        hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
    ylab("Estimated Model Frequency") + xlab("Includes AR(1) choice") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
    scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()

#manual entry from Jon email 19Sep2016
df <- data.frame(model=ordered(mnames), freq=arfreqs$schEf, se=sqrt(diag(arfreqs$schVf)))

#for ggplot with coord_flip, need to reverse
df$m_ordered <- factor(df$model, levels=rev(levels(df$model))) 

pdf("Scho Main BMC v5.pdf", width=6, height=4.3)
ggplot(df, aes(x=m_ordered, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
    annotate("text", x=5, y=0.45, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
        data=data.frame(x=0.8, y=0.43, label=as.character(expression(paste("BOR < ",10^{-37})))),
        hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
    ylab("Estimated Model Frequency") + xlab("Includes Schoenberg choice") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
    scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75"))
dev.off()


#frank TC (replicate v5 above)
frank_bmcef <- readMat("elife_bmc_franktc_frequencies.mat")

df <- data.frame(model=unlist(frank_bmcef$models), freq=frank_bmcef$Ef, se=sqrt(diag(frank_bmcef$Vf)))

df$modelmath <- ordered(df$model, levels=c("K", "K_Lambda", "K_Lambda_Nu", "K_Lambda_Nu_AlphaG",
        "K_Lambda_Nu_AlphaG_AlphaN", "K_Lambda_Nu_AlphaG_AlphaN_Rho", "K_Lambda_Nu_AlphaG_AlphaN_Rho_Epsilon"))

#for ggplot with coord_flip, need to reverse
df$modelmath <- factor(df$modelmath, levels=rev(levels(df$modelmath))) 


pdf("Frank TC BMC v5.pdf", width=6, height=4)
ggplot(df, aes(x=modelmath, y=freq, ymin=freq-se, ymax=freq+se)) + geom_pointrange(stat="identity", size=1.3, fatten=2.5) +
    annotate("text", x=5, y=0.63, label="EP = 1.0", hjust=0, vjust=0.5, size=4.5) + # ylim(-0.05,1.1) +
    geom_label(mapping=aes(x=x,y=y,label=label, ymin=NULL, ymax=NULL), 
        data=data.frame(x=0.8, y=0.60, label=as.character(expression(paste("BOR < ",10^{-35})))),
        hjust=0, vjust=0, parse=TRUE, size=6, label.padding = unit(0.4, "lines")) +
    ylab("Estimated Model Frequency") + xlab("") + coord_flip() +
    theme_bw(base_size=20) + theme(axis.title.x=element_text(margin = margin(t = 15)),
        axis.title.y=element_text(margin = margin(r = 15)),
        panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) +
    scale_y_continuous(breaks=c(0,0.25, 0.5, 0.75), labels=c("0", ".25", ".5", ".75")) +
    scale_x_discrete("Parameter added to TC", labels=rev(expression(K, lambda, nu, alpha[G], alpha[N], rho, epsilon)))

    #scale_x_discrete("test", labels=c(expression(alpha), expression(beta)))
dev.off()
