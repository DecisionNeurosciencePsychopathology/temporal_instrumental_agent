setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "correlates"))
library(gdata)
library(ggplot2)
library(reshape2)
library(R.matlab)
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
#matlab params

transformedpars <- read.table("Transformed_Params.dat", sep="\t", header=TRUE)
transformedpars <- transformedpars[,!grepl("Min|Max|Mean|Std", names(transformedpars))]

#from matlab
#this is an updated version where the prop spread should be fixed to .0125, per testing
nontransformedpars <- read.csv("sceptic_pars_nontransformed.csv") #haha, for some reason, these are definitely *Transformed* into their native scaling...
nontransformedpars$rownum <- 1:nrow(nontransformedpars)

#and the original raw info
fromfmri <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba_fmri/posterior_states_decay_nomultisession_psfixed0p0125.mat")
parmat <- data.frame(fromfmri$pars)[,1:5] #6th column is blank under fixed prop spread
names(parmat) <- c("lunaid", "fmri_F", "fmri_alpha", "fmri_gamma", "fmri_beta")
#funny glitch where the id is coming through with "1010" at the end due to digits in filename during parsing.
parmat$lunaid <- as.numeric(substr(as.character(parmat$lunaid), 1,5)) #only retain first five digits
parmat$fmri_alpha_t <- 1/(1+exp(-parmat$fmri_alpha)) #logistic
parmat$fmri_gamma_t <- 1/(1+exp(-parmat$fmri_gamma)) #logistic
parmat$fmri_beta_t <- exp(parmat$fmri_beta) #exponential

#library(reshape2)
m <- melt(transformedpars, id.variables="Param_and_model", variable.name="subject")
df <- dcast(m, subject ~ Param_and_model)
df$rownum <- 1:nrow(df)
df$subject <- NULL

#idlist
idlist <- read.xls("/Users/michael/Google_Drive/skinner/projects_analyses/SCEPTIC/subject_fitting/id list.xlsx", header=FALSE, col.names="lunaid")
idlist$rownum <- 1:nrow(idlist)

#df <- merge(df, idlist, by="rownum")
df <- merge(nontransformedpars, idlist, by="rownum")

df <- merge(df, parmat, by="lunaid")

#check fmri fit pars versus pars from Jon
#Good news: perfect convergence
cor(df[,c("fixed_decay_beta", "fixed_decay_gamma", "fixed_decay_alpha", "fmri_alpha_t", "fmri_gamma_t", "fmri_beta_t")])

iq <- read.xls("MMY3_RIST_20160505.xlsx")
iq <- plyr::rename(iq, c(LunaID="lunaid"))
df <- merge(df, iq, by="lunaid")

#selfreports (also has age and sex)
selfreports <- read.csv("/Users/michael/Tresors/DEPENd/Projects/SPECC/SelfReports/data/surveys/SelfReportsMerged.csv")
selfreports <- plyr::rename(selfreports, c(LUNA_ID="lunaid"))

df <- merge(df, selfreports, by="lunaid", all=TRUE)
df <- subset(df, !lunaid==11246) #huge head movement -- exclude?

#learningvars <- c("tau_kalman_uv_sum", "beta_kalman_uv_sum", "alpha_fixed_decay", "beta_fixed_decay", "gamma_fixed_decay", 
#    "prop_spread_fixed_decay", "alpha_fixed_uv", "beta_fixed_uv", "prop_spread_fixed_uv", "tau_fixed_uv")

learningvars <- c("kalman_uv_sum_tau", "kalman_uv_sum_beta", "fixed_decay_alpha", "fixed_decay_beta", "fixed_decay_gamma", 
    "fixed_uv_alpha", "fixed_uv_beta", "fixed_uv_tau", "fmri_F", "fmri_alpha", "fmri_gamma", "fmri_beta")

corvars <- c("VerbalTScore", "PerformanceTScore", "RISTIndex")
corwithtarget(df, target=learningvars, withvars=corvars, pmin=.10, digits=3)

#gamma correlation is a bit stronger under transformed variant than Gaussian variant
#look at plots
lattice::splom(df[,c(corvars, "fmri_gamma", "fixed_decay_gamma", "fmri_alpha", "fixed_decay_alpha")])

cor.test(~VerbalTScore + PerformanceTScore, df)
summary(lm(alpha_fixed_decay ~ VerbalTScore + PerformanceTScore, df))

summary(lm(alpha_fixed_decay ~ RISTIndex, df))
summary(rlm(alpha_fixed_decay ~ RISTIndex, df))
summary(lm(alpha_fixed_decay ~ RISTIndex*female, df))
summary(lm(alpha_fixed_decay ~ RISTIndex*age, df))
summary(lm(alpha_fixed_decay ~ RISTIndex*age*female, df))

summary(mfull <- lm(RISTIndex ~ fixed_decay_alpha * fixed_decay_gamma, df))

badsubjs <- c()

#there is some funny shift in the values (very slight) when self reports are merged above. Why??
summary(mfull <- lm(RISTIndex ~ fmri_alpha * fmri_gamma, df))
summary(mdecay <- lm(RISTIndex ~ fmri_gamma, df))
summary(mage <- lm(RISTIndex ~ fmri_gamma * age, df))
summary(msex <- lm(RISTIndex ~ fmri_gamma * female, df))
anova(mdecay, mfull, mage, msex)

car::Anova(mfull)
#no evidence of moderation by sex or age
summary(m2<- lm(RISTIndex ~ fmri_alpha * fmri_gamma * age * female, df))
summary(m2<- lm(RISTIndex ~ fmri_gamma * age * female, df))
summary(m3<- lm(RISTIndex ~ fmri_alpha * fmri_gamma + fmri_alpha*age + fmri_gamma*age + fmri_alpha*female + fmri_gamma*female, df))
anova(mfull, m3, m2)
car::Anova(m2)


#look at models under Gaussian parameters
summary(mfull <- lm(RISTIndex ~ fmri_alpha * fmri_gamma, df))
summary(mfull <- rlm(RISTIndex ~ fmri_alpha * fmri_gamma, df))
covRob(df[,c("RISTIndex", "fmri_alpha", "fmri_gamma")], corr=TRUE)
covRob(df[,c("RISTIndex", "PerformanceTScore", "fixed_decay_alpha", "fixed_decay_gamma")], corr=TRUE)
round(cor(df[,c("RISTIndex", "fmri_alpha", "fmri_gamma")]), 3)

#transformed par
ggplot(df, aes(x=fixed_decay_gamma, y=RISTIndex)) + geom_point() + stat_smooth(method="rlm")

#gaussian counterpart
#gaussian tends to overemphasize influence of low gamma values, which when transformed, are trivially difference
#conclusion: stick with transformed variant
ggplot(df, aes(x=fmri_gamma, y=RISTIndex)) + geom_point() + stat_smooth(method="rlm")

ggplot(df, aes(x=fixed_decay_gamma, y=age)) + geom_point() + stat_smooth()
ggplot(df, aes(x=fixed_decay_alpha, y=age)) + geom_point() + stat_smooth()
ggplot(df, aes(x=fmri_alpha, y=age)) + geom_point() + stat_smooth()

s <- stepAIC(mfull, direction="both")
summary(lm(RISTIndex ~ fixed_decay_alpha:fixed_decay_gamma, df))
summary(m1 <- lm(RISTIndex ~ fixed_decay_alpha, df))
summary(m2 <- lm(RISTIndex ~ fixed_decay_gamma, df))
summary(m3 <- lm(RISTIndex ~ fixed_decay_alpha:fixed_decay_gamma, df))
cor.test(~ fixed_decay_alpha + fixed_decay_gamma, df)
library(leaps)
ll<-regsubsets(RISTIndex ~ fixed_decay_alpha * fixed_decay_gamma,data=df,nbest=10)

df$a.c <- df$alpha_fixed_decay - mean(df$alpha_fixed_decay)
df$g.c <- df$gamma_fixed_decay - mean(df$gamma_fixed_decay)

summary(lm(RISTIndex ~ a.c * g.c, df))

cor.test(df$RISTIndex, df$fixed_decay_gamma)
cor.test(df$PerformanceTScore, df$fixed_decay_gamma)
cor.test(df$PerformanceTScore, df$fixed_decay_alpha)
cor.test(df$PerformanceTScore, df$fmri_alpha)
cor.test(df$VerbalTScore, df$fixed_decay_gamma)

pdf("IQ_Alpha_Corr.pdf", width=8, height=5)
ggplot(df, aes(x=RISTIndex, y=fixed_decay_alpha)) + geom_point() + stat_smooth(method="lm")
dev.off()

pdf("IQ_Alpha_Corr.pdf", width=8, height=5)
ggplot(df, aes(x=RISTIndex, y=alpha_fixed_decay)) + geom_point() + stat_smooth(method="lm")
dev.off()

pdf("IQ_Gamma_Corr.pdf", width=8, height=5)
ggplot(df, aes(x=RISTIndex, y=fixed_decay_gamma)) + geom_point() + stat_smooth(method="lm")
dev.off()

df$performanceIQ <- 100 + (df$PerformanceTScore - 50)/10*15 #rescale T -> IQ
pdf("PerformanceIQ_Gamma_Corr.pdf", width=5, height=3)
ggplot(df, aes(x=performanceIQ, y=fixed_decay_gamma)) + geom_point() + stat_smooth(method="lm", se=FALSE) +
    theme_bw(base_size=15) + ylab(expression(paste("Decay parameter (",gamma, ")"))) + xlab("RIST Nonverbal Intelligence") +
    theme(axis.title.y=element_text(margin=margin(0,10,0,0))) +
    theme(axis.title.x=element_text(margin=margin(10,0,0,0)))
dev.off()

#pdf("IQ_Gamma_Corr.pdf", width=8, height=5)
#ggplot(df, aes(x=RISTIndex, y=gamma_fixed_decay)) + geom_point() + stat_smooth(method="lm")
#dev.off()

summary(lm(gamma_fixed_decay ~ VerbalTScore + PerformanceTScore, df))

summary(lm(gamma_fixed_decay ~ RISTIndex, df))
summary(rlm(gamma_fixed_decay ~ RISTIndex, df))
summary(lm(gamma_fixed_decay ~ RISTIndex*female, df))
summary(lm(gamma_fixed_decay ~ RISTIndex*age, df))

summary(lm(beta_fixed_decay ~ RISTIndex, df))
summary(lm(beta_fixed_decay ~ RISTIndex*female, df))
summary(lm(beta_fixed_decay ~ RISTIndex*age, df))

corvars <- c("neurot", "extra", "open", "agree", "cons", "PosUrg", "NegUrg", "LackPrem", "LackPers", "SenSeek")
corwithtarget(df, target=learningvars, withvars=corvars, pmin=.05, digits=3, prewhiten=FALSE, orderbyr=FALSE)

plot(df$fixed_decay_alpha, df$LackPers)

#summary statistics
df <- merge(df, allfits, by="lunaid")

cor.test(df$alpha_fixed_decay, df$age)
plot(df$age, df$alpha_fixed_decay)
cor.test(df$beta_fixed_decay, df$age)
cor.test(df$gamma_fixed_decay, df$age)

cor.test(df$gamma_fixed_decay, df$female)
cor.test(df$tau_fixed_uv, df$age)

df$female <- factor(df$female, levels=c(0,1), labels=c("male", "female")) 

summary(lm(age ~ alpha_fixed_decay + beta_fixed_decay + gamma_fixed_decay, df))

summary(lm(age ~ alpha_fixed_decay + beta_fixed_decay + gamma_fixed_decay, df))

cor.test(df$extra, df$totreward)
cor.test(df$extra, df$avg_ev)
cor.test(df$extra, df$evhappy)


cor.test(df$extra, df$evhappy)


cor.test(df$tau_kalman_uv_sum, df$avg_ev)
cor.test(df$tau_fixed_uv, df$avg_ev)
plot(df$tau_fixed_uv, df$avg_ev)
cor.test(df$extra, df$tau_fixed_uv)






