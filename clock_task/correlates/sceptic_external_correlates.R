setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "correlates"))
library(gdata)
library(ggplot2)
source(file.path(getMainDir(), "Miscellaneous", "Global_Functions.R"))
#matlab params

transformedpars <- read.table("Transformed_Params.dat", sep="\t", header=TRUE)
transformedpars <- transformedpars[,!grepl("Min|Max|Mean|Std", names(transformedpars))]

library(reshape2)
m <- melt(transformedpars, id.variables="Param_and_model", variable.name="subject")
df <- dcast(m, subject ~ Param_and_model)
df$rownum <- 1:nrow(df)
df$subject <- NULL

#idlist
idlist <- read.xls("/Users/michael/Google_Drive/skinner/SCEPTIC/subject_fitting/id list.xlsx", header=FALSE, col.names="lunaid")
idlist$rownum <- 1:nrow(idlist)

df <- merge(df, idlist, by="rownum")

iq <- read.xls("MMY3_RIST_20160505.xlsx")
iq <- plyr::rename(iq, c(LunaID="lunaid"))
df <- merge(df, iq, by="lunaid")

learningvars <- c("tau_kalman_uv_sum", "beta_kalman_uv_sum", "alpha_fixed_decay", "beta_fixed_decay", "gamma_fixed_decay", 
    "prop_spread_fixed_decay", "alpha_fixed_uv", "beta_fixed_uv", "prop_spread_fixed_uv", "tau_fixed_uv")
corvars <- c("VerbalTScore", "PerformanceTScore", "RISTIndex")
corwithtarget(df, target=learningvars, withvars=corvars, pmin=.10, digits=3)

cor.test(~VerbalTScore + PerformanceTScore, df)
summary(lm(alpha_fixed_decay ~ VerbalTScore + PerformanceTScore, df))

summary(lm(alpha_fixed_decay ~ RISTIndex, df))
summary(rlm(alpha_fixed_decay ~ RISTIndex, df))
summary(lm(alpha_fixed_decay ~ RISTIndex*female, df))
summary(lm(alpha_fixed_decay ~ RISTIndex*age, df))
summary(lm(alpha_fixed_decay ~ RISTIndex*age*female, df))

summary(lm(RISTIndex ~ alpha_fixed_decay * gamma_fixed_decay, df))
cor.test(~ alpha_fixed_decay + gamma_fixed_decay, df)

df$a.c <- df$alpha_fixed_decay - mean(df$alpha_fixed_decay)
df$g.c <- df$gamma_fixed_decay - mean(df$gamma_fixed_decay)

summary(lm(RISTIndex ~ a.c * g.c, df))

pdf("IQ_Alpha_Corr.pdf", width=8, height=5)
ggplot(df, aes(x=RISTIndex, y=alpha_fixed_decay)) + geom_point() + stat_smooth(method="lm")
dev.off()

pdf("IQ_Gamma_Corr.pdf", width=8, height=5)
ggplot(df, aes(x=RISTIndex, y=gamma_fixed_decay)) + geom_point() + stat_smooth(method="lm")
dev.off()

summary(lm(gamma_fixed_decay ~ VerbalTScore + PerformanceTScore, df))

summary(lm(gamma_fixed_decay ~ RISTIndex, df))
summary(rlm(gamma_fixed_decay ~ RISTIndex, df))
summary(lm(gamma_fixed_decay ~ RISTIndex*female, df))
summary(lm(gamma_fixed_decay ~ RISTIndex*age, df))

summary(lm(beta_fixed_decay ~ RISTIndex, df))
summary(lm(beta_fixed_decay ~ RISTIndex*female, df))
summary(lm(beta_fixed_decay ~ RISTIndex*age, df))

#agelist
#agelist <- read.table("/Users/michael/Data_Analysis/clock_analysis/fmri/subinfo_db", header=TRUE)
#df <- merge(df, agelist, by="lunaid")

#selfreports (also has age and sex)
selfreports <- read.csv("/Users/michael/Tresors/DEPENd/Projects/PersonalityRest/data/surveys/SelfReportsMerged.csv")
selfreports <- plyr::rename(selfreports, c(LUNA_ID="lunaid")) 
df <- merge(df, selfreports, by="lunaid", all=TRUE)

corvars <- c("neurot", "extra", "open", "agree", "cons", "PosUrg", "NegUrg", "LackPrem", "LackPers", "SenSeek")
corwithtarget(df, target=learningvars, withvars=corvars, pmin=.05, digits=3, prewhiten=FALSE, orderbyr=FALSE)

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






