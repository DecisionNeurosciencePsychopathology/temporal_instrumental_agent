#try comparing fitted parameters for TC (implemented in fitclock) versus VBA.
#ideally we would see convergence!
library(R.matlab)

source(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "vba_fmri", "clock_functions.R"))
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task"))
fitobjs <- list.files(path=file.path(getMainDir(), "clock_analysis", "fmri", "fmri_fits"), pattern="\\d+_fitinfo\\.RData", full.names=TRUE)
fitarr <- get_fit_array(fitobjs, objname="f_negeps")

pdf("epsilon_n76.pdf", width=8, height=8)
library(ggplot2)
ggplot(fitarr, aes(x=epsilonBeta)) + geom_histogram(binwidth=800) + theme_bw(base_size=24) + geom_vline(xintercept=0)
dev.off()

hist(fitarr$epsilonBeta)

prop.table(table(fitarr$epsilonBeta > 0))
mean(fitarr$epsilonBeta)
median(fitarr$epsilonBeta)

sd(fitarr$epsilonBeta)

#now get the distribution from VBA
#vba <- readMat("/Users/michael/Data_Analysis/temporal_instrumental_agent/clock_task/vba/tc_logevidence.mat")
#
#pars <- data.frame(vba$parameters[,1,]) #drop middle dimension
#names(pars) <- paste0("VBA_", strsplit(unlist(vba$models), "_")[[1]])

#this is the corrected version after fixing the RT versus reward bug
#vba <- readMat("/Users/michael/ics/temporal_instrumental_agent/clock_task/vba/tc_logevidence.mat")
vba <- readMat("/Users/michael/ics/temporal_instrumental_agent/clock_task/vba/tc_logevidence_exponentialepsilon.mat")

pars <- data.frame(vba$parameters[,7,]) #drop middle dimension
names(pars) <- paste0("VBA_", strsplit(unlist(vba$models[[7]]), "_")[[1]])
pars$lunaid <- as.vector(vba$ids)

alldf <- merge(fitarr, pars, by="lunaid")

origpars <- c("K", "lambda", "scale", "alphaG", "alphaN", "rho", "epsilonBeta")
VBApars <- c("VBA_K", "VBA_Lambda", "VBA_Nu", "VBA_AlphaG", "VBA_AlphaN", "VBA_Rho", "VBA_Epsilon")
options(width=135)
round(cor(alldf[,c(origpars, VBApars)]), 2)

corwithtarget(alldf, pmin=.01, target=origpars, withvars=VBApars)

lattice::splom(alldf[,c(origpars, VBApars)])
lattice::splom(alldf[,c("rho", "VBA_Rho", "epsilonBeta", "VBA_Epsilon")])
hist(alldf[,c("VBA_Epsilon")])
lattice::histogram(alldf[,c("epsilonBeta")], nint=15)

cor.test(alldf$VBA_Epsilon, alldf$epsilonBeta)
cor.test(alldf$VBA_Rho, alldf$rho)

#comparison of negative epsilon parameters
fitarr_negeps <- get_fit_array(fitobjs, objname="f_negeps")
fitarr_negeps <- subset(fitarr_negeps, !lunaid==10638) #bug (temporary) with missing pars here

#this was run with just one model (full), not components
vba_negeps <- readMat("/Users/michael/ics/temporal_instrumental_agent/clock_task/vba/tc_logevidence_negeps_gaussian.mat")

pars_negeps <- data.frame(vba_negeps$parameters[,1,]) #drop middle dimension
#names(pars_negeps) <- paste0("VBA_", strsplit(unlist(vba_negeps$models[[1]]), "_")[[1]])
names(pars_negeps) <- paste0("VBA_", c("K", "Sticky", "Rho", "Epsilon", "AlphaG", "AlphaN", "Decay"))
pars_negeps$lunaid <- as.vector(vba_negeps$ids)

alldf_negeps <- merge(fitarr_negeps, pars_negeps, by="lunaid")

origpars <- c("K", "stickyWeight", "stickyDecay", "alphaG", "alphaN", "rho", "epsilonBeta")
VBApars <- c("VBA_K", "VBA_Sticky", "VBA_Decay", "VBA_AlphaG", "VBA_AlphaN", "VBA_Rho", "VBA_Epsilon")
options(width=145)
round(cor(alldf_negeps[,c(origpars, VBApars)], use="pairwise.complete.obs"), 2)