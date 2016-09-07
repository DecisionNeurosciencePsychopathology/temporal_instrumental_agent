library(R.matlab)
library(lattice)
library(reshape2)
#matdir <- "/Users/michael/ics/temporal_instrumental_agent/clock_task/optimality_testing/output"
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "optimality_testing"))
#completed <- list.files(matdir, pattern="optimize_output.*\\.mat", full.names=TRUE)
completed <- list.files("output", pattern="optimize_output.*\\.mat", full.names=TRUE)
contorder <- c("IEV", "DEV", "QUADUP", "IEVLINPROB", "DEVLINPROB", "ALL") #order of contingencies fit in MATLAB optmat

allFits <- list()
for (matname in completed) {
  m <- readMat(matname)
  modelname <- as.vector(m$agents[["name",1,1]])
#  parnames <- switch(modelname,
#      fixedLR_softmax = c("prop_spread", "beta", "alpha"),
#      fixedLR_egreedy = c("prop_spread", "epsilon", "alpha"),
#      fixedLR_egreedy_grw = c("prop_spread", "epsilon", "alpha", "sig_grw"),
#      asymfixedLR_softmax = c("prop_spread", "beta", "alpha", "rho"),
#      kalman_softmax = c("prop_spread", "beta"),
#      kalman_processnoise = c("prop_spread", "beta", "omega"),
#      kalman_sigmavolatility = c("prop_spread", "beta", "phi", "gamma"),
#      kalman_uv_logistic = c("prop_spread", "tradeoff", "discrim"),
#      kalman_uv_sum = c("prop_spread", "beta", "tau"),
#      fixedLR_kl_softmax = c("prop_spread", "beta", "alpha", "kappa", "lambda"),
#      kalman_kl_softmax = c("prop_spread", "beta", "kappa", "lambda"),
#      kalman_processnoise_kl = c("prop_spread", "beta", "omega", "kappa", "lambda"),
#      kalman_uv_sum_kl = c("prop_spread", "beta", "tau", "kappa", "lambda"),
#      franktc = c('lambda', 'epsilon', 'alphaG', 'alphaN', 'K', 'nu', 'rho')
#  )
  
  parnames <- unlist(m$agents[["parnames",1,1]]) #in case a parameters (e.g., beta) was fixed, switch above is cumbersome
  
  allcosts <- drop(m$costs)
  #dimnames(allcosts) <- list(rep=1:nrow(allcosts), contingency=contorder)
  names(allcosts) <- list(rep=1:length(allcosts), contingency="sinusoid")
  
  pars <- data.frame(do.call(rbind, lapply(m$pars, function(el) { unlist(el)})))
  names(pars) <- parnames
  #pars$cont <- rep(contorder, each=100)
  pars$cont <- "sinusoid"
  
  allFits[[modelname]] <- list(modelname = modelname, costs=allcosts, pars=pars)
  
}

pdf("Par_histograms_sinusoid.pdf", width=10, height=10)
lapply(allFits, function(model) {
      mpars <- melt(model$pars, id.vars="cont")
      plot(histogram(~ value | variable + cont, mpars, 
              main=model$modelname, scales=list(x=list(relation="free")), breaks=NULL))
    })
dev.off()

combCosts <- do.call(rbind, lapply(allFits, function(model) {
      d <- data.frame(model$costs)
      d$model <- factor(model$modelname)
      d
    }))

combCosts$model.costs <- -1*combCosts$model.costs #more is better

m <- melt(combCosts, id.vars="model")
library(ggplot2)
pdf("Optimal costs_sinusoid_stackedbars.pdf", width=15, height=10)
ggplot(m, aes(x=value, fill=model)) + facet_wrap(~variable, scales="free") + geom_histogram() + scale_fill_hue("Model") + theme_bw(base_size=16) #scale_fill_brewer("Model", palette="Set3")
dev.off()

pdf("Optimal costs_sinusoid.pdf", width=15, height=10)
ggplot(combCosts, aes(x=model, y=model.costs)) + geom_boxplot() + theme_bw(base_size=16) + theme(axis.text.x=element_text(angle=90)) + ylab("Total points") #scale_fill_brewer("Model", palette="Set3")
dev.off()


#look at mean and median costs across optimizations
#definitely suggests preference for UV sum and UV logistic models
sort(tapply(combCosts$model.costs, combCosts$model, mean))
sort(tapply(combCosts$model.costs, combCosts$model, median))

#I believe this prints the median parameters
lapply(allFits, function(model) {
      mpars <- melt(model$pars, id.vars="cont")
      tapply(mpars$value, list(mpars$cont, mpars$variable), median)
    })

summary(mm <- aov(model.costs ~ model, combCosts))
tuk <- TukeyHSD(mm)
plot(tuk)
tuk


## All EV Equated Contingencies (Aug 2016)
completed <- list.files("output", pattern="optimize_output_allequate.*\\.mat", full.names=TRUE)

allFits <- list()
for (matname in completed) {
  m <- readMat(matname)
  modelname <- as.vector(m$agents[["name",1,1]])
  
  parnames <- unlist(m$agents[["parnames",1,1]]) #in case a parameters (e.g., beta) was fixed, switch above is cumbersome
  
  allcosts <- drop(m$costs)
  #dimnames(allcosts) <- list(rep=1:nrow(allcosts), contingency=contorder)
  #names(allcosts) <- list(rep=1:length(allcosts), contingency="allmono")
  
  pars <- data.frame(do.call(rbind, lapply(m$pars, function(el) { unlist(el)})))
  names(pars) <- parnames
  #pars$cont <- rep(contorder, each=100)
  pars$cont <- "allmono"
  
  allFits[[modelname]] <- list(modelname = modelname, costs=allcosts, pars=pars)
  
}

pdf("Par_histograms_allequated.pdf", width=10, height=10)
lapply(allFits, function(model) {
      mpars <- melt(model$pars, id.vars="cont")
      plot(histogram(~ value | variable + cont, mpars, 
              main=model$modelname, scales=list(x=list(relation="free")), breaks=NULL))
    })
dev.off()

combCosts <- do.call(rbind, lapply(allFits, function(model) {
          d <- data.frame(model$costs)
          d$model <- factor(model$modelname)
          d
        }))

combCosts$model.costs <- -1*combCosts$model.costs #more is better

m <- melt(combCosts, id.vars="model")
library(ggplot2)
pdf("Optimal costs_allequate_stackedbars.pdf", width=15, height=10)
ggplot(m, aes(x=value, fill=model)) + facet_wrap(~variable, scales="free") + geom_histogram() + scale_fill_hue("Model") + theme_bw(base_size=16) #scale_fill_brewer("Model", palette="Set3")
dev.off()

pdf("Optimal costs_allequated.pdf", width=15, height=10)
ggplot(combCosts, aes(x=model, y=model.costs)) + geom_boxplot() + theme_bw(base_size=16) + theme(axis.text.x=element_text(angle=90)) + ylab("Total points") #scale_fill_brewer("Model", palette="Set3")
dev.off()
