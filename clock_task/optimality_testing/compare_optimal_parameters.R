library(R.matlab)
library(lattice)
library(reshape2)
matdir <- "/Users/michael/ics/temporal_instrumental_agent/clock_task/optimality_testing/output"
setwd(file.path(getMainDir(), "temporal_instrumental_agent", "clock_task", "optimality_testing"))
completed <- list.files(matdir, pattern="optimize_output.*\\.mat", full.names=TRUE)
contorder <- c("IEV", "DEV", "QUADUP", "IEVLINPROB", "DEVLINPROB", "ALL") #order of contingencies fit in MATLAB optmat

allFits <- list()
for (matname in completed) {
  m <- readMat(matname)
  modelname <- as.vector(m$agents[["name",1,1]])
  parnames <- switch(modelname,
      fixedLR_softmax = c("prop_spread", "beta", "alpha"),
      fixedLR_egreedy = c("prop_spread", "epsilon", "alpha"),
      fixedLR_egreedy_grw = c("prop_spread", "epsilon", "alpha", "sig_grw"),
      asymfixedLR_softmax = c("prop_spread", "beta", "alpha", "rho"),
      kalman_softmax = c("prop_spread", "beta"),
      kalman_processnoise = c("prop_spread", "beta", "omega"),
      kalman_sigmavolatility = c("prop_spread", "beta", "phi", "gamma"),
      kalman_uv_logistic = c("prop_spread", "tradeoff", "discrim"),
      kalman_uv_sum = c("prop_spread", "beta", "tau"),
      fixedLR_kl_softmax = c("prop_spread", "beta", "alpha", "kappa", "lambda"),
      kalman_kl_softmax = c("prop_spread", "beta", "kappa", "lambda"),
      kalman_processnoise_kl = c("prop_spread", "beta", "omega", "kappa", "lambda"),
      kalman_uv_sum_kl = c("prop_spread", "beta", "tau", "kappa", "lambda"),
      franktc = c('lambda', 'epsilon', 'alphaG', 'alphaN', 'K', 'nu', 'rho')
  )
  
  allcosts <- drop(m$costs)
  dimnames(allcosts) <- list(rep=1:nrow(allcosts), contingency=contorder)
  
  pars <- data.frame(do.call(rbind, lapply(m$pars, function(el) { unlist(el)})))
  names(pars) <- parnames
  pars$cont <- rep(contorder, each=100)
  
  allFits[[modelname]] <- list(modelname = modelname, costs=allcosts, pars=pars)
  
}

pdf("Par_histograms.pdf", width=10, height=10)
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

m <- melt(combCosts, id.vars="model")
library(ggplot2)
pdf("Optimal costs.pdf", width=15, height=10)
ggplot(m, aes(x=value, fill=model)) + facet_wrap(~variable, scales="free") + geom_histogram() + scale_fill_brewer("Model", palette="Set3") + theme_bw(base_size=16)
dev.off()

lapply(allFits, function(model) {
      mpars <- melt(model$pars, id.vars="cont")
      tapply(mpars$value, list(mpars$cont, mpars$variable), median)
    })