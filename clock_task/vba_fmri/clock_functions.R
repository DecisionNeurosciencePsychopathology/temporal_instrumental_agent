getClockGroupData <- function(path, pattern=".*tcExport.csv", idpattern="^.*/fMRIEmoClock_(\\d+)_tc_tcExport\\.csv$") {
  require(plyr)
  tcFiles <- list.files(path=path, pattern=pattern, full.names=TRUE)
  
  allData <- list()
  for (f in tcFiles) {
    subject <- sub(idpattern, "\\1", f, perl=TRUE)
    sdata <- read.csv(f, header=TRUE, comment.char="#")
    sdata$Null <- NULL #delete dummy column
    sdata$Subject <- factor(subject)
    sdata <- ddply(sdata, .(emotion, rewFunc), function(subdf) {
          bestRew <- -1
          bestEV <- -1
          bestRewRT <- subdf[1, "rt"]
          bestEVRT <- subdf[1, "rt"]
          subdf$bestRewRT <- NA_real_
          subdf$bestEVRT <- NA_real_
          for (i in 1:nrow(subdf)) {
            if (subdf[i,"score"] >= bestRew) {
              bestRewRT <- subdf[i,"bestRewRT"] <- subdf[i,"rt"]
              bestRew <- subdf[i,"score"]
            } else {
              subdf[i,"bestRewRT"] <- bestRewRT
            }
            if (subdf[i,"score"] > 0 && subdf[i,"ev"] >= bestEV) {
              bestEVRT <- subdf[i,"bestEVRT"] <- subdf[i,"rt"]
              bestEV <- subdf[i,"ev"]
            } else {
              subdf[i,"bestEVRT"] <- bestEVRT
            }
          }
          subdf
        })
    allData[[f]] <- sdata
  }
  
  ##rearrange column headers for readability
  allData <- do.call(rbind, allData)
  row.names(allData) <- NULL
  allData <- allData[,c("Subject", "run", "trial", "rewFunc", "emotion", "magnitude", "probability", "score", "ev", "rt", "bestRewRT", "bestEVRT", "image")]
  allData <- plyr::rename(allData, c(Subject="LunaID"))
  #allData$trialRel <- #unlist(lapply(split(allData, f=list(allData$LunaID, allData$rewFunc, allData$emotion)), function(l) { return(1:nrow(l)) } ))
  allData <- ddply(allData, .(LunaID, rewFunc, emotion), function(subdf) {
        subdf$trial_abs <- subdf$trial
        subdf$trial <- 1:nrow(subdf)
        return(subdf)
      })
  return(allData)
}

#function to extract id, run condition, and outcome statistics
get_fit_array <- function(fitobjs) {
  allfits <- c()
  for (i in 1:length(fitobjs)) {
    cat("Loading object: ", fitobjs[i], "\n")
    loc <- local({load(fitobjs[i]); environment()})$f_poseps #time-clock fit object
    emo <- loc$run_condition
    rew <- loc$rew_function
    pars <- loc$theta[,"cur_value"]
    
    rw <- local({load(fitobjs[i]); environment()})$f_value #rescorla-wagner fit object
    rewhappy <- sum(rw$Reward[which(rw$run_condition == "happy"),])
    rewfear <- sum(rw$Reward[which(rw$run_condition == "fear"),])
    rewscram <- sum(rw$Reward[which(rw$run_condition == "scram"),])
    
    evhappy <- sum(rw$ev[which(rw$run_condition == "happy"),])
    evfear <- sum(rw$ev[which(rw$run_condition == "fear"),])
    evscram <- sum(rw$ev[which(rw$run_condition == "scram"),])
    
    pars <- c(pars, rw_alphaV=rw$theta["alphaV", "cur"], rw_betaV=rw$theta["betaV", "cur"], avg_ev=mean(rw$ev), totreward=sum(rw$Reward),
        rewhappy=rewhappy, rewfear=rewfear, rewscram=rewscram,
        evhappy=evhappy, evfear=evfear, evscram=evscram)
    
    lunaid <- as.integer(sub("[^\\d]*(\\d+)_fitinfo.RData$", "\\1", fitobjs[i], perl=TRUE))
    bpdsub <- lunaid < 10000
    allfits <- rbind(allfits, c(lunaid=lunaid, bpd=bpdsub, pars))
  }
  
  allfits <- data.frame(allfits)
}