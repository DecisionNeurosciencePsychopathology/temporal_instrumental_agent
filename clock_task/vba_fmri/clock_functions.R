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

