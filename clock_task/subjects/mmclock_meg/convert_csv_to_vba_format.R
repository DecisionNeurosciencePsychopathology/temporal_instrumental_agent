csvs <- list.files(pattern="MEG.*_tc\\.csv", path=getwd())

library(dplyr)
library(tidyselect)
for (cc in csvs) {
  new <- read.csv(cc)
  new <- new %>% select(-X, -block, -null) %>%
    rename(rt=RT, score=scoreinc, probability=freq, rewFunc=function., magnitude=mag) %>%
    select(run, trial, rewFunc, emotion, magnitude, probability, score, ev, rt, starttime)
  write.csv(new, file=sub("_tc.csv", "_concat.csv", cc, fixed=TRUE), row.names=FALSE)
  
}
