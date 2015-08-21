##simulate positive and negative correlation data for rbffit
library(mvtnorm)
#correlated data at r=0.9
xy <- rmvnorm(200, mean=c(0,0), sigma=rbind(c(1,0.9),c(0.9,1)))
cor(xy)
plot(xy)
write.table(xy, file="xy_cor0.9.txt", row.names=FALSE, col.names=FALSE, sep="\t")
