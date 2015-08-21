t <- 0:5000

#pure quad
yt <- function(t) {
    mag <- 0.00002*(t-1800)^2+20
    freq <- 1-.62*t/5000
    return(list(mag=mag, freq=freq, ev=mag*freq))
}

#yt <- function(t) {
#    return(.0001*(t-1500)^3 +  0.00001*(t-2500)^2+50)
#}

func <- yt(t)
#range(yt(t)[["mag"]])

png("quadup.png", width=8, height=6, units="in", res=200)
par(mfrow=(c(3,1)))
par(mar=c(4,4,2,3))
plot(t, yt(t)[["mag"]], type="l", xlab="timestep", ylab="magnitude")
plot(t, yt(t)[["freq"]], type="l", xlab="timestep", ylab="frequency")
plot(t, yt(t)[["ev"]], type="l", xlab="timestep", ylab="ev")
dev.off()
