
#general form of logistic function (including variation of base of exponent)
#c sets the upper bound of the function (1 is a good default for us)
#b controls
logistic <- function(x, a, b, c=1, base=exp(1)) {
    #c is max (asymptote) of function
    return(c/(1+a*base^-(x*b)))
    #return(c/(a*base^(-(x*b))))
}

x <- -15:250

#u <- rev(seq(-.25, 0, by=.001))
#epsilon <- .15
#x <- epsilon*(1/(-u + .001))
plot(1:length(x), x, type="l")

plot(x, logistic(x,a=1,b=-.5, base=exp(1)), type="l")
plot(x, logistic(x,a=100000,b=0.5), type="l")

#two parameter logistic model
twoPL <- function(x, A, B) {
    return(1/(1+exp(-A*(x-B))))
}



x <- seq(0,-0.25, by=-.01)
plot(x, twoPL(x,A=100,B=-.13), type="l")
lines(x, twoPL(x,A=100,B=-.05), type="l", lty=5)

indifferencepoint <- -.13
library(animation)
#ani.options(interval=1/5)
saveVideo({
    for (i in 1:length(x)) {
        plot(x, twoPL(x,A=100,B=indifferencepoint), type="l", xlab="uncertainty", ylab="p(explore)")
        points(x[i], twoPL(x[i], A=100,B=indifferencepoint))
        text(x=-.20, y=1, labels="Indiff (b) = -.13")
    }
    for (i in 1:length(x)) {
        plot(x, twoPL(x,A=100,B=indifferencepoint), type="l", xlab="uncertainty", ylab="p(explore)")
        text(x=-.20, y=1, labels="Indiff (b) = -.13")
        
        lines(x, twoPL(x,A=100,B=-.06), type="l", xlab="uncertainty", ylab="p(explore)", lty=5)
        points(x[i], twoPL(x[i], A=100,B=-.06))
        text(x=-.20, y=0.9, labels="Indiff (b) = -.06")
    }
}
          , video.name="uncertaintyWalk.mp4", interval=1/5)


rasch <- function(x, difficulty, discrimination) {
    return(1/(1+exp(-(x-difficulty))))
}

#x <- 400:600

uncertainty <- seq(-10, 10, by=.01)

plot(uncertainty, rasch(uncertainty, difficulty = indifferencepoint), type="l")

#a is discrimination (steepness)
#b is difficulty (location)
#c is guessing (lower asymptote of function)

threePL <- function(x, a=1, b, c) {
    return(c + (1-c)/(1+exp(-a*(x - b))))
}

