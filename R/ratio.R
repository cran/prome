
## dratio(1,mu=c(1,2),sigma=c(2.5,3.5))
## pratio(0.30,mu=c(1.1,6.5),sigma=c(1.5,2))

## use numerical method to compute the mean and SD
pratio <- function(x,mu,sigma, conf.level=0.95){
    .dratio <- function(x,mu,sigma){
        v <- sigma^2
        A <- sqrt(x^2/v[1]+1/v[2])
        B <- mu[1]/v[1]*x + mu[2]/v[2]
        C <- sum(mu^2/v)
        D <- exp(0.5*((B/A)^2-C))
        dx <- B*D/A^3/sqrt(2*pi)/prod(sigma) *
            (pnorm(B/A,0,1)-pnorm(-B/A,0,1)) +
            1/A^2/pi/prod(sigma)*exp(-C/2)
        return(dx)
    }
    stopifnot(conf.level>0&&conf.level<1)
    stopifnot(length(mu)==2)
    stopifnot(length(sigma)==2)
    stopifnot(all(sigma>0))
    alp2 <- 0.5 - 0.5*conf.level
    xbar <- mu[1]/mu[2]
    s <- abs(xbar)*sqrt(sum(sigma^2/mu^2))
    x0 <- seq(xbar-4*s, xbar+4*s, length=4001)
    d <- diff(x0)[1]
    f0 <- dratio(x0,mu,sigma)
    print(c(d,sum(f0*d)))
    mu <- sum(f0*d*x0)
    F0 <- cumsum(f0*d)
    F0 <- F0/max(F0)*(1-min(F0))
    y <- abs(F0 - 0.5)
    i <- which(y==min(y))[1]
    med <- x0[i]
    y <- abs(F0 - alp2)
    i <- which(y==min(y))[1]
    ll <- x0[i]
    y <- abs(F0 - 1 + alp2)
    i <- which(y==min(y))[1]
    ul <- x0[i]
    y <- abs(x0 - x)
    i <- which(y==min(y))[1]
    Fx <- 1-F0[i]
    list(x=x0,y=f0,mean=mu,median=med,prob=Fx, ci=c(ll,ul))
}
    
