MeanHM <- function(x,group,sigma){
    stopifnot(sigma>0)
    k <- length(x)
    if(length(group) != k)
        stop("lengths of 'x' and 'group' differ")
    fgrp <- as.factor(group)
    xbar <- tapply(x,fgrp,mean,na.rm=TRUE)
    s <- tapply(x,fgrp,sd,na.rm=TRUE)
    n <- tapply(x,fgrp,length)
    n.na <- tapply(is.na(x),fgrp,sum)
    
    y <- data.frame(
        size=n,na=n.na,
        mean=round(xbar,4),
        stdev=round(s,4))

    prop_dat <- list(N=length(xbar),sigma=sigma,x=xbar,sig=s/sqrt(n-n.na))
    fit <- rstan::sampling(stanmodels$hmmean,
                           data = prop_dat,
                           refresh=0)
    la <- rstan::extract(fit, permuted = TRUE)
    theta <- mean(la$theta); 
    mu <- apply(la$mu,2,mean); 
    y$estimate <- round(mu,4)
    out <- list(data=y,theta=theta,sigma=sigma)
    out
}


