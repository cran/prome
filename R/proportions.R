PropHM <- function(x,n,weights){
    n <- round(n)
    x <- round(x)
    j <- length(n)
    if(j<2) stop("< 2 groups")
    stopifnot(length(x) == j)
    stopifnot(all(x >= 0))
    stopifnot(all(n >= x))
    if(missing(weights)) weights <- rep(1,j)
    else{
        stopifnot(length(weights) == j)
        stopifnot(all(weights > 0 & weights <= 1))
    }
    
    prop_dat <- list(J = j, y = x, n = n, pwr=weights)
    fit <- rstan::sampling(stanmodels$proportion,
                           data = prop_dat,
                           refresh=0)
    la <- rstan::extract(fit, permuted = TRUE)
    out1 <- apply(la$theta,2,mean); 
    out2 <- mean(la$alpha);
    y <- data.frame(size=n,event=x,
                    weights=weights,
                    proportion=round(x/n,4),
                    estimate=round(out1,4))
    out <- list(data=y,alpha=out2,beta=2-out2)
    out
}

