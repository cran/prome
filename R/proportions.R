PropHM <- function(x,n,kappa){
    n <- round(n)
    x <- round(x)
    j <- length(n)
    stopifnot(length(x) == j)
    stopifnot(all(x >= 0))
    stopifnot(all(n >= x))

    if(j<2) stop("< 2 groups")
    else if(j==2){
        if(missing(kappa))
            stop("'kappa' must be given")
        else
            stopifnot(kappa > 1)
    }

    if(missing(kappa)){
        if(j==2){
            out <- NULL
            warning("'kappa' must be specified when there are two groups")
        }else{
            prop_dat <- list(J = j, y = x, n = n)
            fit <- rstan::sampling(stanmodels$proportions,
                                   data = prop_dat, refresh=0,
                                   control=list(stepsize=0.01, adapt_delta=0.99))
            la <- rstan::extract(fit, permuted = TRUE)
            out1 <- apply(la$theta,2,mean); 
            out2 <- mean(la$kappa);
            out3 <- mean(la$phi); 
            y <- data.frame(size=n,event=x,
                            proportion=round(x/n,4),
                            estimate=round(out1,4))
            out <- list(data=y,alpha=out2*out3,beta=out2*(1-out3))
        }
    }else{
        prop_dat <- list(J = j, y = x, n = n, kappa=kappa)
        fit <- rstan::sampling(stanmodels$proportion,
                               data = prop_dat,
                               refresh=0)
        la <- rstan::extract(fit, permuted = TRUE)
        out1 <- apply(la$theta,2,mean);
        out2 <- kappa;
        out3 <- mean(la$phi); 
        y <- data.frame(size=n,event=x,
                        proportion=round(x/n,4),
                        estimate=round(out1,4))
        out <- list(data=y,alpha=out2*out3,beta=out2*(1-out3))
    }
    
    out
}

