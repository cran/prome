MeanHM <- function(x,sigma){
    if(missing(sigma)){
        out <- .meanhm1(x)
    }else{
        out <- .meanhm2(x,sigma=sigma)
    }
    out
}

.meanhm1 <- function(x){
    x <- as.matrix(x)
    nr <- nrow(x)
    nc <- ncol(x)
    mydat <- list(N=nr,M=nc,x=t(x))
    fit <- rstan::sampling(stanmodels$hm2mean,
                               data = mydat,
                               refresh=0)
    la <- rstan::extract(fit, permuted = TRUE)
    theta.mean <- mean(la$theta); 
    theta.median <- median(la$theta); 
    sigma.median <- sqrt(median(la$sigma)); 
    sigma.mean <- sqrt(mean(la$sigma)); 
    out <- list(theta.mean=theta.mean,
                theta.median=theta.median,
                sigma.mean=sigma.mean,
                sigma.median=sigma.median)
}

.meanhm2 <- function(x,sigma){
    x <- as.matrix(x)
    nr <- nrow(x)
    nc <- ncol(x)
    if(sigma <= 0) stop("invalid 'sigma'")
    mydat <- list(N=nr,M=nc,x=t(x),sigma=sigma)
    fit <- rstan::sampling(stanmodels$hmmean,
                               data = mydat,
                               refresh=0)
    la <- rstan::extract(fit, permuted = TRUE)
    theta.mean <- mean(la$theta); 
    theta.median <- median(la$theta); 
    out <- list(theta.mean=theta.mean,
                theta.median=theta.median)
}


