memixed <- function(x0,x1,group,method="ATT"){
    n <- length(group)
    fgrp <- as.factor(group)
    ngrp <- as.numeric(fgrp)-1
    if(nlevels(fgrp) != 2)
        stop("current version support two groups only")
    x0 <- as.matrix(x0)
    x1 <- as.matrix(x1)
    if(nrow(x0) != n)
        stop("different numbers of subjects in 'x0' and 'group'")
    if(nrow(x1) != n)
        stop("different numbers of subjects in 'x1' and 'group'")

    type <- match.arg(tolower(method),c("att","ate"))
    if(type=="att"){
        sampledat <- list(N=n, group = ngrp,
                          r0 = ncol(x0), r1 = ncol(x1),
                          w0 = t(x0), w1=t(x1))
        fit3 <- rstan::sampling(stanmodels$att,
                                data = sampledat,
                                iter=5000*2,
                                refresh=0)
        la <- rstan::extract(fit3, permuted = TRUE)
        es0 <- mean(la$es0); 
        es1 <- mean(la$es1); 
        mu0 <- mean(la$mu0) 
        mu1 <- mean(la$mu1) 
        sig0 <- mean(la$sigma0) 
        sig1 <- mean(la$sigma1) 
        tau <- mean(la$tau)
        s1 <- mean(la$sig_active)
        s0 <- mean(la$sig_sham)
        out <- structure(
            list(xfit=fit3,mean0.t0=mu0,mean1.t0=mu1,
                 sig0.t0=sig0,sig1.t0=sig1,
                 sig.me=tau,mu.active=es1,sig.active=s1,
                 mu.sham=es0,sig.sham=s0,groups=levels(fgrp)),
            class='memix')
    }else if(type=="ate"){
        sampledat <- list(N=n, group = ngrp,
                          r0 = ncol(x0), r1 = ncol(x1),
                          w0 = t(x0), w1=t(x1))
        fit3 <- rstan::sampling(stanmodels$ate,
                                data = sampledat,
                                iter=5000*2,
                                refresh=0)
        la <- rstan::extract(fit3, permuted = TRUE)
        es0 <- mean(la$es0); 
        es1 <- mean(la$es1); 
        mu0 <- mean(la$mu0) 
        mu1 <- mean(la$mu1) 
        sig0 <- mean(la$sigma0) 
        sig1 <- mean(la$sigma1) 
        tau <- mean(la$tau)
        s1 <- mean(la$sig_active)
        s0 <- mean(la$sig_sham)
        out <- structure(
            list(xfit=fit3,mean0.t0=mu0,mean1.t0=mu1,
                 sig0.t0=sig0,sig1.t0=sig1,
                 sig.me=tau,mu.active=es1,sig.active=s1,
                 mu.sham=es0,sig.sham=s0,groups=levels(fgrp)),
            class='memix')
    }else{
        out <- NULL
        warning("method not supported yet")
    }
    out
}

print.memix <- function(x,...){
    if(class(x) != "memix"){
        print(x)
    }else{
        out <- data.frame(mean=c(x$mean0.t0,x$mean1.t0,x$mu.active,x$mu.sham,0),
                          sd=c(x$sig0.t0,x$sig1.t0,x$sig.active,x$sig.sham,x$sig.me))
        nam1 <- paste(x$groups,".baseline",sep='')
        rownames(out) <- c(nam1,"active","placebo","ME")
        print(out)
    }
    invisible(NULL)
}

ResponderAnalysis <- function(x,mcid,type="absolute",conf.level=0.95){
    if(class(x) != "memix") stop("'x' must be 'ATT/ATE' by 'memixed'")
    type <- match.arg(tolower(type),c("absolute","relative"))
    if(conf.level <= 0||conf.level >= 1)
        stop("invalid 'mcid': must in [0,1]")
    alp2 <- 0.5-0.5*conf.level
    la <- rstan::extract(x$xfit, permuted = TRUE)
    if(type=="relative"){
        if(mcid<0||mcid>1) stop("invalid 'mcid': must in [0,1]")
        es1 <- abs(la$es_rel);
    }else{
        es1 <- abs(la$es_abs);
    }
    r <- round(mean(es1>mcid),4)
    ci <- as.numeric(quantile(es1, prob=c(alp2,1-alp2)))
    cat("\nPoint estimates:\n  mean =",mean(es1),"\n  median =", median(es1))
    cat("\nPosterior predictive probability:\n  P(theta >",mcid,") =", r, "\n")
    cat(100*conf.level,"% credible interval:\n  (", ci[1], ",",
        ci[2],")\n\n")
    invisible(NULL)
}

dratio <- function(x,mu,sigma){
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

## dratio(1,mu=c(1,2),sigma=c(2.5,3.5))

## use numerical method to compute the mean and SD
pratio <- function(x,mu,sigma, conf.level=0.95){
    stopifnot(conf.level>0&&conf.level<1)
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
    list(mean=mu,median=med,rate=Fx, ci=c(ll,ul))
}
## pratio(0.30,mu=c(1.6,6.5),sigma=c(1.5,2))
