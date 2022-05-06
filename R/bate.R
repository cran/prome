## 04/25/2022: using lognormal SD results in similar results. 
## Unconstrained SD doesn't improve the convengence much. We 
## use inv_gamma for SDs and changed the default stepsize to 5.

## 04/08/2022: (1) change function name from 'memixed' to 'bate'. (2)
## apply the bayesian method on the subject level mean values to
## reduce the variance of the estimates. This will also reduce the
## computing time.


## 03/27/2022 In ver x.11, ATE and ATT are for two
## different parameterizations of the treatment effects. For RCTs,
## there is little difference between the two methods. ATT is more
## flexible as the variance of control arm could be larger than the
## variance of the treatment arm. We keep the two version in the new
## version, but the argument of the memixed function will be changed
## to avoid confusion.

##     control = list(adapt_delta = 0.8,
##               stepsize = 5,
##               max_treedepth = 10),


bate <- function(x0,x1,group,z,x.range,...){
    if(!missing(x.range)){
        if(!is.numeric(x.range))
            stop("'x.range' must be numeric")
        if(length(x.range) != 2)
            stop("'x.range' must have length 2")
        a <- x.range[1]; b <- x.range[2]
        xrng <- b - a
        if(xrng <= 0) stop("invalid range of data")
        if(any(x0<a)||any(x1<a)||any(x0>b)||any(x1>b)){
            warning("some measures fall beyond the data range")
            x0[x0 < a] <- a
            x1[x1 < a] <- a
            x0[x0 > b] <- b
            x1[x1 > b] <- b
        }
    }
    x0 <- as.matrix(x0); 
    x1 <- as.matrix(x1);
    xbar0 <- apply(x0,1,mean,na.rm=TRUE)
    xbar1 <- apply(x1,1,mean,na.rm=TRUE)
    mu0 <- mean(xbar0)
    s0 <- sd(c(x0))
    n <- nrow(x0)
    if(nrow(x1) != n)
        stop("different numbers of subjects in 'x1' and 'x0'")

    if(missing(group)){
        narm <- 1
    }else{
        fgrp <- as.factor(group)
        narm <- nlevels(fgrp)
        ngrp <- as.numeric(fgrp)-1
        if(length(fgrp) != n)
            stop("different numbers of subjects in 'x0' and 'group'")
    }
    
    out <- NULL
    if(narm == 1){
        stop("Algorithm for 1-arm data is not supported in this version")
    }else if(narm == 2){
        if(missing(z)){
            sampledat <- list(N=n,group=ngrp,
                              r1=ncol(x1),r0=ncol(x0),
                              w0=t(x0),w1=t(x1),
                              mu=mu0,sigma=s0)
            fit3 <- rstan::sampling(stanmodels$ate2m,
                                    data = sampledat,
                                    iter=5000*2,
                                    refresh=0,...)
            
            la <- rstan::extract(fit3, permuted = TRUE)
            es0 <- mean(la$es0); 
            es1 <- mean(la$es1); 
            mu0 <- mean(la$mu0) 
            sig0 <- mean(la$sigma0) 
            tau0 <- mean(la$tau) 
            s1 <- mean(la$sa)
            s0 <- mean(la$ss)
            out <- structure(
                list(xfit=fit3,
                     mean0.t0=mu0,sig0.t0=sqrt(sig0),
                     sig.me=sqrt(tau0),
                     mu.active=es1,sig.active=sqrt(s1),
                     mu.sham=es0,sig.sham=sqrt(s0),
                     groups=levels(fgrp)),
                class='ate')
        }else{
            znam <- names(z)
            z <- as.matrix(z)
            dv0 <- (x0 - xbar0)^2
            tau <- sqrt(mean(dv0))
            if(nrow(z) != n)
                stop("different numbers of subjects in 'z' and 'group'")
            sampledat <- list(N=n,group=ngrp,
                              r1=ncol(x1),rz=ncol(z),
                              w0=xbar0,w1=xbar1,tau=tau,z=t(z))
            fit4 <- rstan::sampling(stanmodels$zate,
                                    data = sampledat,
                                    iter=5000*2,
                                    refresh=0,...)
            la <- rstan::extract(fit4, permuted = TRUE)
            es0 <- mean(la$es0); 
            es1 <- mean(la$es1); 
            mu0 <- mean(la$mu0) 
            sig0 <- median(la$sigma0) 
                                        #tau <- mean(la$tau)
            s1 <- median(la$sig_active)
            s0 <- median(la$sig_sham)
            beta.hat <- apply(la$beta,2,mean,na.rm=TRUE)
            k <- length(beta.hat)
            ci <- NULL
            for(i in 1:k){
                tmp <- as.numeric(quantile(la$beta[,i],
                                           prob=c(0.025,0.975)))
                ci <- rbind(ci, tmp)
            }
            out.coef <- data.frame(beta.hat,ci)
            names(out.coef) <- c("coef","95.CI.ll","95.CI.ul")
            rownames(out.coef) <- znam
            out <- structure(
                list(xfit=fit4,mean0.t0=mu0,
                     sig0.t0=sig0,
                     sig.me=tau,mu.active=es1,sig.active=s1,
                     mu.sham=es0,sig.sham=s0,groups=levels(fgrp),
                     coef=out.coef),
                class='ate')
        }
    }else{
        stop("Algorithm doesn't support data from trials with >2 arms")
    }
    out
}

print.att <- function(x,...){
    out <- data.frame(mean=c(x$mean0.t0,x$mu.active,x$mu.sham,0),
                      sd=c(x$sig0.t0,x$sig.active,x$sig.sham,x$sig.me))
    rownames(out) <- c("baseline","active","placebo","ME")
    print(out)
    invisible(out)
}

print.ate <- function(x,...){
    tmp1 <- c(x$mean0.t0,x$mu.active,x$mu.sham,0)
    tmp2 <- c(x$sig0.t0,x$sig.active,x$sig.sham,x$sig.me)
    out <- data.frame(mean=tmp1,sd=tmp2)
    rownames(out) <- c("baseline","treat","control","M.E.")
    print(out)
    if(!is.null(x$coef)){
        cat("\ncoefficient table:\n")
        print(x$coef)
    }
    invisible(out)
}

ResponderAnalysis <- function(x,mcid,type="absolute",conf.level=0.95,show=TRUE){
    type <- match.arg(tolower(type),c("absolute","relative"))
    if(conf.level <= 0||conf.level >= 1)
        stop("invalid 'mcid': must in [0,1]")
    alp2 <- 0.5-0.5*conf.level
    la <- rstan::extract(x$xfit, permuted = TRUE)
    if(type=="relative"){
        if(mcid<0||mcid>1) stop("invalid 'mcid': must in [0,1]")
        es1 <- la$es_rel;
    }else{
        es1 <- la$es_abs;
    }
    r <- round(mean(abs(es1)>mcid),4)
    ci <- as.numeric(quantile(es1, prob=c(alp2,1-alp2)))
    if(show){
        cat("\nPoint estimates:\n  mean =",mean(es1),"\n  median =", median(es1))
        cat("\n",100*conf.level,"% credible interval:\n  (", ci[1], ",",
            ci[2],")")
        cat("\nPosterior probability:\n  P(|theta| > ",mcid,") = ", r, "\n\n",sep='')
    }
    out <- list(mean=mean(es1),median=median(es1),ll=ci[1],ul=ci[2])
    invisible(out)
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
