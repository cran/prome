## 06/16/2023: To analyze data from AB/BA crossover design.
## remove all missing values. will be updated later.

.xosum <- function(x){
    tout <- t.test(x)
    sout <- shapiro.test(x)
    tmp <- c(mean(x),sd(x),tout$conf.int[1],
             tout$conf.int[1],sout$p.value)
    c(length(x),round(tmp,3))
}

xover <- function(group,y1,y2,y0,...){
    if(missing(y0)){
        out <- xover1(group=group,y1=y1,y2=y2,...)
    }else{
        out <- xover0(group=group,y1=y1,y2=y2,y0=y0,...)
    }
    out
}

xover1 <- function(group,y1,y2,...){
    n <- length(group)
    if(length(y1) != n)
        stop("'group' and 'y1' lengths differ")
    if(length(y2) != n)
        stop("'group' and 'y2' lengths differ")
    ## two levels for AB/BA trials ######################
    group <- as.factor(group)
    if(nlevels(group) != 2)
        stop("Only two groups (treatment sequences) are allowed")
    sele.na <- is.na(y1) | is.na(y2) | is.na(group)
    if(any(sele.na)){
        if(sum(!sele.na) < 6){
            stop("sample too small, too many missing values")
        }else{
            warning("missing value(s) removed!")
            y1 <- y1[!sele.na]
            y2 <- y2[!sele.na]
            group <- as.character(group[!sele.na])
            group <- as.factor(group)
        }
    }
    group1 <- group == levels(group)[1]
    n1 <- sum(group1); n2 <- sum(!group1)
    if(n1 < 3) stop("group 1 too small")
    if(n2 < 3) stop("group 2 too small")
    y11 <- y1[group1]
    y12 <- y2[group1]
    y21 <- y1[!group1]
    y22 <- y2[!group1]

    x <- c(y11,y12,y21,y22)
    mu <- mean(x); sigma <- sd(x);

    xout1 <- NULL
    tmp <- .xosum(y11); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y12); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y21); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y22); xout1 <- rbind(xout1,tmp)
    xout1 <- as.data.frame(xout1)
    rownames(xout1) <- c("group(1), period=1",
                         "group(1), period=2",
                         "group(2), period=1",
                         "group(2), period=2")
    names(xout1) <- c("n","mean","sd","ll.95","ul.95","Shapiro")
    
    sampledat <- list(n1=n1, n2=n2,
                      mu = mu, sig = sigma,
                      y11=y11,y12=y12,
                      y21=y21,y22=y22)
    fit4 <- rstan::sampling(stanmodels$cross1,
                            data = sampledat,
                            iter=5000*2,
                            refresh=0,...)
    la <- rstan::extract(fit4, permuted = TRUE)
    out <- NULL; out2 <- NULL

    tmp <- summary(la$tau_d);
    qtl <- quantile(la$tau_d, prob=c(0.025,0.975))
    out <- rbind(out,tmp)
    out2 <- rbind(out2,qtl)

    tmp <- summary(la$pi_d);
    qtl <- quantile(la$pi_d, prob=c(0.025,0.975))
    out <- rbind(out,tmp)
    out2 <- rbind(out2,qtl)

    tmp <- summary(la$lambda_d);
    qtl <- quantile(la$lambda_d, prob=c(0.025,0.975))
    out <- rbind(out,tmp)
    out2 <- rbind(out2,qtl)
    xout2 <- cbind(out, out2)
    rownames(xout2) <- c("tau(2-1)","period(2-1)","lambda(2-1)")

    out <- structure(
        list(stat = xout1, best=xout2, xfit = fit4), class='xover')
}

xover0 <- function(group,y1,y2,y0,...){
    n <- length(y0)
    if(length(y1) != n)
        stop("'y0' and 'y1' lengths differ")
    if(length(y2) != n)
        stop("'y0' and 'y2' lengths differ")
    if(length(group) != n)
        stop("'y0' and 'group' lengths differ")
    ## two levels for AB/BA trials ######################
    group <- as.factor(group)
    if(nlevels(group) != 2)
        stop("Only two groups (treatment sequences) are allowed")
    sele.na <- is.na(y0) | is.na(y1) | is.na(y2) | is.na(group)
    if(any(sele.na)){
        if(sum(!sele.na) < 6){
            stop("sample too small, too many missing values")
        }else{
            warning("missing value(s) removed!")
            y0 <- y0[!sele.na]
            y1 <- y1[!sele.na]
            y2 <- y2[!sele.na]
            group <- group[!sele.na]
        }
    }
    group1 <- group == levels(group)[1]
    n1 <- sum(group1); n2 <- sum(!group1)
    if(n1 < 3) stop("group 1 too small")
    if(n2 < 3) stop("group 2 too small")
    y10 <- y0[group1]
    y11 <- y1[group1]
    y12 <- y2[group1]
    y20 <- y0[!group1]
    y21 <- y1[!group1]
    y22 <- y2[!group1]

    xout1 <- NULL
    tmp <- .xosum(y10); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y11); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y12); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y20); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y21); xout1 <- rbind(xout1,tmp)
    tmp <- .xosum(y22); xout1 <- rbind(xout1,tmp)
    xout1 <- as.data.frame(xout1)
    rownames(xout1) <- c("group(1), period=0","group(1), period=1",
                         "group(1), period=2","group(2), period=0",
                         "group(2), period=1","group(2), period=2")
    names(xout1) <- c("n","mean","sd","ll.95","ul.95","Shapiro")
    
    sampledat <- list(n1=n1, n2=n2,y10=y10,y20=y20,
                      y11=y11,y12=y12,y21=y21,y22=y22)
    fit4 <- rstan::sampling(stanmodels$cross,
                            data = sampledat,
                            iter=5000*2,
                            refresh=0,...)
    la <- rstan::extract(fit4, permuted = TRUE)
    out <- NULL; out2 <- NULL

    tmp <- summary(la$mu);
    qtl <- quantile(la$mu, prob=c(0.025,0.975))
    out <- rbind(out,tmp)
    out2 <- rbind(out2,qtl)

    tmp <- summary(la$tau_d);
    qtl <- quantile(la$tau_d, prob=c(0.025,0.975))
    out <- rbind(out,tmp)
    out2 <- rbind(out2,qtl)

    tmp <- summary(la$pi_d);
    qtl <- quantile(la$pi_d, prob=c(0.025,0.975))
    out <- rbind(out,tmp)
    out2 <- rbind(out2,qtl)

    tmp <- summary(la$lambda_d);
    qtl <- quantile(la$lambda_d, prob=c(0.025,0.975))
    out <- rbind(out,tmp)
    out2 <- rbind(out2,qtl)
    xout2 <- cbind(out, out2)
    rownames(xout2) <- c("mean(0)","tau(2-1)","period(2-1)","lambda(2-1)")

    out <- structure(
        list(stat = xout1, best=xout2, xfit = fit4), class='xover')
}

print.xover <- function(x,...){
    cat("\nSummary statistics:\n")
    print(x$stat)
    cat("\n\nBayesian estimates:\n")
    print(round(x$best,3))
}

