
blinding.test <- function(x, group, guess, mu0 = 0, s0 = 1,...){
    y <- x
    if (length(y) != length(group) || length(y) != length(guess)) {
        stop("y, group, and guess must have the same length.")
    }
    if (!all(group %in% c(0, 1))) {
        if (is.logical(group)) group <- as.integer(group)
        if (!all(group %in% c(0, 1))) {
            stop("group must be a 0/1 vector (or coercible to 0/1).")
        }
    }
    if (!all(guess %in% c(0, 1))) {
        if (is.logical(guess)) guess <- as.integer(guess)
        if (!all(guess %in% c(0, 1))) {
            stop("guess must be a 0/1 vector (or coercible to 0/1).")
        }
    }
    if (s0 <= 0) stop("s0 must be > 0 for the Normal prior on theta.")
    mu12 <- tapply(y,group,mean)
    tau.hat <- diff(mu12)
    y2 <- y - tau.hat * group
    
    stan_data <- list(
        N     = length(y),
        guess = as.integer(guess),
        group = as.integer(group),
        y     = as.vector(y),
        z     = as.vector(y2),
        mu0   = as.numeric(mu0),
        s0    = as.numeric(s0)
    )
    
    fit <- rstan::sampling(stanmodels$blindlogit,
        data    = stan_data,
        refresh = 0, ...)
    la <- rstan::extract(fit, permuted = TRUE)
    theta0 <- mean(la$theta > 0)
    cdist <- get(load(system.file("extdata/ctest.rda",
                                  package = "prome")))
    ddist <- get(load(system.file("extdata/dtest.rda",
                                  package = "prome")))
    edist <- get(load(system.file("extdata/etest.rda",
                                  package = "prome")))
    p0 <- mean(theta0 < cdist)
    p1 <- mean(theta0 < ddist)
    p2 <- mean(theta0 < edist)
    beta2hat <- mean(la$gamma1 > 0)
    out2 <- .BItest(x=x,group=group,guess=guess,iter=999)
    out <- list(fit = fit, p.value = (p0+p1+p2)/3,
                p1.5 = p0, p2.5 = p1, p3.5 = p2,
                BI.test = out2,
                test.theta = theta0, test.beta2=beta2hat)
    structure(out,class="ctest")
}

##iter = 2000,
##warmup = floor(iter / 2),
##chains = 4,
##seed = 1234,
##adapt_delta = 0.95,
##max_treedepth = 12
##) {

##iter    = iter,
##warmup  = warmup,
##chains  = chains,
##seed    = seed,
##control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
##refresh max(1L, iter %/% 10L) # print progress every ~10%

.BIEst <- function(group,guess){
  arm <- group
  arm[arm==0] = "B:Placebo"
  arm[arm==1] = "A:Treatment"
  blind = guess
  blind[blind!=1] <- "2:Placebo"
  blind[blind==1] <- "1:Treatment"
  xtbl <- table(blind, arm)
  xtbl2 <- rbind(xtbl, "3:DNK"=c(0,0))
  res <- BI::BI(xtbl2)
}

.subptbBI <- function(x,db){
  if(x==1){
    db$guess <- sample(db$guess,replace=FALSE)
  }
  BIres = .BIEst(group=db$group,guess=db$guess)
  lm0 <- lm(y ~ guess + group, data=db)
  uadj <- lm0$coef[[3]]
  usham <- lm0$coef[[2]]

  res <- c(BIres$JamesBI[1,1],
           BIres$BangBI[2,1],
           BIres$BangBI[1,1],
           uadj, usham)
  return(res)
}

.BItest <- function(x, group, guess,iter=200){
  grp <- group
  sele <- grp == 0
  grp[sele] <- "control"
  grp[!sele] <- "treatment"
  gus <- guess
  sele <- gus == 0
  gus[sele] <- "control"
  gus[!sele] <- "treatment"
  
  tbl <- table(group=grp,guess=gus)
  iter <- round(iter)
  res <- NULL
  if(iter > 30){
      mydata = data.frame(y = x, group = group, guess = guess)
      res1 = .subptbBI(0,db=mydata)
      z = matrix(1,nrow=iter,ncol=1)
      out = apply(z,1,FUN=.subptbBI,db=mydata)
      res2 = apply(out,1,quantile,probs=c(0.0125,0.025,0.5,0.975,0.9875))
      nams = rownames(res2)
      res <- data.frame(estimate=res1,t(res2))
      rownames(res) <- c("JBI","BBI.ctrl","BBI.treat","EffSize","Sham")
      names(res) <- c("estimate",nams)
  }

  list(data = tbl, BI = res)
}

print.ctest <- function(x,...){
    cat("\nData summary:\n")
    print(x$BI.test$data)
    cat("\nFisher randomization test:\n")
    print(round(x$BI.test$BI,3))
    cat("\nLatent-shift test of contamination:\n  p-value = ",
        x$p.value,"\n\n")
}

