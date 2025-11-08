#' @import Rcpp
#' @importFrom RcppArmadillo armadillo_throttle_cores
#'
blinding.cpe <- function(x, group, guess){
    if (!(is.numeric(x) && is.null(dim(x))))
        stop("'x' must be a numeric vector, not a matrix or other type.")
    n <- length(x)
    if(length(group) != n)
        stop("'x' and 'group' lengths differ")
    if(length(guess) != n)
        stop("'x' and 'guess' lengths differ")
    group <- as.factor(group)
    if(nlevels(group) != 2)
        stop("'group' must have two levels")
    group <- as.numeric(group) - 1

    guess <- as.factor(guess)
    if(nlevels(guess) != 2)
        stop("'guess' must have two levels")
    guess <- as.numeric(guess) - 1

    lm0 <- lm(x ~ guess + group)
    uadj <- lm0$coef[[3]]
    lmout1 <- lm(x ~ group)
    uclass <- lmout1$coef[[2]]

    out <- .cpe2(x=x,group=group,guess=guess)
    out$classic <- uclass
    out$adjusted <- uadj

    structure(out,class="cpe")
}

.cpe2 <- function(x, group, guess){
    mydata  <-  data.frame(y = x, group = group, guess = guess)
    lmout1 <- lm(y ~ group, data = mydata)
    res1 <- lmout1$coef[[2]]
    lmout2 <- lm(y ~ group + guess, data = mydata)
    res2 <- lmout2$coef[[2]]
    usham <- max(abs(res1),abs(res2),abs(res1+lmout1$coef[[1]]))

    out <- NULL
    m <- 400
    d <- seq(0,2.5,length=m)*usham
    b1 <- sapply(d,FUN=.cpe2dev,db=mydata)
    out$x <- d
    out$y <- b1
    usham2 <- .cpfind(d,b1)[2]
    y2 <- x - usham2 * guess
    lm0 <- lm(y2~group)
    tau.hat <- lm0$coef[[2]]; 
    ctrl.hat <- lm0$coef[[2]];
    out$control <- ctrl.hat
    out$active <- ctrl.hat + tau.hat
    out$sham <- usham2
    out$CPE <- tau.hat
    out
}

.cpe2dev <- function(d,db){
    db$y <- db$y - d * db$guess
    options(warn = -1)
    g1 <- glm(guess ~ y, family = binomial(), data=db)
    options(warn = 0)
    p <- fitted(g1)
    ## Clip extreme probabilities to avoid log(0)
    eps <- 1e-12
    p_clipped <- p
    p_clipped[p < eps] <- eps
    p_clipped[p > 1 - eps] <- 1 - eps
    ## Compute stabilized deviance
    dev_stable <- -2 * sum(db$guess * log(p_clipped)
                           + (1 - db$guess) * log(1 - p_clipped))
    dev_stable
}

.cpfind <- function(x,y){
  n <- length(x)
  xm <- which.max(y)
  k2 = max(30, xm)
  x = x[1:xm]
  y = y[1:xm]

  n <- length(x)
  xm <- which.max(y)
  sele = y <= y[xm]
  if(sum(sele) > 30){
    if(n-xm > 10){
      x = x[1:(xm+10)]
      y = y[1:(xm+10)]
    }else{
      x = x[1:xm]
      y = y[1:xm]
    }
  }else{
    if(n-xm > 20){
      x = x[1:(xm+20)]
      y = y[1:(xm+20)]
    }
  }
  n <- length(x)
  if(length(y) != n) stop("'x' and 'y' length differ")
  if(n<10) stop("too few data points to find changing point")
  ## starting point cannot be < 4
  lm2 <- lm(y~x)
  RSS2 <- sum((lm2$resi)^2)
  k <- max(round(n*0.1),4)
  y2 <- NULL
  x2 <- NULL
  for(i in k:(n-k-1)){
    y0 <- y[1:i]
    x0 <- x[1:i]
    x02 <- x0^2
    lm0 <- lm(y0~x0+x02)
    RSS0 <- sum((lm0$resi)^2)
    y1 <- y[(i+1):n]
    x1 <- x[(1+i):n]
    x12 <- x1^2
    lm1 <- lm(y1~x1+x12)
    RSS1 <- sum((lm1$resi)^2)
    RSS <- RSS0 + RSS1
    x2 <- c(x2,i)
    y2 <- c(y2,RSS)
  }
  xid <- which.min(y2)[1]
  
  c(x2[xid],x[xid])
}

print.cpe <- function(x,...){
    cat("\nMean effect (treatment) = ", x$active,
      "\nMean effect (control) = ", x$control,
      "\nEffect size = ", x$CPE,
      "\nSham effect = ", x$sham,
      "\nClassic estimate = ", x$classic,
      "\nUnblinding adjusted estimate = ", x$adjusted, "\n\n")
}

