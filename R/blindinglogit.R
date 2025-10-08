
.EstimateEffects <- function(fit, y, group, guess) {
    stopifnot(length(y) == length(group), length(y) == length(guess))
    N <- length(y)
    g0_idx <- which(group == 0L)
    g1_idx <- which(group == 1L)
    n0 <- length(g0_idx); n1 <- length(g1_idx)
    if (n0 == 0 || n1 == 0)
        stop("Both groups (0 and 1) must be present.")
    mu12 <- tapply(y,group,mean)
    tau.hat <- diff(mu12)
    y2 <- y - tau.hat * group

    draws <- as_draws_matrix(
        fit,
        variable = c("beta0","beta1","beta_group","beta_y",
                     "gamma0","gamma1","theta")
    )
    S <- nrow(draws)
  
    out <- matrix(NA_real_, nrow = S, ncol = 5)
    colnames(out) <- c("mean_g0","sd_g0","mean_g1","sd_g1","mean_diff")
  
    logit <- function(x) 1/(1+exp(-x))
    options(warn=-1)
    for (s in seq_len(S)) {
        beta0 <- draws[s, "beta0"]
        beta1 <- draws[s, "beta1"]
        beta_g <- draws[s, "beta_group"]
        beta_y <- draws[s, "beta_y"]
        gamma0 <- draws[s, "gamma0"]
        gamma1 <- draws[s, "gamma1"]
        th     <- draws[s, "theta"]
        
        eta0 <- beta0 + c(beta_g) * group + c(beta_y) * y
        eta1 <- beta1 + c(beta_g) * group + c(beta_y) * (y - th)
        pz   <- logit(gamma0 + c(gamma1) * y2)
    
        p0 <- logit(eta0)
        p1 <- logit(eta1)
        l0 <- dbinom(guess, size = 1, prob = p0)
        l1 <- dbinom(guess, size = 1, prob = p1)
        
        r  <- (pz * l1) / (pz * l1 + (1 - pz) * l0)
    
        m  <- y - th * r                     # E[y* | ...]
        v  <- (th^2) * r * (1 - r)           # Var[y* | ...] from latent z
    
        mean0 <- mean(m[g0_idx])
        mean1 <- mean(m[g1_idx])
    
        bvar0 <- if (n0 > 1) stats::var(m[g0_idx]) else 0
        bvar1 <- if (n1 > 1) stats::var(m[g1_idx]) else 0
        wvar0 <- mean(v[g0_idx])
        wvar1 <- mean(v[g1_idx])
        
        sd0 <- sqrt(bvar0 + wvar0)
        sd1 <- sqrt(bvar1 + wvar1)
        
        out[s, "mean_g0"]   <- mean0
        out[s, "sd_g0"]     <- sd0
        out[s, "mean_g1"]   <- mean1
        out[s, "sd_g1"]     <- sd1
        out[s, "mean_diff"] <- mean1 - mean0
    }
    options(warn=0)
    
    draws_df <- as.data.frame(out)
    
    summ <- function(x)
        c(mean = mean(x),
          `2.5%` = quantile(x, 0.025)[[1]],
          `97.5%` = quantile(x, 0.975)[[1]])
    summary <- rbind(
        mean_g0   = summ(draws_df$mean_g0),
        sd_g0     = summ(draws_df$sd_g0),
        mean_g1   = summ(draws_df$mean_g1),
        sd_g1     = summ(draws_df$sd_g1),
        mean_diff = summ(draws_df$mean_diff)
    )
    list(draws = draws_df, summary = as.data.frame(summary))
}

blinding.lslogit <- function(x, group, guess, mu0 = 0, s0 = 1, ...){
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
    p0 <- mean(la$theta > 0)
    p1 <- mean(la$gamma1 > 0)
    
    res <- .EstimateEffects(fit, y, group, guess)
    list(fit = fit, summary = res$summary, test.theta = p0, test.beta2=p1)
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
