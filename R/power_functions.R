logit <- function(p) log(p/(1-p))

gwas_power_cc <- function(Ncase, Nctrl, maf, OR, alpha=5e-8) {
  if(any(maf <= 0 | maf >= 0.5)) stop("maf must be in (0,0.5)")
  Neff <- 4 / (1/Ncase + 1/Nctrl)
  varG <- 2*maf*(1-maf)
  se <- sqrt(1 / (Neff * varG))
  beta <- log(OR)
  zcrit <- qnorm(1 - alpha/2)
  ncp <- beta / se
  power <- pnorm(-zcrit - ncp) + (1 - pnorm(zcrit - ncp))
  list(power=power, Neff=Neff, se=se, ncp=ncp)
}

gwas_mde_or <- function(Ncase, Nctrl, maf, alpha=5e-8, target_power=0.8,
                        OR_grid = exp(seq(log(1.02), log(2.0), length.out=2000))) {
  pw <- sapply(OR_grid, function(or) gwas_power_cc(Ncase,Nctrl,maf,or,alpha)$power)
  idx <- which(pw >= target_power)
  if(length(idx)==0) return(NA)
  OR_grid[min(idx)]
}

meta_neff <- function(Ncase_vec, Nctrl_vec) {
  sum(4 / (1/Ncase_vec + 1/Nctrl_vec))
}

meta_mde_or <- function(Neff, maf, alpha=5e-8, target_power=0.8) {
  varG <- 2*maf*(1-maf)
  se <- sqrt(1/(Neff*varG))
  zcrit <- qnorm(1-alpha/2)
  f <- function(beta) {
    ncp <- beta/se
    pw <- pnorm(-zcrit - ncp) + (1 - pnorm(zcrit - ncp))
    pw - target_power
  }
  beta_hat <- uniroot(f, interval=c(log(1.01), log(3.0)))$root
  exp(beta_hat)
}

mr_power_placeholder <- function(Neff_outcome, R2_total, OR_causal, alpha=0.05) {
  zcrit <- qnorm(1-alpha/2)
  z <- abs(log(OR_causal)) * sqrt(Neff_outcome * R2_total)
  power <- pnorm(z - zcrit)
  list(power=power, z=z)
}

t_test_power <- function(n1, n2, d, alpha=0.05) {
  se <- sqrt(1/n1 + 1/n2)
  zcrit <- qnorm(1-alpha/2)
  ncp <- d / se
  pnorm(-zcrit - ncp) + (1 - pnorm(zcrit - ncp))
}

events_required_cox <- function(HR, alpha=0.05, power=0.8, p_exposed=0.5) {
  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)
  v <- p_exposed*(1-p_exposed)
  (z_alpha + z_beta)^2 / ((log(HR))^2 * v)
}

cox_power_from_events <- function(events, HR, alpha=0.05, p_exposed=0.5) {
  z_alpha <- qnorm(1 - alpha/2)
  v <- p_exposed*(1-p_exposed)
  z <- sqrt(events * v) * abs(log(HR))
  pnorm(z - z_alpha)
}

ledd_slope_power <- function(n1, n2, delta_slope, sigma_slope, alpha=0.05) {
  se <- sigma_slope * sqrt(1/n1 + 1/n2)
  zcrit <- qnorm(1-alpha/2)
  ncp <- abs(delta_slope)/se
  pnorm(-zcrit - ncp) + (1 - pnorm(zcrit - ncp))
}

interaction_power_two_strata <- function(events_s1, events_s2, HR_s1, HR_s2,
                                         alpha=0.05, p_exposed_s1=0.5, p_exposed_s2=0.5) {
  v1 <- p_exposed_s1*(1-p_exposed_s1)
  v2 <- p_exposed_s2*(1-p_exposed_s2)
  se1 <- 1 / sqrt(events_s1 * v1)
  se2 <- 1 / sqrt(events_s2 * v2)
  delta <- log(HR_s1) - log(HR_s2)
  se_delta <- sqrt(se1^2 + se2^2)
  zcrit <- qnorm(1 - alpha/2)
  ncp <- abs(delta) / se_delta
  power <- pnorm(-zcrit - ncp) + (1 - pnorm(zcrit - ncp))
  list(power=power, delta=delta, se_delta=se_delta, ncp=ncp)
}
