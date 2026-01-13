# Minimal runner: reads params.yml and prints results
if(!requireNamespace("yaml", quietly = TRUE)) install.packages("yaml")
library(yaml)

source("R/power_functions.R")
p <- yaml::read_yaml("config/params.yml")

# AIM 1
maf_vec <- unlist(p$aim1$maf_vec)
alpha_gwas <- p$aim1$alpha_gwas
target_power <- p$aim1$target_power
coh <- p$aim1$cohorts

aim1_tbl <- do.call(rbind, lapply(names(coh), function(nm){
  x <- coh[[nm]]
  Ncase_en <- round(x$Ncase_total * x$subtype_frac)
  Ncase_non <- x$Ncase_total - Ncase_en
  out <- expand.grid(cohort=nm, maf=maf_vec, comparison=c("PD_met_enriched_vs_ctrl","PD_non_enriched_vs_ctrl"))
  out$Ncase <- ifelse(out$comparison=="PD_met_enriched_vs_ctrl", Ncase_en, Ncase_non)
  out$Nctrl <- x$Nctrl
  out$MDE_OR_80pct <- mapply(function(nc,n0,mf) gwas_mde_or(nc,n0,mf,alpha_gwas,target_power),
                            out$Ncase,out$Nctrl,out$maf)
  out
}))
print(aim1_tbl)

# Meta Neff & MDE
Ncase_meta <- sapply(coh, function(x) round(x$Ncase_total * x$subtype_frac))
Nctrl_meta <- sapply(coh, function(x) x$Nctrl)
Neff_meta <- meta_neff(Ncase_meta, Nctrl_meta)
cat("\nMeta Neff (enriched vs ctrl) =", round(Neff_meta), "\n")
meta_mde <- sapply(maf_vec, function(m) meta_mde_or(Neff_meta, m, alpha_gwas, target_power))
print(data.frame(maf=maf_vec, MDE_OR_80pct=meta_mde))

# AIM 1 MR placeholder
mr <- p$aim1_mr
print(mr_power_placeholder(Neff_meta, mr$R2_total, mr$OR_causal, mr$alpha_mr))

# AIM 2
a2 <- p$aim2
pw2 <- sapply(unlist(a2$d_vec), function(d) t_test_power(a2$n_enriched, a2$n_non_enriched, d, a2$alpha_prot))
print(data.frame(d=unlist(a2$d_vec), power=pw2))

# AIM 3 incident
a3i <- p$aim3_incident
pw3i <- cox_power_from_events(a3i$events_incident_pd, a3i$HR_target, a3i$alpha, a3i$p_exposed)
need3i <- events_required_cox(a3i$HR_target, a3i$alpha, 0.8, a3i$p_exposed)
print(list(power_incident=pw3i, events_required_for_80pct=need3i))

# AIM 3 LEDD
a3l <- p$aim3_ledd
pw3l <- ledd_slope_power(a3l$n_exposed_pd, a3l$n_comp_pd, a3l$delta_slope, a3l$sigma_slope, a3l$alpha)
print(list(power_ledd=pw3l))

# AIM 3 interaction
a3x <- p$aim3_interaction
pw3x <- interaction_power_two_strata(a3x$events_stratum1, a3x$events_stratum2, a3x$HR_s1, a3x$HR_s2,
                                     a3x$alpha, a3x$p_exposed_s1, a3x$p_exposed_s2)
print(pw3x)
