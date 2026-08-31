## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)
has_survival <- requireNamespace("survival", quietly = TRUE)

## ----sim-data-----------------------------------------------------------------
library(unsurv)

set.seed(2026)
n <- 150
Q <- 60
sim_times <- seq(0, 5, length.out = Q)

group <- sample(1:3, n, TRUE, prob = c(0.35, 0.4, 0.25))
haz   <- c(0.18, 0.45, 0.8)[group]

S <- sapply(sim_times, function(t) exp(-haz * t))
S <- S + matrix(rnorm(n * Q, 0, 0.02), nrow = n)
S[S < 0] <- 0
S[S > 1] <- 1

## ----sim-fit------------------------------------------------------------------
sim_fit <- unsurv(S, sim_times, K = NULL, K_max = 6, distance = "L2",
                  enforce_monotone = TRUE, smooth_median_width = 5,
                  standardize_cols = TRUE, eps_jitter = 0.0005, seed = 1)
sim_fit

## ----sim-plot, fig.cap = "Cluster medoid survival curves (simulated data)."----
plot(sim_fit)

## ----sim-predict--------------------------------------------------------------
predict(sim_fit, S[1:5, ])

## ----metabric-load------------------------------------------------------------
f <- system.file("extdata", "metabric_survdnn_curves.rds", package = "unsurv")
mb <- readRDS(f)

str(mb, max.level = 1)
cat(mb$provenance)

## ----metabric-fit-------------------------------------------------------------
fit <- unsurv(mb$S_partition, mb$times, K = 3, distance = "L2",
              enforce_monotone = TRUE, smooth_median_width = 3,
              standardize_cols = FALSE, eps_jitter = 0, seed = 20260615)
fit

## ----metabric-medoids, fig.cap = "Predicted survival curves coloured by unsurv cluster, with medoid prototypes (thick lines).", eval = has_survival----
library(ggplot2)

grid_n  <- length(mb$times)
curf <- data.frame(
  id       = rep(seq_len(nrow(mb$S_partition)), each = grid_n),
  time     = rep(mb$times, times = nrow(mb$S_partition)),
  survival = as.vector(t(mb$S_partition)),
  cluster  = factor(rep(fit$clusters, each = grid_n))
)
medf <- data.frame(
  time     = rep(mb$times, times = nrow(fit$medoids)),
  survival = as.vector(t(fit$medoids)),
  cluster  = factor(rep(seq_len(nrow(fit$medoids)), each = grid_n))
)

ggplot(curf, aes(time, survival, group = id, colour = cluster)) +
  geom_line(linewidth = 0.25, alpha = 0.20) +
  geom_line(data = medf, aes(group = cluster), linewidth = 1.1) +
  labs(x = "Months", y = "Predicted survival probability", colour = "Cluster") +
  theme_minimal(base_size = 11)

## ----metabric-predict---------------------------------------------------------
val_clusters <- predict(fit, mb$S_validation)
table(validation = val_clusters)

## ----metabric-baseline--------------------------------------------------------
h <- which.min(abs(mb$times - stats::median(mb$times)))
risk_p <- 1 - mb$S_partition[, h]
risk_v <- 1 - mb$S_validation[, h]

pam_scalar <- cluster::pam(as.matrix(risk_p), k = 3)
centres <- tapply(risk_p, pam_scalar$clustering, mean)
scalar_v <- vapply(risk_v, function(v) which.min(abs(v - centres)), integer(1))

## ----metabric-compare-partition, eval = has_survival--------------------------
cmp_part <- unsurv_compare(
  list(`unsurv curve` = fit$clusters, `scalar risk` = pam_scalar$clustering),
  mb$os_time$partition, mb$os_event$partition,
  reference = "unsurv curve"
)
cmp_part$summary
cmp_part$cluster_summary

## ----metabric-compare-validation, eval = has_survival-------------------------
cmp_val <- unsurv_compare(
  list(`unsurv curve` = val_clusters, `scalar risk` = scalar_v),
  mb$os_time$validation, mb$os_event$validation,
  reference = "unsurv curve"
)
cmp_val$cluster_summary

## ----metabric-km, fig.cap = "Kaplan-Meier curves by unsurv cluster on the validation set.", eval = has_survival----
autoplot(unsurv_compare(
  list(`unsurv curve` = val_clusters),
  mb$os_time$validation, mb$os_event$validation
))

## ----metabric-stability-------------------------------------------------------
stab <- unsurv_stability(
  mb$S_partition, mb$times, fit,
  B = 30, frac = 0.7, mode = "subsample",
  jitter_sd = 0.005, weight_perturb = 0.1, eps_jitter = 0,
  return_distribution = TRUE
)
stab$mean

