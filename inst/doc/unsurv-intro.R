## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

## ----data---------------------------------------------------------------------
library(unsurv)

set.seed(2026)
n <- 150
Q <- 60
times <- seq(0, 5, length.out = Q)

group <- sample(1:3, n, TRUE, prob = c(0.35, 0.4, 0.25))
haz   <- c(0.18, 0.45, 0.8)[group]

S <- sapply(times, function(t) exp(-haz * t))
S <- S + matrix(rnorm(n * Q, 0, 0.02), nrow = n)
S[S < 0] <- 0
S[S > 1] <- 1

## ----fit----------------------------------------------------------------------
fit <- unsurv(S, times, K = NULL, K_max = 6, distance = "L2",
              enforce_monotone = TRUE, smooth_median_width = 5,
              standardize_cols = TRUE, eps_jitter = 0.0005)
fit

## ----medoids-plot, fig.cap = "Cluster medoid survival curves."----------------
plot(fit)

## ----autoplot, fig.cap = "Medoid curves via ggplot2 autoplot."----------------
library(ggplot2)
autoplot(fit)

## ----predict------------------------------------------------------------------
new_curves <- S[1:5, ]
predict(fit, new_curves)

## ----stability----------------------------------------------------------------
stab <- unsurv_stability(
  S, times, fit,
  B = 20, frac = 0.7,
  mode = "subsample",
  jitter_sd = 0.01,
  weight_perturb = 0.05,
  return_distribution = TRUE
  )
stab$mean

