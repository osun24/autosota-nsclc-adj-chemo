#!/usr/bin/env Rscript

# Locked R bridge for one outer-fold/seed causal and prognostic forest pair.
# Python owns data access, fit-only transforms, folds, and scoring. This bridge
# receives numeric matrices with deliberately generic x* column names so that
# neither outcomes nor treatment can accidentally enter the effect-modifier X.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 14) {
  stop("Expected 14 arguments: fit assess predictions importance seed trees sample_fraction min_node_size honesty_fraction alpha imbalance_penalty threads horizon fit_prognostic")
}

suppressPackageStartupMessages(library(grf))
required_grf <- "2.6.1"
if (as.character(packageVersion("grf")) != required_grf) {
  stop(sprintf("Locked grf version is %s; found %s", required_grf, packageVersion("grf")))
}

fit_path <- args[[1]]
assess_path <- args[[2]]
prediction_path <- args[[3]]
importance_path <- args[[4]]
seed <- as.integer(args[[5]])
num_trees <- as.integer(args[[6]])
sample_fraction <- as.numeric(args[[7]])
min_node_size <- as.integer(args[[8]])
honesty_fraction <- as.numeric(args[[9]])
alpha <- as.numeric(args[[10]])
imbalance_penalty <- as.numeric(args[[11]])
num_threads <- as.integer(args[[12]])
horizon <- as.numeric(args[[13]])
fit_prognostic <- as.integer(args[[14]]) == 1L

fit <- read.csv(fit_path, check.names = FALSE)
assess <- read.csv(assess_path, check.names = FALSE)
feature_names <- grep("^x[0-9]+$", names(fit), value = TRUE)
if (length(feature_names) < 1 || !identical(feature_names, names(assess)[grepl("^x[0-9]+$", names(assess))])) {
  stop("Fit and assessment feature matrices are empty or misaligned")
}

X_fit <- data.matrix(fit[, feature_names, drop = FALSE])
X_assess <- data.matrix(assess[, feature_names, drop = FALSE])
Y <- as.numeric(fit$Y)
W <- as.numeric(fit$W)
D <- as.numeric(fit$D)
W_hat <- as.numeric(fit$W_hat)
if (any(!is.finite(X_fit)) || any(!is.finite(X_assess)) ||
    any(!is.finite(Y)) || any(!is.finite(W_hat))) {
  stop("Non-finite input reached the locked grf bridge")
}
if (!all(W %in% c(0, 1)) || !all(D %in% c(0, 1))) {
  stop("W and D must be binary")
}
if (num_trees > 1000 || num_trees < 2 || num_trees %% 2 != 0) {
  stop("Tree count must be even and in [2, 1000]")
}

csf <- causal_survival_forest(
  X = X_fit,
  Y = Y,
  W = W,
  D = D,
  W.hat = W_hat,
  target = "RMST",
  horizon = horizon,
  failure.times = seq(0, horizon, length.out = 121),
  num.trees = num_trees,
  sample.fraction = sample_fraction,
  mtry = ncol(X_fit),
  min.node.size = min_node_size,
  honesty = TRUE,
  honesty.fraction = honesty_fraction,
  honesty.prune.leaves = TRUE,
  alpha = alpha,
  imbalance.penalty = imbalance_penalty,
  stabilize.splits = TRUE,
  ci.group.size = 2,
  tune.parameters = "none",
  compute.oob.predictions = FALSE,
  num.threads = num_threads,
  seed = seed
)
tau <- as.numeric(predict(csf, newdata = X_assess, num.threads = num_threads)$predictions)
importance <- as.numeric(variable_importance(csf))

# C-index is explicitly secondary, so only the first forest seed fits this
# separate model. A fixed grid makes its risk scale comparable across folds.
if (fit_prognostic) {
  grid <- seq(horizon / 60, horizon, length.out = 60)
  X_fit_prognostic <- cbind(X_fit, observed_ACT = W)
  X_assess_prognostic <- cbind(X_assess, observed_ACT = as.numeric(assess$W))
  sf <- survival_forest(
    X = X_fit_prognostic,
    Y = Y,
    D = D,
    failure.times = grid,
    num.trees = num_trees,
    sample.fraction = sample_fraction,
    mtry = ncol(X_fit_prognostic),
    min.node.size = min_node_size,
    honesty = TRUE,
    honesty.fraction = honesty_fraction,
    honesty.prune.leaves = TRUE,
    alpha = alpha,
    prediction.type = "Kaplan-Meier",
    compute.oob.predictions = FALSE,
    num.threads = num_threads,
    seed = seed + 100000L
  )
  survival <- predict(
    sf,
    newdata = X_assess_prognostic,
    failure.times = grid,
    prediction.times = "curve",
    prediction.type = "Kaplan-Meier",
    num.threads = num_threads
  )$predictions
  risk <- -rowMeans(survival)
} else {
  risk <- rep(NA_real_, nrow(X_assess))
}

write.csv(data.frame(tau = tau, risk = risk), prediction_path, row.names = FALSE)
write.csv(data.frame(feature_index = seq_along(importance), importance = importance), importance_path, row.names = FALSE)
