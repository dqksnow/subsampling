make_rf_logit_data <- function(seed = 1, N = 3000) {
  set.seed(seed)
  Z1 <- rbinom(N, 1, 0.06)
  Z2 <- rbinom(N, 1, 0.08)
  X1 <- rnorm(N)
  X2 <- rnorm(N)
  eta <- -0.4 + 0.1 * Z1 + 0.08 * Z2 + 0.15 * X1 - 0.1 * X2
  Y <- rbinom(N, 1, plogis(eta))
  data.frame(Y = Y, Z1 = Z1, Z2 = Z2, X1 = X1, X2 = X2)
}

check_basic_rf_fit <- function(object, p, info_prefix) {
  expect_true(inherits(object, "list"),
              info = paste(info_prefix, "returns a list."))
  expect_true(inherits(object, "ssp.glm.rF"),
              info = paste(info_prefix, "has class 'ssp.glm.rF'."))
  expect_equal(length(object$coef.ssp), p,
               info = paste(info_prefix, "returns one coefficient per model column."))
  expect_true(all(is.finite(object$coef.ssp)),
              info = paste(info_prefix, "has finite subsample coefficients."))
  expect_true(all(is.finite(diag(object$cov.ssp))),
              info = paste(info_prefix, "has finite subsample variances."))
  expect_true(all(is.finite(object$coef.cmb)),
              info = paste(info_prefix, "has finite union-combined coefficients."))
  expect_true(all(is.finite(diag(object$cov.cmb))),
              info = paste(info_prefix, "has finite union-combined variances."))
}

data_rf <- make_rf_logit_data(seed = 41)
formula_rf <- Y ~ .
n.plt <- 600
n.ssp <- 700
family_rf <- "quasibinomial"
p_model <- ncol(model.matrix(formula_rf, data_rf))

## Documentation-style examples and default BL-Uni behavior.
set.seed(2)
fit_default <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  rareFeature.index = 1:2
)
check_basic_rf_fit(fit_default, p_model, "Default BL-Uni")

expect_equal(fit_default$subsample.size.expect, n.plt + n.ssp,
             info = "One-step BL-Uni should use expected sample size n.plt + n.ssp.")
expect_equal(fit_default$rareFeature.index, 2:3,
             info = "Numeric rareFeature.index should shift by one with an intercept.")
expect_null(fit_default$stage.time,
            info = "Stage timings should be omitted unless requested.")
expect_true(is.null(fit_default$coef.cmb.samples) &&
              is.null(fit_default$coef.cmb.estimators),
            info = "Only union combination output should be returned.")

summary_output <- capture.output(summary(fit_default))
expect_true(any(grepl("Model Summary", summary_output)),
            info = "summary.ssp.glm.rF should print a model summary.")
expect_true(all(c("Pilot Coefficients:",
                  "Second-Step Coefficients:",
                  "Combined-Union Coefficients:") %in% summary_output),
            info = "summary.ssp.glm.rF should print all coefficient blocks.")

## Character rare-feature names should work without intercept arithmetic.
set.seed(3)
fit_names <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "BL-Uni",
  rareFeature.index = c("Z1", "Z2")
)
expect_equal(fit_names$rareFeature.index, 2:3,
             info = "Character rare-feature names should map to design-matrix columns.")

## All criteria should run. Two-step methods should always compute union combine.
two_step_criteria <- c("Lopt", "Aopt", "R-Lopt", "BL-Lopt")
for (crit in two_step_criteria) {
  set.seed(match(crit, two_step_criteria) + 10)
  fit_crit <- ssp.glm.rF(
    formula_rf,
    data = data_rf,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = crit,
    balance.X.plt = TRUE,
    rareFeature.index = c("Z1", "Z2")
  )
  check_basic_rf_fit(fit_crit, p_model, paste("Criterion", crit))
  expect_true(identical(fit_crit$index.cmb, fit_crit$index.cmb.union),
              info = paste("Criterion", crit, "should use union combine."))
  expect_true(length(fit_crit$index.cmb.union) >= length(fit_crit$index.ssp),
              info = paste("Criterion", crit, "union size should be at least second-step size."))
}

set.seed(20)
fit_uni <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "Uni",
  rareFeature.index = c("Z1", "Z2")
)
expect_equal(fit_uni$subsample.size.expect, n.plt + n.ssp,
             info = "Uni should use expected sample size n.plt + n.ssp.")

## Automatic rare-feature detection should find the two rare binary columns.
set.seed(21)
fit_auto <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "BL-Uni",
  rareFeature.index = NULL,
  rareThreshold = 0.09
)
expect_true(all(c(2L, 3L) %in% fit_auto$rareFeature.index),
            info = "Automatic detection should include rare binary columns.")

## Response-balancing modes.
set.seed(22)
fit_y_all <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "BL-Uni",
  balance.Y.all = TRUE,
  rareFeature.index = c("Z1", "Z2")
)
expect_true(all(which(data_rf$Y == 1) %in% fit_y_all$index.ssp),
            info = "balance.Y.all should include all Y = 1 observations.")

set.seed(23)
fit_y_ssp <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "Uni",
  balance.Y.ssp = TRUE,
  rareFeature.index = c("Z1", "Z2")
)
expect_true(fit_y_ssp$Y.count.ssp > 0 &&
              fit_y_ssp$Y.count.ssp < length(fit_y_ssp$index.ssp),
            info = "balance.Y.ssp should sample both response classes.")

expect_warning(
  ssp.glm.rF(
    formula_rf,
    data = data_rf,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = "R-Lopt",
    balance.Y.ssp = TRUE,
    balance.X.plt = TRUE,
    rareFeature.index = c("Z1", "Z2")
  ),
  pattern = "balance.Y.ssp",
  info = "balance.Y.ssp should warn and be ignored for two-step criteria."
)

## Objective weighting restrictions and valid unweighted one-step fitting.
set.seed(24)
fit_unweighted <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "BL-Uni",
  objective.weight = "unweighted",
  rareFeature.index = c("Z1", "Z2")
)
check_basic_rf_fit(fit_unweighted, p_model, "Unweighted one-step")

expect_error(
  ssp.glm.rF(
    formula_rf,
    data = data_rf,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = "R-Lopt",
    objective.weight = "unweighted",
    balance.X.plt = TRUE,
    rareFeature.index = c("Z1", "Z2")
  ),
  pattern = "objective.weight",
  info = "Two-step criteria should require weighted second-step objective."
)

expect_error(
  ssp.glm.rF(
    formula_rf,
    data = data_rf,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = "BL-Uni",
    objective.weight = "unweighted",
    balance.Y.ssp = TRUE,
    rareFeature.index = c("Z1", "Z2")
  ),
  pattern = "objective.weight",
  info = "Unweighted one-step objective should error when probabilities depend on Y."
)

expect_error(
  ssp.glm.rF(
    formula_rf,
    data = data_rf,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = "R-Lopt",
    objective.weight.plt = "unweighted",
    balance.Y.plt = TRUE,
    rareFeature.index = c("Z1", "Z2")
  ),
  pattern = "objective.weight.plt",
  info = "Unweighted pilot objective should error when pilot sampling depends on Y."
)

## Control options and stage timing.
set.seed(25)
fit_estimated <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "R-Lopt",
  balance.X.plt = TRUE,
  rareFeature.index = c("Z1", "Z2"),
  control = list(alpha = 0.1, poi.method = "exact", b = 2),
  record.stage.time = TRUE
)
expect_true(!is.null(fit_estimated$stage.time) &&
              all(c("balance_score", "pilot.estimate", "subsampling",
                    "subsample.estimate", "combining.union") %in%
                    names(fit_estimated$stage.time)),
            info = "record.stage.time should return major stage timings.")

## subset, contrasts, and ... passthrough to glm.fit().
set.seed(26)
fit_subset <- ssp.glm.rF(
  formula_rf,
  data = data_rf,
  subset = 1:1200,
  n.plt = 120,
  n.ssp = 180,
  family = family_rf,
  criterion = "BL-Uni",
  rareFeature.index = 1:2
)
expect_equal(fit_subset$N, 1200,
             info = "Returned N should reflect the subsetted data.")

set.seed(27)
data_factor <- data_rf
data_factor$F1 <- sample(c("A", "B", "C"), nrow(data_factor), replace = TRUE)
fit_contrast <- ssp.glm.rF(
  Y ~ .,
  data = data_factor,
  n.plt = n.plt,
  n.ssp = n.ssp,
  family = family_rf,
  criterion = "BL-Uni",
  contrasts = list(F1 = "contr.sum"),
  rareFeature.index = c("Z1", "Z2"),
  maxit = 50
)
expect_true(all(is.finite(fit_contrast$coef.ssp)),
            info = "Contrast fit should return finite coefficients.")

## Input validation for rare features.
expect_error(
  ssp.glm.rF(
    formula_rf,
    data = data_rf,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = "BL-Uni",
    rareFeature.index = "X1"
  ),
  pattern = "not binary",
  info = "Supplying a continuous column as a rare feature should error."
)

data_nonrare <- data_rf
set.seed(28)
data_nonrare$B_hi <- rbinom(nrow(data_nonrare), 1, 0.3)
expect_warning(
  ssp.glm.rF(
    Y ~ B_hi + Z1 + Z2 + X1 + X2,
    data = data_nonrare,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = "BL-Uni",
    rareFeature.index = "B_hi"
  ),
  pattern = "prevalence",
  info = "High-prevalence binary columns supplied as rare features should warn."
)

expect_error(
  ssp.glm.rF(
    formula_rf,
    data = data_rf,
    n.plt = n.plt,
    n.ssp = n.ssp,
    family = family_rf,
    criterion = "BL-Uni",
    rareFeature.index = "not_a_column"
  ),
  pattern = "unknown rare feature",
  info = "Unknown rare-feature names should error."
)

## Non-binomial GLMs.
set.seed(29)
N_g <- 1200
Z_g <- rbinom(N_g, 1, 0.05)
X_g <- rnorm(N_g)
y_g <- 1 + 0.5 * Z_g + 0.2 * X_g + rnorm(N_g)
data_g <- data.frame(y = y_g, Z = Z_g, X = X_g)
fit_gaussian <- ssp.glm.rF(
  y ~ .,
  data = data_g,
  n.plt = 120,
  n.ssp = 180,
  family = "gaussian",
  criterion = "BL-Uni",
  rareFeature.index = "Z"
)
expect_equal(fit_gaussian$family.name, "gaussian",
             info = "Gaussian fit should record family name.")
check_basic_rf_fit(fit_gaussian, ncol(model.matrix(y ~ ., data_g)), "Gaussian BL-Uni")
gaussian_summary <- capture.output(summary(fit_gaussian))
expect_false(any(grepl("Y_sum|Y_mean", gaussian_summary)),
             info = "Non-binomial summary should not print response sums or means.")

set.seed(30)
N_p <- 1200
Z_p <- rbinom(N_p, 1, 0.05)
X_p <- rnorm(N_p)
mu_p <- exp(0.1 + 0.4 * Z_p + 0.1 * X_p)
y_p <- rpois(N_p, mu_p)
data_p <- data.frame(y = y_p, Z = Z_p, X = X_p)
expect_warning(
  ssp.glm.rF(
    y ~ .,
    data = data_p,
    n.plt = 120,
    n.ssp = 180,
    family = "poisson",
    criterion = "BL-Uni",
    balance.Y.all = TRUE,
    rareFeature.index = "Z"
  ),
  pattern = "only consider balancing Y",
  info = "Y balancing should warn and be ignored for non-binomial families."
)
fit_poisson <- suppressWarnings(
  ssp.glm.rF(
    y ~ .,
    data = data_p,
    n.plt = 120,
    n.ssp = 180,
    family = "poisson",
    criterion = "BL-Uni",
    balance.Y.all = TRUE,
    rareFeature.index = "Z"
  )
)
expect_equal(fit_poisson$family.name, "poisson",
             info = "Poisson fit should record family name.")
check_basic_rf_fit(fit_poisson, ncol(model.matrix(y ~ ., data_p)), "Poisson BL-Uni")

## Cleanup
rm(list = ls())
gc()
