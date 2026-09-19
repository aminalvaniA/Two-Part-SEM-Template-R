################################################################################
# TEMPLATE: Two-Part (Hurdle) Semi-Continuous Moderated-Mediation SEM
# for zero-inflated continuous/count outcomes
#
# Use case: any outcome where a meaningful subgroup of respondents scores
# exactly zero (no symptoms/behavior) and, among those who score above
# zero, severity is continuous and skewed. Examples: NSSI frequency,
# problem-gambling symptoms, screen-time overuse, dissociation episodes,
# binge-eating frequency, panic-attack counts.
#
# WHY NOT ORDINARY REGRESSION?
# A single linear or Poisson-type model on a zero-inflated outcome
# conflates two distinct questions: (a) does the behavior/symptom occur
# at all, and (b) how severe is it once present. Forcing one model to
# answer both usually misspecifies the outcome distribution and biases
# coefficients. This template splits the outcome into:
#   Part A (entry):    probit SEM, WLSMV estimator, "any symptoms?" (0/1)
#   Part B (severity):  continuous SEM, MLR estimator, log(Y) | Y > 0
# and layers on continuous moderators for the mediator -> outcome paths.
#
# HOW TO USE THIS TEMPLATE
# 1. Fill in the CONFIGURATION block below with your own variable names.
# 2. Point data_path at your .sav file (or adapt read_data() for .csv).
# 3. Run the whole script. It prints structural paths, conditional
#    indirect effects, indices of moderated mediation, Johnson-Neyman
#    regions, and writes diagnostic figures to fig_dir.
# 4. Run once against toy_data.csv first to confirm your setup works
#    before pointing it at real data.
################################################################################

# ==============================================================================
# 0. CONFIGURATION -- edit this block for your own study
# ==============================================================================

# Predictor (X): the distal cause feeding into the mediators
X_var <- "trauma_or_predictor"

# Mediator(s) (M): one or more parallel mediators between X and Y.
# The template supports 1-3 mediators out of the box; add more by
# extending build_structural_model() below.
M_vars <- c("attachment_anx", "attachment_avo")

# Outcome (Y): the zero-inflated variable. Must be numeric, >= 0.
# 0 = no symptoms/behavior; > 0 = count or severity among those affected.
Y_zero_inflated <- "NSSI_frequency_or_screen_time"

# Moderator(s) (W): continuous variable(s) moderating each M -> Y path.
# The template supports 1-2 moderators out of the box.
W_moderators <- c("dissociation_or_alexithymia")

# Upper bound of the outcome's plausible range, if known (used only for a
# sanity check on out-of-range values). Set to Inf if there is no fixed
# maximum (e.g., open-ended counts).
Y_max_plausible <- Inf

# Path to the data file (.sav via haven, or .csv -- see read_data()).
data_path <- "toy_data.csv"

# Output directory for figures.
fig_dir <- "figures"

# ==============================================================================
# PACKAGES
# ==============================================================================
required_packages <- c("haven", "psych", "lavaan", "ggplot2", "pscl")
for (p in required_packages) {
  if (!requireNamespace(p, quietly = TRUE)) stop(sprintf("Package '%s' is required. Run install.packages('%s').", p, p))
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

set.seed(12345)
ARROW <- "->"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================
read_data <- function(path) {
  if (!file.exists(path)) stop(sprintf("Data file not found at:\n  %s", path))
  ext <- tolower(tools::file_ext(path))
  if (ext == "sav") {
    as.data.frame(haven::read_sav(path))
  } else if (ext == "csv") {
    utils::read.csv(path, stringsAsFactors = FALSE)
  } else {
    stop("Unsupported file type. Use .sav or .csv, or extend read_data().")
  }
}

raw <- read_data(data_path)

analysis_vars <- c(X_var, M_vars, Y_zero_inflated, W_moderators)
missing_vars <- setdiff(analysis_vars, names(raw))
if (length(missing_vars) > 0) {
  stop(sprintf("Missing required variables in data: %s", paste(missing_vars, collapse = ", ")))
}

df <- raw[, analysis_vars]
df[analysis_vars] <- lapply(df[analysis_vars], as.numeric)

range_bad <- !is.na(df[[Y_zero_inflated]]) & (df[[Y_zero_inflated]] < 0 | df[[Y_zero_inflated]] > Y_max_plausible)
if (any(range_bad)) {
  warning(sprintf("%d case(s) had %s outside [0, %s] and were removed.", sum(range_bad), Y_zero_inflated, Y_max_plausible))
  df <- df[!range_bad, ]
}

n_before <- nrow(df)
df <- df[stats::complete.cases(df), ]
n_after <- nrow(df)

# ==============================================================================
# 2. DERIVED VARIABLES
# ==============================================================================
Y_bin_col <- paste0(Y_zero_inflated, "_bin")
Y_bin_ord_col <- paste0(Y_zero_inflated, "_bin_ord")
df[[Y_bin_col]] <- ifelse(df[[Y_zero_inflated]] > 0, 1, 0)
df[[Y_bin_ord_col]] <- factor(df[[Y_bin_col]], levels = c(0, 1), ordered = TRUE)

# Grand-mean center all mediators and moderators (required before forming
# interaction terms, to reduce multicollinearity and ease interpretation).
center_var <- function(x) x - mean(x, na.rm = TRUE)
M_c_vars <- paste0(M_vars, "_c")
W_c_vars <- paste0(W_moderators, "_c")
for (i in seq_along(M_vars)) df[[M_c_vars[i]]] <- center_var(df[[M_vars[i]]])
for (i in seq_along(W_moderators)) df[[W_c_vars[i]]] <- center_var(df[[W_moderators[i]]])

# Interaction terms: every mediator x every moderator
interaction_grid <- expand.grid(m = seq_along(M_c_vars), w = seq_along(W_c_vars))
int_names <- mapply(function(m, w) paste0("int_M", m, "_W", w), interaction_grid$m, interaction_grid$w)
for (k in seq_len(nrow(interaction_grid))) {
  m_i <- interaction_grid$m[k]; w_i <- interaction_grid$w[k]
  df[[int_names[k]]] <- df[[M_c_vars[m_i]]] * df[[W_c_vars[w_i]]]
}

zero_diag <- function(x) c(n = length(x), pct_zero = mean(x == 0) * 100, mean = mean(x), sd = sd(x),
                            skew = psych::skew(x), kurtosis = psych::kurtosi(x))
zdiag_full <- zero_diag(df[[Y_zero_inflated]])

cat("\n==================== SAMPLE ====================\n")
cat(sprintf("N analyzed = %d (%d cases removed for missing/out-of-range data)\n", n_after, n_before - n_after))
cat(sprintf("%s: %.1f%% zeros, M = %.2f, SD = %.2f, skew = %.2f, kurtosis = %.2f\n",
            Y_zero_inflated, zdiag_full["pct_zero"], zdiag_full["mean"], zdiag_full["sd"], zdiag_full["skew"], zdiag_full["kurtosis"]))

# ==============================================================================
# 3. JOHNSON-NEYMAN HELPERS
# Region of significance for a continuous moderator, via the delta method
# on the model-implied covariance matrix (Johnson & Neyman, 1936;
# Bauer & Curran, 2005).
# ==============================================================================
johnson_neyman <- function(b, bw, var_b, var_bw, cov_b_bw, w_range, z = 1.959964) {
  A <- bw^2 - z^2 * var_bw
  B <- 2 * (b * bw - z^2 * cov_b_bw)
  C <- b^2 - z^2 * var_b
  if (abs(A) < 1e-12) {
    if (abs(B) < 1e-12) return(numeric(0))
    roots <- -C / B
  } else {
    disc <- B^2 - 4 * A * C
    if (disc < 0) return(numeric(0))
    roots <- c((-B + sqrt(disc)) / (2 * A), (-B - sqrt(disc)) / (2 * A))
  }
  roots <- sort(roots[is.finite(roots)])
  roots[roots >= w_range[1] & roots <= w_range[2]]
}

cond_effect_curve <- function(b, bw, var_b, var_bw, cov_b_bw, w_seq, z = 1.959964) {
  theta <- b + bw * w_seq
  se    <- sqrt(var_b + (w_seq^2) * var_bw + 2 * w_seq * cov_b_bw)
  data.frame(W = w_seq, theta = theta, se = se, lo = theta - z * se, hi = theta + z * se,
             sig = (theta - z * se > 0) | (theta + z * se < 0))
}

# ==============================================================================
# 4. STRUCTURAL MODEL BUILDER
# Generalized for K mediators and up to 2 moderators. Builds:
#   - a-paths:  X -> each mediator
#   - b-paths:  each mediator (+ interactions) -> Y
#   - defined parameters: conditional indirect effects at representative
#     values of each moderator, and indices of moderated mediation
# ==============================================================================
build_structural_model <- function(y_var) {
  lines <- character(0)
  a_labels <- paste0("a", seq_along(M_c_vars))
  for (i in seq_along(M_c_vars)) {
    lines <- c(lines, sprintf("  %s ~ %s*%s\n", M_c_vars[i], a_labels[i], X_var))
  }

  b_labels <- paste0("b", seq_along(M_c_vars))
  bw_labels <- matrix(NA_character_, nrow = length(M_c_vars), ncol = length(W_c_vars))
  for (m in seq_along(M_c_vars)) for (w in seq_along(W_c_vars)) bw_labels[m, w] <- sprintf("b%dw%d", m, w)

  rhs_terms <- c(sprintf("cp*%s", X_var))
  for (i in seq_along(M_c_vars)) rhs_terms <- c(rhs_terms, sprintf("%s*%s", b_labels[i], M_c_vars[i]))
  for (w in seq_along(W_c_vars)) rhs_terms <- c(rhs_terms, sprintf("w%d*%s", w, W_c_vars[w]))
  for (k in seq_len(nrow(interaction_grid))) {
    m_i <- interaction_grid$m[k]; w_i <- interaction_grid$w[k]
    rhs_terms <- c(rhs_terms, sprintf("%s*%s", bw_labels[m_i, w_i], int_names[k]))
  }
  lines <- c(lines, sprintf("  %s ~ %s\n", y_var, paste(rhs_terms, collapse = " + ")))

  # Conditional indirect effects at the 16th/50th/84th percentile of each
  # moderator, and at mean +/- 1 SD, for every mediator x moderator pair.
  for (m in seq_along(M_c_vars)) {
    for (w in seq_along(W_c_vars)) {
      w_col <- W_c_vars[w]
      q16 <- stats::quantile(df[[w_col]], .16); q50 <- stats::quantile(df[[w_col]], .50); q84 <- stats::quantile(df[[w_col]], .84)
      sdw <- sd(df[[w_col]])
      lines <- c(lines,
        sprintf("  ind%d_W%d_lo16  := %s*(%s + %s*(%.6f))\n", m, w, a_labels[m], b_labels[m], bw_labels[m, w], q16),
        sprintf("  ind%d_W%d_med50 := %s*(%s + %s*(%.6f))\n", m, w, a_labels[m], b_labels[m], bw_labels[m, w], q50),
        sprintf("  ind%d_W%d_hi84  := %s*(%s + %s*(%.6f))\n", m, w, a_labels[m], b_labels[m], bw_labels[m, w], q84),
        sprintf("  ind%d_W%d_loSD  := %s*(%s + %s*(%.6f))\n", m, w, a_labels[m], b_labels[m], bw_labels[m, w], -sdw),
        sprintf("  ind%d_W%d_mean  := %s*%s\n", m, w, a_labels[m], b_labels[m]),
        sprintf("  ind%d_W%d_hiSD  := %s*(%s + %s*(%.6f))\n", m, w, a_labels[m], b_labels[m], bw_labels[m, w], sdw),
        sprintf("  imm%d_W%d := %s*%s\n", m, w, a_labels[m], bw_labels[m, w])
      )
    }
  }
  paste0(lines, collapse = "")
}

# ==============================================================================
# 5. FIT MODELS
# Part A: probit latent-response SEM (WLSMV) -- transition to any symptoms
# Part B: continuous SEM (MLR) -- severity among those affected
# ==============================================================================
model_A <- build_structural_model(Y_bin_ord_col)
fit_A <- tryCatch(
  lavaan::sem(model_A, data = df, ordered = Y_bin_ord_col, estimator = "WLSMV"),
  error = function(e) { message("Part A SEM failed: ", conditionMessage(e)); NULL }
)

df_pos <- df[df[[Y_zero_inflated]] > 0, ]
log_y_col <- paste0("log_", Y_zero_inflated)
df_pos[[log_y_col]] <- log(df_pos[[Y_zero_inflated]])
model_B <- build_structural_model(log_y_col)
fit_B <- tryCatch(
  lavaan::sem(model_B, data = df_pos, estimator = "MLR"),
  error = function(e) { message("Part B SEM failed: ", conditionMessage(e)); NULL }
)

# Robustness companion to Part B: Gamma(log) GLM on raw (untransformed) severity
gamma_rhs <- paste(c(X_var, M_c_vars, W_c_vars, int_names), collapse = " + ")
frm_gamma <- stats::as.formula(paste(Y_zero_inflated, "~", gamma_rhs))
fit_B_gamma <- tryCatch(stats::glm(frm_gamma, data = df_pos, family = Gamma(link = "log")),
                         error = function(e) { message("Gamma(log) model failed: ", conditionMessage(e)); NULL })

extract_param_tab <- function(fit) {
  if (is.null(fit)) return(NULL)
  pe <- lavaan::parameterEstimates(fit, standardized = FALSE)
  reg <- pe[pe$op == "~", ]
  data.frame(Path = paste0(reg$lhs, " ", ARROW, " ", reg$rhs, ifelse(reg$label != "", paste0(" (", reg$label, ")"), "")),
             b = sprintf("%.4f", reg$est), SE = sprintf("%.4f", reg$se), z = sprintf("%.3f", reg$z),
             p = ifelse(reg$pvalue < .001, "< .001", sprintf("%.3f", reg$pvalue)), stringsAsFactors = FALSE)
}

extract_defined_tab <- function(fit) {
  if (is.null(fit)) return(NULL)
  pe <- lavaan::parameterEstimates(fit, standardized = FALSE)
  ind <- pe[pe$op == ":=", ]
  data.frame(Effect = ind$label, Estimate = sprintf("%.4f", ind$est), SE = sprintf("%.4f", ind$se),
             LLCI = sprintf("%.4f", ind$ci.lower), ULCI = sprintf("%.4f", ind$ci.upper),
             Significant = ifelse(ind$ci.lower > 0 | ind$ci.upper < 0, "Yes", "No"), stringsAsFactors = FALSE)
}

cat("\n==================== PART A: PROBIT SEM (WLSMV) ====================\n")
cat(sprintf("Converged: %s | N = %d\n", if (!is.null(fit_A)) lavaan::lavInspect(fit_A, "converged") else NA, nrow(df)))
cat("\n-- Structural paths --\n"); if (!is.null(fit_A)) print(extract_param_tab(fit_A), row.names = FALSE)
cat("\n-- Conditional indirect effects / indices of moderated mediation --\n"); if (!is.null(fit_A)) print(extract_defined_tab(fit_A), row.names = FALSE)

cat("\n==================== PART B: CONTINUOUS SEM (MLR, log severity) ====================\n")
cat(sprintf("Converged: %s | n with %s > 0 = %d of %d (%.1f%%)\n",
            if (!is.null(fit_B)) lavaan::lavInspect(fit_B, "converged") else NA,
            Y_zero_inflated, nrow(df_pos), nrow(df), 100 * nrow(df_pos) / nrow(df)))
cat("\n-- Structural paths --\n"); if (!is.null(fit_B)) print(extract_param_tab(fit_B), row.names = FALSE)
cat("\n-- Conditional indirect effects / indices of moderated mediation --\n"); if (!is.null(fit_B)) print(extract_defined_tab(fit_B), row.names = FALSE)

cat("\n-- Robustness check: Gamma(log) GLM on raw severity --\n")
if (!is.null(fit_B_gamma)) print(summary(fit_B_gamma)$coefficients) else cat("Model did not converge.\n")

# ==============================================================================
# 6. ROBUSTNESS CHECK: ZERO-INFLATED NEGATIVE BINOMIAL (ZINB)
# Cross-checks the pattern of effects under a competing zero-generating
# process (a structurally-always-zero latent class) against the hurdle
# mechanism used above.
# ==============================================================================
zinb_rhs <- paste(c(X_var, M_c_vars, W_c_vars, int_names), collapse = " + ")
zinb_frm <- stats::as.formula(paste(Y_zero_inflated, "~", zinb_rhs, "|", zinb_rhs))
fit_zinb <- tryCatch(pscl::zeroinfl(zinb_frm, data = df, dist = "negbin", link = "logit"),
                      error = function(e) { message("ZINB model failed: ", conditionMessage(e)); NULL })

cat("\n==================== ROBUSTNESS CHECK: ZINB ====================\n")
if (!is.null(fit_zinb)) {
  zs <- summary(fit_zinb)
  cat(sprintf("Log-likelihood = %.2f, AIC = %.2f, BIC = %.2f, theta = %.3f\n",
              as.numeric(logLik(fit_zinb)), AIC(fit_zinb), stats::BIC(fit_zinb), fit_zinb$theta))
  cat("\n-- Count part (log link) --\n"); print(zs$coefficients$count)
  cat("\n-- Zero-inflation part (logit link) --\n"); print(zs$coefficients$zero)
} else {
  cat("Model did not converge (expected with very small samples, e.g. the toy dataset).\n")
}

# ==============================================================================
# 7. JOHNSON-NEYMAN REGIONS OF SIGNIFICANCE (first moderator, first mediator)
# Extend the loop below to cover every mediator x moderator pair you need.
# ==============================================================================
jn_for_fit <- function(fit, m_idx = 1, w_idx = 1) {
  if (is.null(fit)) return(NULL)
  vc <- lavaan::vcov(fit)
  pe <- lavaan::parameterEstimates(fit)
  get_est <- function(lab) pe$est[pe$label == lab][1]
  w_col <- W_c_vars[w_idx]
  w_range <- range(df[[w_col]])
  w_seq   <- seq(w_range[1], w_range[2], length.out = 200)
  b_lab <- sprintf("b%d", m_idx); bw_lab <- sprintf("b%dw%d", m_idx, w_idx)
  b <- get_est(b_lab); bw <- get_est(bw_lab)
  vb <- vc[b_lab, b_lab]; vbw <- vc[bw_lab, bw_lab]; cvb <- vc[b_lab, bw_lab]
  jn_pts_c <- johnson_neyman(b, bw, vb, vbw, cvb, w_range)
  curve <- cond_effect_curve(b, bw, vb, vbw, cvb, w_seq)
  curve$W_raw <- curve$W + mean(df[[W_moderators[w_idx]]])
  list(mediator = M_vars[m_idx], moderator = W_moderators[w_idx],
       jn_points_raw = jn_pts_c + mean(df[[W_moderators[w_idx]]]), curve = curve)
}

jn_summary_line <- function(jn, part_label) {
  if (is.null(jn)) return(sprintf("%s: model did not converge.", part_label))
  if (length(jn$jn_points_raw) == 0) {
    sprintf("%s (%s moderated by %s): no boundary within the observed range.", part_label, jn$mediator, jn$moderator)
  } else {
    sprintf("%s (%s moderated by %s): boundary at %s.", part_label, jn$mediator, jn$moderator,
            paste(sprintf("%.2f", jn$jn_points_raw), collapse = ", "))
  }
}

jn_A <- jn_for_fit(fit_A, m_idx = 1, w_idx = 1)
jn_B <- jn_for_fit(fit_B, m_idx = 1, w_idx = 1)

cat("\n==================== JOHNSON-NEYMAN REGIONS ====================\n")
cat(jn_summary_line(jn_A, "Part A (entry, probit)"), "\n")
cat(jn_summary_line(jn_B, "Part B (severity, log)"), "\n\n")

# ==============================================================================
# 8. FIGURES
# ==============================================================================
theme_pub <- ggplot2::theme_minimal(base_size = 12, base_family = "serif") +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), plot.title = ggplot2::element_text(face = "bold", size = 13))

fig1_path <- file.path(fig_dir, "Figure1_ZeroInflated_Distribution.png")
ggplot2::ggsave(fig1_path, width = 8, height = 5, dpi = 300, plot =
  ggplot2::ggplot(df, ggplot2::aes(x = .data[[Y_zero_inflated]])) +
    ggplot2::geom_bar(fill = "#4C72B0", color = "white", linewidth = 0.2) +
    ggplot2::labs(title = paste("Distribution of", Y_zero_inflated),
                  subtitle = sprintf("%.1f%% of the sample scored zero (N = %d)", zdiag_full["pct_zero"], n_after),
                  x = Y_zero_inflated, y = "Frequency") + theme_pub)

make_jn_plot <- function(jn, y_lab, title_txt, path) {
  if (is.null(jn)) return(NULL)
  ggplot2::ggsave(path, width = 8, height = 5, dpi = 300, plot =
    ggplot2::ggplot(jn$curve, ggplot2::aes(x = W_raw, y = theta)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, fill = "#4C72B0") +
      ggplot2::geom_line(linewidth = 1, color = "#4C72B0") +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::labs(title = title_txt, x = paste(jn$moderator, "(raw score)"), y = y_lab) +
      theme_pub)
  path
}
make_jn_plot(jn_A, "Conditional effect on P(Y > 0) [probit]", "Johnson-Neyman: Part A (Entry)", file.path(fig_dir, "Figure2a_JohnsonNeyman_PartA.png"))
make_jn_plot(jn_B, "Conditional effect on log(Y) | Y > 0", "Johnson-Neyman: Part B (Severity)", file.path(fig_dir, "Figure2b_JohnsonNeyman_PartB.png"))

cat("\n==================== DONE ====================\n")
cat("Figures written to:", fig_dir, "\n")
