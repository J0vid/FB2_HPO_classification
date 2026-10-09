# Evaluate HPO-informed priors on the cross-validated HDRDA posteriors.
#
#   Rscript evaluate_priors.R [n_pcs=80]
#
# Re-weights each face's out-of-fold posterior (no refitting):
#   shape              HDRDA alone (training class proportions as prior)
#   fb2_measured       FB2 90/10 rule with the face's own measured present terms
#   calibrated         calibrated measurement likelihood of every measured call
#   fb2_clinical_k     FB2 90/10 rule with k clinical (non-facial) terms sampled
#                      from the true syndrome's annotations (each kept with its
#                      annotated frequency), as in the FB2 dissertation
#   fb2_clinical_k+cal both together
# FB2 rule per term: classes annotated with it (true path) plus Non-syndromic
# share 90% of the prior; the rest share 10%. Unaffected faces get no
# clinical terms (nothing to report).

args <- commandArgs(trailingOnly = TRUE)
n_pcs <- if (length(args) >= 1) as.integer(args[1]) else 80L
here <- "/Users/jovid/Documents/Hallgrimsson/gnm_classifier"
bench <- "/Users/jovid/Documents/Hallgrimsson/fb_cohort_export/benchmark"
set.seed(20261009)

cv <- readRDS(file.path(here, sprintf("cv_posteriors_%dpc.rds", n_pcs)))
P0 <- cv$post; TP <- cv$train_prior; CALLS <- cv$calls; Fmat <- cv$Fmat
cls <- colnames(P0)
truth <- cv$meta$Syndrome
clin <- jsonlite::fromJSON(file.path(here, "class_clinical_terms.json"), simplifyVector = FALSE)
covered <- lapply(clin, function(x) unlist(x$covered))
fm <- jsonlite::fromJSON(file.path(bench, "full", "face_model.json"), simplifyVector = FALSE)

# Returns the prior to use, given the training prior `tp`. With no
# informative term the model's own prior is kept (as in the FB2 app).
fb2_prior <- function(terms_present, annotated_for, tp) {
  p <- rep(1, length(cls))
  used <- FALSE
  for (t in terms_present) {
    a <- vapply(cls, function(k) k == "Non-syndromic" || isTRUE(annotated_for(k, t)), logical(1))
    if (all(a) || !any(a)) next
    p <- p * ifelse(a, 0.9 / sum(a), 0.1 / sum(!a))
    used <- TRUE
  }
  if (!used) return(tp)
  p / sum(p)
}
measured_annotated <- function(k, t) Fmat[k, t] > 0
clinical_annotated <- function(k, t) !is.null(covered[[k]]) && t %in% covered[[k]]

calibrated_lik <- function(call_row, syndrome) {
  m <- if (!is.null(fm$loso[[syndrome]])) fm$loso[[syndrome]] else fm$global
  logp <- rep(0, length(cls))
  for (t in names(call_row)[!is.na(call_row)]) {
    pr <- unlist(m[[call_row[[t]]]])
    logp <- logp + log(Fmat[, t] * pr[2] + (1 - Fmat[, t]) * pr[1])
  }
  exp(logp - max(logp))
}

sample_clinical <- function(syndrome, k) {
  pool <- clin[[syndrome]]$pool
  if (is.null(pool) || !length(pool)) return(character(0))
  kept <- names(pool)[runif(length(pool)) < unlist(pool)]
  if (length(kept) > k) kept <- sample(kept, k)
  kept
}

norm <- function(p) p / sum(p)
rank_of <- function(p, y) match(y, names(sort(p, decreasing = TRUE)))
arms <- c("shape", "fb2_measured", "calibrated", "fb2_clinical_1", "fb2_clinical_3", "fb2_clinical_1+cal", "fb2_clinical_3+cal")
R <- matrix(NA_integer_, nrow(P0), length(arms), dimnames = list(NULL, arms))
TOP <- matrix(NA_character_, nrow(P0), length(arms), dimnames = list(NULL, arms))
for (i in seq_len(nrow(P0))) {
  p0 <- P0[i, ]; tp <- TP[i, ]; y <- truth[i]
  present <- names(CALLS[i, ])[which(CALLS[i, ] == "present")]
  cal <- calibrated_lik(CALLS[i, ], y)
  c1 <- if (y == "Non-syndromic") character(0) else sample_clinical(y, 1)
  c3 <- if (y == "Non-syndromic") character(0) else sample_clinical(y, 3)
  pc1 <- fb2_prior(c1, clinical_annotated, tp) / tp
  pc3 <- fb2_prior(c3, clinical_annotated, tp) / tp
  ps <- list(
    shape = p0,
    fb2_measured = p0 * fb2_prior(present, measured_annotated, tp) / tp,
    calibrated = p0 * cal,
    fb2_clinical_1 = p0 * pc1,
    fb2_clinical_3 = p0 * pc3,
    `fb2_clinical_1+cal` = p0 * pc1 * cal,
    `fb2_clinical_3+cal` = p0 * pc3 * cal)
  for (a in arms) {
    p <- norm(ps[[a]])
    R[i, a] <- rank_of(p, y)
    TOP[i, a] <- names(which.max(p))
  }
}

summ <- function(r) round(c(top1 = mean(r <= 1), top3 = mean(r <= 3), top10 = mean(r <= 10), median = median(r)), 3)
syn <- truth != "Non-syndromic"
out <- list(
  syndromic = apply(R[syn, ], 2, summ),
  nonsyndromic_correct = round(colMeans(TOP[!syn, ] == "Non-syndromic"), 3),
  syndromic_called_nonsyndromic = round(colMeans(TOP[syn, ] == "Non-syndromic"), 3)
)
per <- do.call(rbind, lapply(split(seq_along(truth), truth), function(ix) data.frame(
  syndrome = truth[ix[1]], n = length(ix), t(colMeans(R[ix, , drop = FALSE] <= 3)), check.names = FALSE)))
write.csv(per[order(-per$n), ], file.path(here, sprintf("priors_per_class_top3_%dpc.csv", n_pcs)), row.names = FALSE)
writeLines(jsonlite::toJSON(out, digits = 3, pretty = TRUE), file.path(here, sprintf("priors_summary_%dpc.json", n_pcs)))
print(out)
