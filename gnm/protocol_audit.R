# Audit of the FB2 HPO-prior simulation protocol on the GNM classifier.
#
#   Rscript protocol_audit.R [n_pcs=80] [reps=20]
#
# Reuses the out-of-fold HDRDA posteriors and applies one HPO term per face
# under different simulation protocols, to separate what the published
# numbers measured from what a clinic would see:
#   prev1_exact         every syndromic face gets one of its syndrome's terms
#                       (prevalence 1); favoured set = classes annotated with
#                       exactly that term (FB2 loocv_HPO_sim.R)
#   prev1_truepath      as above, favoured set by ontology (any class with the
#                       term or a more specific one)
#   prev_inverted_exact term applied with probability 1 - frequency
#                       (prevalence_simulation_job.R as written)
#   prev_exact          term applied with probability = frequency
#   prev_truepath       as above, ontology-aware favoured sets
#   noise_<e>           prev_truepath, plus every face (unaffected included)
#                       gets a wrong term (annotated to another class, not its
#                       own) with probability e
# Metrics: per-face (micro) and by-syndrome mean (macro, as in chapter 4,
# over syndromes and the non-syndromic class).

args <- commandArgs(trailingOnly = TRUE)
n_pcs <- if (length(args) >= 1) as.integer(args[1]) else 80L
reps <- if (length(args) >= 2) as.integer(args[2]) else 20L
here <- "/Users/jovid/Documents/Hallgrimsson/gnm_classifier"
set.seed(20261009)

cv <- readRDS(file.path(here, sprintf("cv_posteriors_%dpc.rds", n_pcs)))
P0 <- cv$post; TP <- cv$train_prior
cls <- colnames(P0)
truth <- cv$meta$Syndrome
clin <- jsonlite::fromJSON(file.path(here, "class_clinical_terms.json"), simplifyVector = FALSE)
pool <- lapply(clin, function(x) unlist(x$pool_all))
covered <- lapply(clin, function(x) unlist(x$covered))
all_terms <- unique(unlist(lapply(pool, names)))

# favoured-set membership matrices (classes x terms)
exact <- sapply(all_terms, function(t) vapply(cls, function(k) k == "Non-syndromic" || (!is.null(pool[[k]]) && t %in% names(pool[[k]])), logical(1)))
truep <- sapply(all_terms, function(t) vapply(cls, function(k) k == "Non-syndromic" || (!is.null(covered[[k]]) && t %in% covered[[k]]), logical(1)))

prior_from <- function(terms, member, tp) {
  terms <- intersect(terms, colnames(member))
  if (!length(terms)) return(tp)
  p <- rep(1, length(cls))
  used <- FALSE
  for (t in terms) {
    a <- member[, t]
    if (all(a) || !any(a)) next
    p <- p * ifelse(a, 0.9 / sum(a), 0.1 / sum(!a))
    used <- TRUE
  }
  if (!used) tp else p / sum(p)
}

rank_of <- function(p, y) match(y, names(sort(p, decreasing = TRUE)))
wrong_term <- function(y) {
  others <- setdiff(names(pool)[lengths(pool) > 0], y)
  k <- sample(others, 1)
  cand <- setdiff(names(pool[[k]]), c(names(pool[[y]]), covered[[y]]))
  if (!length(cand)) return(character(0))
  sample(cand, 1)
}

protocols <- c("shape", "prev1_exact", "prev1_truepath", "prev_inverted_exact", "prev_exact", "prev_truepath",
               "noise_0.1", "noise_0.25", "noise_0.5")
ranks <- array(NA_integer_, c(nrow(P0), length(protocols), reps), dimnames = list(NULL, protocols, NULL))
tops <- array(NA_character_, c(nrow(P0), length(protocols), reps), dimnames = list(NULL, protocols, NULL))

for (r in seq_len(reps)) {
  for (i in seq_len(nrow(P0))) {
    y <- truth[i]; p0 <- P0[i, ]; tp <- TP[i, ]
    pl <- if (y == "Non-syndromic") NULL else pool[[y]]
    t <- if (length(pl)) sample(names(pl), 1) else character(0)
    f <- if (length(t)) pl[[t]] else 0
    on <- runif(1) < f
    inv <- runif(1) < 1 - f
    post <- list(
      shape = p0,
      prev1_exact = p0 * prior_from(t, exact, tp) / tp,
      prev1_truepath = p0 * prior_from(t, truep, tp) / tp,
      prev_inverted_exact = p0 * prior_from(if (inv) t else character(0), exact, tp) / tp,
      prev_exact = p0 * prior_from(if (on) t else character(0), exact, tp) / tp,
      prev_truepath = p0 * prior_from(if (on) t else character(0), truep, tp) / tp)
    for (e in c(0.1, 0.25, 0.5)) {
      terms <- c(if (on) t, if (runif(1) < e) wrong_term(if (y == "Non-syndromic") "Non-syndromic" else y))
      post[[paste0("noise_", e)]] <- p0 * prior_from(terms, truep, tp) / tp
    }
    for (a in protocols) {
      p <- post[[a]] / sum(post[[a]])
      ranks[i, a, r] <- rank_of(p, y)
      tops[i, a, r] <- names(which.max(p))
    }
  }
  message("replicate ", r)
}

syn <- truth != "Non-syndromic"
micro <- function(k) apply(ranks[syn, , , drop = FALSE] <= k, 2, mean)
macro <- function(k) apply(sapply(split(seq_along(truth), truth), function(ix) apply(ranks[ix, , , drop = FALSE] <= k, 2, mean)), 1, mean)
out <- rbind(
  micro_top1 = micro(1), micro_top3 = micro(3), micro_top10 = micro(10),
  macro_top1 = macro(1), macro_top3 = macro(3), macro_top10 = macro(10),
  nonsyndromic_correct = apply(tops[!syn, , , drop = FALSE] == "Non-syndromic", 2, mean),
  syndromic_called_nonsyndromic = apply(tops[syn, , , drop = FALSE] == "Non-syndromic", 2, mean))
write.csv(round(out, 3), file.path(here, sprintf("protocol_audit_%dpc.csv", n_pcs)))
print(round(t(out), 3))
