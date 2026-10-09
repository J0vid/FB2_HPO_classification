# HDRDA syndrome classifier on the GNM topology, with HPO-informed priors.
#
#   Rscript fit_hdrda.R [n_pcs=80] [folds=10]
#
# Follows the FB2 dissertation pipeline (Procrustes shape -> PCA -> PC scores
# adjusted for Sex + poly(Age, 3) -> HDRDA, sparsediscrim::rda_high_dim), refitting PCA, the
# age/sex model and HDRDA inside each cross-validation fold.
#
# Posteriors are then re-weighted by each face's own measured HPO calls,
# without refitting (posterior_k ∝ posterior_k · prior_new_k / prior_train_k):
#   fb2:        for each present term, the classes annotated with it plus
#               Non-syndromic share 90% of the prior and the rest share 10%
#               (Aponte dissertation, FB2_HPO_classification).
#   calibrated: every measured call (present, borderline, not called,
#               excluded) multiplies the prior by its calibrated likelihood
#               F_k·P(call | present) + (1 − F_k)·P(call | absent), using the
#               FaceBase measurement model fitted leaving out the syndrome.
# Participant-level outputs stay in this directory and out of git.

args <- commandArgs(trailingOnly = TRUE)
n_pcs <- if (length(args) >= 1) as.integer(args[1]) else 80L
k_folds <- if (length(args) >= 2) as.integer(args[2]) else 10L
set.seed(20261008)
here <- "/Users/jovid/Documents/Hallgrimsson/gnm_classifier"
bench <- "/Users/jovid/Documents/Hallgrimsson/fb_cohort_export/benchmark"

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- data ---------------------------------------------------------------
meta <- read.csv(file.path(here, "gnm_meta.csv"), stringsAsFactors = FALSE)
meta$Syndrome[meta$group == "nonsyndromic"] <- "Non-syndromic"
nv <- 4503L  # measured GNM vertices (measured_vertex_ids.npy)
con <- file(file.path(here, "gnm_measured_f32.bin"), "rb")
raw <- readBin(con, "numeric", n = nrow(meta) * nv * 3, size = 4)
close(con)
A <- aperm(array(raw, dim = c(3, nv, nrow(meta))), c(2, 1, 3))  # nv x 3 x n
rm(raw); gc()

ct <- jsonlite::fromJSON(file.path(here, "class_term_freq.json"), simplifyVector = FALSE)
classes <- unlist(ct$classes)
keep <- meta$qc_status != "fail" & meta$Syndrome %in% classes & !is.na(meta$Age) & meta$Sex %in% c("F", "M")
meta <- meta[keep, ]
A <- A[, , keep]
y <- factor(meta$Syndrome, levels = c("Non-syndromic", setdiff(classes, "Non-syndromic")))
message(sprintf("%d faces, %d classes, %d measured vertices", nrow(meta), nlevels(y), nv))

# ---- Procrustes shape (label-free, so done once on all faces) ------------
cache <- file.path(here, "gpa_cache.rds")  # scratch; delete when done
if (file.exists(cache) && nrow(readRDS(cache)) == nrow(meta)) {
  X <- readRDS(cache)
} else {
  gpa <- Morpho::ProcGPA(A, scale = TRUE, CSinit = TRUE, silent = TRUE)
  X <- t(apply(gpa$rotated, 3, as.vector))  # n x (3 nv)
  saveRDS(X, cache, compress = FALSE)
}
rm(A); gc()

# ---- HPO evidence per face ---------------------------------------------
calls <- jsonlite::fromJSON(file.path(bench, "face_calls.json"), simplifyVector = FALSE)
names(calls) <- vapply(calls, function(x) as.character(x$image), "")
fm_all <- jsonlite::fromJSON(file.path(bench, "full", "face_model.json"), simplifyVector = FALSE)
terms <- unlist(ct$measured_terms)
freq <- ct$freq
cls <- levels(y)
Fmat <- matrix(NA_real_, length(cls), length(terms), dimnames = list(cls, terms))
for (k in cls[-1]) {
  f <- freq[[k]]
  if (!is.null(f)) Fmat[k, ] <- vapply(terms, function(t) if (is.null(f[[t]])) NA_real_ else f[[t]], numeric(1))
}
annotated <- rowSums(!is.na(Fmat)) > 0
Fmat[is.na(Fmat)] <- 0
Fmat["Non-syndromic", ] <- 0

face_call_vector <- function(img) {
  x <- calls[[as.character(img)]]
  out <- setNames(rep(NA_character_, length(terms)), terms)
  if (is.null(x)) return(out)
  for (t in intersect(unlist(x$measured), terms)) out[t] <- "indeterminate"
  for (t in intersect(unlist(x$excluded), terms)) out[t] <- "excluded"
  for (t in intersect(unlist(x$borderline), terms)) out[t] <- "borderline"
  for (t in intersect(unlist(x$present), terms)) out[t] <- "present"
  out
}
CALLS <- t(vapply(meta$Image_Name, face_call_vector, character(length(terms))))

prior_fb2 <- function(call_row) {
  p <- rep(1, length(cls))
  for (t in names(call_row)[which(call_row == "present")]) {
    a <- (Fmat[, t] > 0) | cls == "Non-syndromic"
    if (!any(a) || all(a)) next
    w <- ifelse(a, 0.9 / sum(a), 0.1 / sum(!a))
    p <- p * w
  }
  p / sum(p)
}

prior_calibrated <- function(call_row, syndrome) {
  m <- fm_all$loso[[syndrome]] %||% fm_all$global
  logp <- rep(0, length(cls))
  for (t in names(call_row)[!is.na(call_row)]) {
    pr <- unlist(m[[call_row[[t]]]])  # (P(call | absent), P(call | present))
    lik <- Fmat[, t] * pr[2] + (1 - Fmat[, t]) * pr[1]
    # classes without annotations keep F = 0, i.e. the Non-syndromic baseline
    logp <- logp + log(lik)
  }
  p <- exp(logp - max(logp))
  p / sum(p)
}
# ---- one fold -----------------------------------------------------------
fit_predict <- function(train, test) {
  mu <- colMeans(X[train, ])
  Xc <- sweep(X[train, ], 2, mu)
  G <- tcrossprod(Xc)
  e <- eigen(G, symmetric = TRUE)
  V <- crossprod(Xc, e$vectors[, 1:n_pcs]) %*% diag(1 / sqrt(e$values[1:n_pcs]))  # loadings
  S_tr <- Xc %*% V
  S_te <- sweep(X[test, , drop = FALSE], 2, mu) %*% V
  d_tr <- data.frame(Age = meta$Age[train], Sex = meta$Sex[train])
  d_te <- data.frame(Age = meta$Age[test], Sex = meta$Sex[test])
  adj <- lm(S_tr ~ Sex + poly(Age, 3), data = d_tr)
  R_tr <- S_tr - fitted(adj)
  R_te <- S_te - predict(adj, newdata = d_te)
  colnames(R_tr) <- colnames(R_te) <- paste0("PC", seq_len(n_pcs))
  # sparsediscrim >= 0.3 renamed hdrda() to rda_high_dim(); defaults
  # (lambda = 1, gamma = 0, ridge) match the FB2 fits.
  mod <- sparsediscrim::rda_high_dim(R_tr, y[train])
  post <- as.matrix(predict(mod, R_te, type = "prob"))[, cls, drop = FALSE]
  list(post = post, train_prior = as.numeric(table(y[train]) / length(train)),
       model = list(mu = mu, V = V, adj = adj, hdrda = mod))
}

rank_of <- function(p, truth) match(truth, names(sort(p, decreasing = TRUE)))
fold <- sample(rep_len(seq_len(k_folds), nrow(meta)))
res <- vector("list", nrow(meta))
POST <- matrix(NA_real_, nrow(meta), length(cls), dimnames = list(meta$Image_Name, cls))
TRAIN_PRIOR <- matrix(NA_real_, nrow(meta), length(cls), dimnames = list(meta$Image_Name, cls))
t0 <- Sys.time()
for (f in seq_len(k_folds)) {
  test <- which(fold == f)
  r <- fit_predict(which(fold != f), test)
  for (j in seq_along(test)) {
    i <- test[j]
    p0 <- r$post[j, ]
    POST[i, ] <- p0
    TRAIN_PRIOR[i, ] <- r$train_prior
    p_fb2 <- p0 * prior_fb2(CALLS[i, ]) / r$train_prior
    p_cal <- p0 * prior_calibrated(CALLS[i, ], as.character(y[i]))
    norm <- function(p) p / sum(p)
    res[[i]] <- data.frame(
      image = meta$Image_Name[i], syndrome = as.character(y[i]), age = meta$Age[i],
      n_present = sum(CALLS[i, ] == "present", na.rm = TRUE),
      rank_shape = rank_of(norm(p0), as.character(y[i])),
      rank_fb2 = rank_of(norm(p_fb2), as.character(y[i])),
      rank_calibrated = rank_of(norm(p_cal), as.character(y[i])),
      top_shape = names(which.max(p0)), top_fb2 = names(which.max(p_fb2)), top_calibrated = names(which.max(p_cal)),
      stringsAsFactors = FALSE)
  }
  message(sprintf("fold %d/%d done (%.1f min)", f, k_folds, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, res)
saveRDS(list(post = POST, train_prior = TRAIN_PRIOR, meta = meta[, c("Image_Name", "Syndrome", "Age", "Sex")],
             fold = fold, calls = CALLS, Fmat = Fmat, annotated = annotated),
        file.path(here, sprintf("cv_posteriors_%dpc.rds", n_pcs)))
write.csv(res, file.path(here, sprintf("cv_predictions_%dpc.csv", n_pcs)), row.names = FALSE)

# ---- summary ------------------------------------------------------------
summ <- function(r) c(n = length(r), top1 = mean(r <= 1), top3 = mean(r <= 3), top10 = mean(r <= 10), median = median(r))
syn <- res$syndrome != "Non-syndromic"
ann_cls <- names(annotated)[annotated]
out <- list(
  n_pcs = n_pcs, folds = k_folds, n = nrow(res), classes = nlevels(y),
  all = sapply(res[, c("rank_shape", "rank_fb2", "rank_calibrated")], summ),
  syndromic = sapply(res[syn, c("rank_shape", "rank_fb2", "rank_calibrated")], summ),
  syndromic_annotated = sapply(res[syn & res$syndrome %in% ann_cls, c("rank_shape", "rank_fb2", "rank_calibrated")], summ),
  nonsyndromic_correct = c(shape = mean(res$top_shape[!syn] == "Non-syndromic"),
                           fb2 = mean(res$top_fb2[!syn] == "Non-syndromic"),
                           calibrated = mean(res$top_calibrated[!syn] == "Non-syndromic")),
  syndromic_called_nonsyndromic = c(shape = mean(res$top_shape[syn] == "Non-syndromic"),
                                    fb2 = mean(res$top_fb2[syn] == "Non-syndromic"),
                                    calibrated = mean(res$top_calibrated[syn] == "Non-syndromic"))
)
per_class <- do.call(rbind, lapply(split(res, res$syndrome), function(d) data.frame(
  syndrome = d$syndrome[1], n = nrow(d), annotated = d$syndrome[1] %in% c("Non-syndromic", ann_cls),
  top1_shape = mean(d$rank_shape == 1), top1_fb2 = mean(d$rank_fb2 == 1), top1_calibrated = mean(d$rank_calibrated == 1),
  top3_shape = mean(d$rank_shape <= 3), top3_fb2 = mean(d$rank_fb2 <= 3), top3_calibrated = mean(d$rank_calibrated <= 3))))
write.csv(per_class[order(-per_class$n), ], file.path(here, sprintf("cv_per_class_%dpc.csv", n_pcs)), row.names = FALSE)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = 4, pretty = TRUE), file.path(here, sprintf("cv_summary_%dpc.json", n_pcs)))
print(out)
