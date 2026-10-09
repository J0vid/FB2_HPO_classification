# HDRDA syndrome classifier on GNM v3: cross-validated results (2026-10-09)

## Setup

- **Data:** 5,053 FaceBase faces on GNM v3. Only the 4,503 measured face vertices are used; model-filled vertices are left out.
  - QC failures are excluded.
  - There are 78 classes: syndromes with at least 10 faces, plus Non-syndromic (n = 2,365).
- **Pipeline** (as in FB2):
  - Procrustes shape, then PCA with 80 PCs;
  - PC scores adjusted for Sex + poly(Age, 3);
  - HDRDA (`sparsediscrim::rda_high_dim`, lambda = 1, gamma = 0).
  - 10-fold cross-validation, with the PCA, the age/sex model and HDRDA refitted in every fold.
- **Priors:** posteriors are re-weighted without refitting (`evaluate_priors.R`).
  - FB2 rule: for each term, the classes annotated with it plus Non-syndromic share 90% of the prior and the rest share 10%. With no informative term, the model's own prior is kept.

## Syndromic faces (n = 2,688)

| Arm | Top-1 | Top-3 | Top-10 | Syndromic called Non-syndromic |
|---|---|---|---|---|
| Shape only | 34.9% | 57.1% | 76.9% | 34.5% |
| + FB2 rule, the face's own measured present terms | 27.7% | 53.5% | 75.4% | 51.3% |
| + calibrated likelihood of all measured calls | 35.1% | 55.9% | 76.8% | 32.4% |
| + FB2 rule, 1 sampled clinical term | 45.6% | 67.9% | 83.1% | 26.1% |
| + FB2 rule, 3 sampled clinical terms | **51.6%** | **80.6%** | **89.8%** | 35.3% |
| + 3 clinical terms + calibrated face | 51.5% | 80.5% | 90.0% | 34.3% |

Unaffected faces are correctly called Non-syndromic in 96–97% of cases in every arm.

Clinical terms are non-facial annotations of the true syndrome, each sampled with its annotated frequency, as in the FB2 dissertation. That design is optimistic.

Largest top-3 gains with 3 clinical terms:

| Syndrome | Shape only | 3 clinical terms |
|---|---|---|
| Fetal alcohol syndrome | 7% | 87% |
| Hypophosphatasia | 26% | 87% |
| Bardet–Biedl | 17% | 75% |
| Beckwith–Wiedemann | 19% | 75% |
| Spondyloepiphyseal dysplasia | 22% | 78% |
| Rett | 37% | 88% |
| Stickler | 41% | 91% |
| Loeys–Dietz | 41% | 89% |

## Reading

- **The FB2 strategy holds on GNM.** HPO priors from findings that carry information independent of face shape raise top-3 from 57% to 81%, with no loss of specificity for unaffected faces.
- **HPO terms measured from the same face add nothing to a shape classifier.** They are redundant with the shape it already sees.
  - Under the FB2 rule they hurt, because Non-syndromic is always in the favoured set and measured present terms are often not annotated for the true syndrome. That pushes syndromic faces toward Non-syndromic.
  - The calibrated likelihood is neutral.
- **The value of facial HPO terms is elsewhere:** as the bridge from the face to about 7,000 annotated diseases (the differential), for exclusions, and for explaining the call. It is not as a booster for a classifier that already sees the whole face.

## Audit of the FB2 simulation protocol (`protocol_audit.R`, 20 replicates)

One HPO term per face, on the same out-of-fold posteriors. Rates are per face (micro) and by-syndrome mean (macro, as in chapter 4).

| Protocol | Macro top-1 | Macro top-3 | Micro top-1 | Unaffected correct |
|---|---|---|---|---|
| Shape only | 31.0% | 47.7% | 34.9% | 96.8% |
| Prevalence 1, exact-match favoured sets (loocv_HPO_sim.R) | 40.4% | 58.4% | 46.4% | 96.8% |
| Prevalence 1, ontology-aware favoured sets | 40.0% | 57.9% | 46.0% | 96.8% |
| Term applied with probability 1 − frequency (prevalence_simulation_job.R as written) | 36.4% | 53.7% | 41.5% | 96.8% |
| Term applied with probability = frequency (corrected) | 35.1% | 52.2% | 39.8% | 96.8% |
| Corrected, plus a wrong term for 10% of faces | 34.6% | 51.8% | 39.0% | 96.6% |
| Corrected, plus a wrong term for 25% of faces | 34.0% | 51.5% | 38.1% | 96.2% |
| Corrected, plus a wrong term for 50% of faces | 33.1% | 51.1% | 36.5% | 95.6% |

Unaffected faces get wrong terms at the same rate. A wrong term is one annotated to another class and not to the face's own.

## Files

`fit_hdrda.R` also scores two prior arms internally. Its FB2 arm does not keep the model prior when a face has no informative term, so use `evaluate_priors.R` and `protocol_audit.R` for the prior results. The scripts use absolute paths to a local working directory.

The code is `fit_hdrda.R`, `evaluate_priors.R` and `export_measured.py`. The other outputs contain participant-level data and are not shareable:
- `cv_posteriors_80pc.rds`, `cv_predictions_80pc.csv`: participant-level.
- `cv_per_class_80pc.csv`, `priors_per_class_top3_80pc.csv`, `cv_summary_80pc.json`, `priors_summary_80pc.json`: aggregate.
