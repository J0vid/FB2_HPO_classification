# FB2_HPO_classification
Syndrome classification with dynamic priors

# Here are the important files for constructing figures from the simulations
- permuted_results_ranked_plots.Rdata: this file has the top 1, 3, and 10 results for 1000 permutations of syndrome X hpo term uses. All simulations used a term prevalence of .5 if there was no known prevalence for the term for a given syndrome. This was run from Dense_HPO_analysis.R

- hpo_results_loocv_full.Rdata: this file has the main syndrome X hpo term calculation. I ended up using a prevalence of 1 (assigning each term X syndrome) here and simulating prevalence before plotting. It's much more computationally efficient that way. This was run as a job using loocv_prevalence_sim.R

- loocv_hdrda.Rdata: this file contains the predictions and posteriors for each iteration of the face shape only HDRDA model. This was run as a job using loocv_training_job.R

- adjusted_data_combined.Rdata: this file has everything you need to get going. Here are it's contents from the original save call -- # save(age.sex.lm, atlas, d.meta, front.face, PC.eigenvectors, PC.scores, synd.mshape, phenotype.df.synd, hpo, hpo.pos, hdrda.mod, hdrda.df, official.names, file = "adjusted_data_combined.Rdata"). The PC scores used to train models were adjusted with a linear model PCs[,1:200] ~ sex + poly(age, 3). Initial processing to get that starting data file was done in data_processing.R

- data_processing.R: this is where I combined the non-syndromic and syndromic dense landmark data, combined the metadata, did a PCA and adjusted the PC scores with the age and sex model described above.

