library(sparsediscrim)
library(caret)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#let's define frequencies
phenotype.df.synd2 <- phenotype_2022[phenotype_2022$DiseaseName %in% official.names, ]

standardized.freqs <- data.frame(synd = phenotype.df.synd2$DiseaseName, term = phenotype.df.synd2$HPO_ID, Frequency = phenotype.df.synd2$Frequency)
standardized.freqs$Frequency <- as.character(standardized.freqs$Frequency)

standardized.freqs <- data.frame(standardized.freqs[is.na(standardized.freqs$Frequency) == F,])

# how bad is it to throw out NA freqs? 3 synds don't have any frequency info for any terms
# View(data.frame(table(standardized.freqs$synd[is.na(standardized.freqs$Frequency) == F])))

#for terms with ranges, let's uniformly draw from that range.
standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040280"] <- 1
for(i in which(standardized.freqs$Frequency == "HP:0040281")) standardized.freqs$Frequency[i] <- runif(1, 80, 99)/100 #(80+99)/200
for(i in which(standardized.freqs$Frequency == "HP:0040282")) standardized.freqs$Frequency[i] <- runif(1, 30, 79)/100 #(30+79)/200
for(i in which(standardized.freqs$Frequency == "HP:0040283")) standardized.freqs$Frequency[i] <- runif(1, 5, 29)/100 #(5+29)/200
for(i in which(standardized.freqs$Frequency == "HP:0040284")) standardized.freqs$Frequency[i] <- runif(1, 1, 4)/100 #(1+4)/200
standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040285"] <- 0
#deal with fractions: fractions <- grep("/", standardized.freqs$Frequency)
for(i in grep("/", standardized.freqs$Frequency)) standardized.freqs$Frequency[i] <- eval(parse(text = standardized.freqs$Frequency[i]))
#deal with percentages: percs <- grep("%", standardized.freqs$Frequency)
for(i in grep("%", standardized.freqs$Frequency)) standardized.freqs$Frequency[i] <- as.numeric(substr(standardized.freqs$Frequency[i], 1, nchar(standardized.freqs$Frequency[i])-1))/100

# View(standardized.freqs)

#Now a dataframe where we exclusively deal with NA frequencies####
na.freqs <- data.frame(synd = phenotype.df.synd2$DiseaseName, term = phenotype.df.synd2$HPO_ID, Frequency = phenotype.df.synd2$Frequency)
na.freqs$Frequency <- as.character(na.freqs$Frequency)
na.freqs <- data.frame(na.freqs[is.na(na.freqs$Frequency) == T,])

#NAs are defined as .3####
na.freqs$Frequency[is.na(na.freqs$Frequency)] <- .3

View(na.freqs)

final.freqs <- rbind.data.frame(standardized.freqs, na.freqs)

colnames(hpo.distribution) <- levels(hdrda.df$synd)
ultimate.bunduru <- data.frame(hpo.meta, pred = hpo.pred, hpo.distribution)
colnames(ultimate.bunduru)[1:3] <- c("synd", "hpo.id", "hpo.term")

permuted.results <- matrix(NA, nrow = nrow(ultimate.bunduru), ncol = 1000)

for(k in 1:1000){
  sim.synd.hpo.results <- NULL
  for(i in 1:length(unique(ultimate.bunduru$synd))){
    tmp.terms <- unique(ultimate.bunduru$hpo.term[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    for(j in 1:length(unique(tmp.terms))){
      #grab synd i, term j and simulate .5 hpo term prevalence
      tmp.hpo.pred <- ultimate.bunduru$pred[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i] & ultimate.bunduru$hpo.term == tmp.terms[j]]
      hpo.prevalence <- final.freqs$term[final.freqs$synd == unique(ultimate.bunduru$synd)[i]] ####need to solve this correspondence
      tmp.prevalence <- sample(1:length(tmp.hpo.pred), length(tmp.hpo.pred) * .5)
      tmp.hpo.pred[tmp.prevalence] <- loocv.pred[hdrda.df$synd == unique(ultimate.bunduru$synd)[i]][tmp.prevalence]
      sim.synd.hpo.results <- c(sim.synd.hpo.results, tmp.hpo.pred)
    }
  }
  permuted.results[,k] <- sim.synd.hpo.results
  print(k)
}



