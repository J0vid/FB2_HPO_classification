#loocv prevalence simulation job script####

library(sparsediscrim)
library(caret)

load("hpo_results_loocv_full.Rdata")

#only 26 syndromes have frequency data associated with HPO terms: unique(phenotype.df.synd[is.na(phenotype.df.synd$Frequency) == F, "DiseaseName"])
#that tells me that we need to run a simulation on what it looks like if terms have any of the frequency hpo ranges
#let's define this behavior
standardized.freqs <- data.frame(synd = phenotype.df.synd$DiseaseName, term = phenotype.df.synd$HPO_ID, Frequency = phenotype.df.synd$Frequency)
standardized.freqs$Frequency <- as.character(standardized.freqs$Frequency)

standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040280"] <- 1
standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040281"] <- (80+99)/200
standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040282"] <- (30+79)/200
standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040283"] <- (5+29)/200
standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040284"] <- (1+4)/200
standardized.freqs$Frequency[standardized.freqs$Frequency == "HP:0040285"] <- 0
#deal with fractions
for(i in grep("/", standardized.freqs$Frequency)) standardized.freqs$Frequency[i] <- eval(parse(text = standardized.freqs$Frequency[grep("/", standardized.freqs$Frequency)[i]]))
#deal with percentages
for(i in grep("%", standardized.freqs$Frequency)) standardized.freqs$Frequency[i] <- as.numeric(substr(standardized.freqs$Frequency[i], 1, nchar(standardized.freqs$Frequency[i])-1))/100

# View(standardized.freqs)

#NAs are defined as obligate####
standardized.freqs$Frequency[is.na(standardized.freqs$Frequency)] <- 1

#testing the method with one HPO term at a time with simulated term prevalence####
#for all the people with syndrome i, let's look at the sensitivity with HPO term j
#how to deal with differing number of terms for each synd? Save only the mean top1,3,10 sens and the min/max, and the number of HPO terms associated

hdrda.orig <- hdrda(synd ~ ., data = hdrda.df)
# synd.hpo.result <- NULL
# hpo.pred <- NULL
# hpo.distribution <- NULL
# hpo.meta <- NULL
for(i in 19 : length(unique(hdrda.df$synd))){
  
  N.hpo <- length(unique(hpo.pos[hpo.pos$V3 == official.names[i],5]))
  # N.hpo <- length(unique(standardized.freqs$term[standardized.freqs$synd == official.names[i]]))
  if(N.hpo > 0){
    if(N.hpo >= 5){hpo.sampler <- sample(1:N.hpo, 5);N.hpo <- 5} else{hpo.sampler <- 1:N.hpo}
    for(j in 1 : N.hpo){
      
      #make vector to store perf for each term
      hpo.term <- as.character(unique(hpo.pos[hpo.pos$V3 == official.names[i],5])[hpo.sampler[j]])
      #see all hpo names for current synd: 
      # View(data.frame(hpo$name[names(hpo$name) %in% unique(standardized.freqs$term[standardized.freqs$synd == official.names[i]])]))
      # View(data.frame(hpo$name[names(hpo$name) %in% unique(hpo.pos[hpo.pos$V3 == official.names[i],5])]))
      #see current hpo name: hpo$name[names(hpo$name) == hpo.term]
      
      for(k in 1:length(in.hpo)) in.hpo[k] <- length(grep(official.names[k], x = hpo.pos[hpo.pos[,5] == hpo.term,3])) > 0
      
      in.hpo[official.names == "Non-syndromic"] <- TRUE
      #debug: levels(hdrda.df$synd)[in.hpo]
      
      #get frequency of term with current syndrome####
      tmp.hpo.frequency <- as.numeric(standardized.freqs$Frequency[standardized.freqs$term == hpo.term][grep(official.names[i], standardized.freqs[standardized.freqs$term == hpo.term,1], ignore.case = T)])#what's the frequency of the selected term? phenotype.df$Frequency[phenotype.df$HPO_ID == hpo.term][grep(official.names[i], phenotype.df[phenotype.df$HPO_ID == hpo.term,2], ignore.case = T)]
      if(length(tmp.hpo.frequency) == 0) tmp.hpo.frequency <- 1
      
      #priors are adjusted to uniformly sharing 20% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list
      updated.priors <- rep(NA, length(official.names))
      names(updated.priors) <- official.names
      
      updated.priors[in.hpo == F] <- 0.1 / length(official.names[in.hpo == F])
      updated.priors[in.hpo == T] <- 0.9 / length(official.names[in.hpo == T])
      
      #round updated.priors to avoid floating point error in summing to 1
      priorsequal1 <- 0
      z <- 9
      while(priorsequal1 != 1){
        z <- z +1
        updated.priors <- round(updated.priors, digits = z)
        updated.priors <- updated.priors/sum(updated.priors)
        priorsequal1 <- sum(updated.priors)
      }
      
      #calculate an index of which observations get the HPO term boost
      synd.count <- 1:nrow(hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],])
      
      current.synd.index <- which(hdrda.df$synd == levels(hdrda.df$synd)[i])
      updated.posterior.sample <- sample(current.synd.index, length(synd.count)*tmp.hpo.frequency)
      
      #create matrices to store updated loocv preds
      posterior.distribution <- loocv.post[hdrda.df$synd == levels(hdrda.df$synd)[i],]
      posterior.class <- loocv.pred[hdrda.df$synd == levels(hdrda.df$synd)[i]]
      
      for(l in 1:length(updated.posterior.sample)){ #only loop through selected inds
        hdrda.updated <- hdrda(synd ~ ., data = hdrda.df[-updated.posterior.sample[l],], prior = updated.priors)
        tmp.updated.pred <- predict(hdrda.updated, newdata = hdrda.df[updated.posterior.sample[l],-1], type = "prob")
        posterior.class[which(current.synd.index == updated.posterior.sample[l])] <- as.character(tmp.updated.pred$class)
        posterior.distribution[which(current.synd.index == updated.posterior.sample[l]),] <- tmp.updated.pred$posterior
      }
      
      # View(posterior.distribution)
      
      posterior.class <- factor(posterior.class, levels = levels(hdrda.df$synd))
      tmp.gs <- factor((hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i], 1]), levels = levels(posterior.class))
      
      #what information do I want to keep about this synd * term combination? syndrome, term name, identifier, Sensitivity, first 5 rows of the mistaken table?
      if(length(hpo$name[names(hpo$name) == hpo.term]) == 0) tmp.hpo.name <- "NO NAME" else{tmp.hpo.name <- hpo$name[names(hpo$name) == hpo.term]}
      
      tmp.result <- data.frame(synd = levels(hdrda.df$synd)[i], hpo.id = hpo.term, hpo.name = tmp.hpo.name, sensitivity = confusionMatrix(posterior.class, tmp.gs)$byClass[which(levels(tmp.gs) == tmp.gs[1]),1])
      synd.hpo.result <- rbind.data.frame(synd.hpo.result, tmp.result)#, predicted.classes = t(sort(confusionMatrix(posterior.class, tmp.gs)$table[,which(official.names == official.names[i])], decreasing = T)[1:5])))
      hpo.pred <- c(hpo.pred, as.character(posterior.class))
      hpo.distribution <- rbind(hpo.distribution, posterior.distribution)
      hpo.meta <- rbind(hpo.meta, cbind(as.character(rep(levels(hdrda.df$synd)[i], nrow(posterior.distribution))), rep(hpo.term, nrow(posterior.distribution)), rep(tmp.hpo.name, nrow(posterior.distribution))))
      
      print(tmp.result)
    }
    
  }
  save(synd.hpo.result, hpo.pred, hpo.distribution, hpo.meta, file = "hpo_results_loocv_full.Rdata")
}

save(synd.hpo.result, hpo.pred, hpo.distribution, hpo.meta, file = "hpo_results_loocv_full.Rdata")
# 