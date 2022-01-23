
library(sparsediscrim)
library(caret)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
levels(hdrda.df$synd) <- levels(d.meta$Syndrome)
in.hpo <- rep(NA, length(unique(hdrda.df$synd)))

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

#NAs are defined as obligate####
na.prevalence <- 1
standardized.freqs$Frequency[is.na(standardized.freqs$Frequency)] <- na.prevalence

#testing the method with one HPO term at a time with simulated term prevalence####
#for all the people with syndrome i, let's look at the sensitivity with HPO term j
#how to deal with differing number of terms for each synd? Save only the mean top1,3,10 sens and the min/max, and the number of HPO terms associated

hdrda.orig <- hdrda(synd ~ ., data = hdrda.df)
synd.hpo.result <- NULL
for(i in 1 : length(unique(hdrda.df$synd))){
  
  N.hpo <- length(unique(hpo.pos[hpo.pos$V3 == official.names[i],5]))
  if(N.hpo > 0){
    
    for(j in 1 : N.hpo){
      
      #make vector to store perf for each term
      hpo.term <- as.character(unique(hpo.pos[hpo.pos$V3 == official.names[i],5])[j])
      #see hpo name: hpo$name[names(hpo$name) == hpo.term]
      for(k in 1:length(in.hpo)) in.hpo[k] <- length(grep(official.names[k], x = hpo.pos[hpo.pos[,5] == hpo.term,3])) > 0
      
      in.hpo[official.names == "Non-syndromic"] <- TRUE
      
      #get frequency of term with current syndrome####
      
      tmp.hpo.frequency <- as.numeric(standardized.freqs$Frequency[standardized.freqs$term == hpo.term][grep(official.names[i], standardized.freqs[standardized.freqs$term == hpo.term,1], ignore.case = T)])#what's the frequency of the selected term? phenotype.df$Frequency[phenotype.df$HPO_ID == hpo.term][grep(official.names[i], phenotype.df[phenotype.df$HPO_ID == hpo.term,2], ignore.case = T)]
      if(length(tmp.hpo.frequency) == 0) tmp.hpo.frequency <- na.prevalence
      if(tmp.hpo.frequency > 0){
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
        
        hdrda.updated <- hdrda(synd ~ ., data = hdrda.df, prior = updated.priors)
        
        #calculate an index of which observations get the HPO term boost
        synd.count <- 1:nrow(hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],])
        updated.posterior.sample <- sample(synd.count, length(synd.count)*tmp.hpo.frequency)
        
        posterior.distribution <- predict(hdrda.updated, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$post
        orig.posterior.distribution <- predict(hdrda.mod, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$post
        
        #debug: View(data.frame(posterior.distribution[1,], orig.posterior.distribution[1,]))
        
        posterior.class <- predict(hdrda.updated, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$class
        levels.to.keep <- levels(posterior.class)
        old.posterior.class <- predict(hdrda.orig, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$class
        
        posterior.class <- as.character(posterior.class)
        old.posterior.class <- as.character(old.posterior.class)
        
        #merge simulated hpo data w/original data
        posterior.distribution[-updated.posterior.sample,] <- orig.posterior.distribution[-updated.posterior.sample,]
        posterior.class[-updated.posterior.sample] <- old.posterior.class[-updated.posterior.sample]
        
        posterior.class <- factor(posterior.class, levels = levels.to.keep)
        tmp.gs <- factor((hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i], 1]), levels = levels(posterior.class))
      } else {
        
        posterior.distribution <- predict(hdrda.mod, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")
        posterior.class <- predict(hdrda.mod, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1])$class
        tmp.gs <- factor((hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i], 1]), levels = levels(posterior.class))
      }
      
      #what information do I want to keep about this synd * term combination? syndrome, term name, identifier, Sensitivity, first 5 rows of the mistaken table?
      if(length(hpo$name[names(hpo$name) == hpo.term]) == 0) tmp.hpo.name <- "NO NAME" else{tmp.hpo.name <- hpo$name[names(hpo$name) == hpo.term]}
      tmp.result <- data.frame(synd = levels(hdrda.df$synd)[i], hpo.id = hpo.term, hpo.name = tmp.hpo.name, sensitivity = confusionMatrix(posterior.class, tmp.gs)$byClass[which(levels(tmp.gs) == tmp.gs[1]),1])
      synd.hpo.result <- rbind.data.frame(synd.hpo.result, tmp.result)#, predicted.classes = t(sort(confusionMatrix(posterior.class, tmp.gs)$table[,which(official.names == official.names[i])], decreasing = T)[1:5])))
      
      print(tmp.result)
      
    }
    
  }
  save(synd.hpo.result, file = "adjusted_hpo_results_NA_1.Rdata")
}

save(synd.hpo.result, file = "adjusted_hpo_results_NA_1.Rdata")








