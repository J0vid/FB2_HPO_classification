library(sparsediscrim)
library(caret)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#begin prevalence sim####
#top choice####
load("hpo_results_loocv_full_80PC.Rdata")
colnames(hpo.distribution) <- levels(hdrda.df$synd)
ultimate.bunduru <- data.frame(hpo.meta, pred = hpo.pred, hpo.distribution)
colnames(ultimate.bunduru)[1:3] <- c("synd", "hpo.id", "hpo.term")
bunduru.final.key <- official.names[which(levels(hdrda.df$synd) %in% unique(ultimate.bunduru$synd))] #name correspondences for frequency info set and loocv training

permuted.results <- matrix(NA, nrow = nrow(ultimate.bunduru), ncol = 1000)

for(k in 1:1000){
  sim.synd.hpo.results <- NULL
  for(i in 1:length(unique(ultimate.bunduru$synd))){
    tmp.terms <- unique(ultimate.bunduru$hpo.term[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    tmp.terms.id <- unique(ultimate.bunduru$hpo.id[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    for(j in 1:length(unique(tmp.terms))){
      #grab synd i, term j and simulate hpo term prevalence based on the known value
      tmp.hpo.pred <- ultimate.bunduru$pred[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i] & ultimate.bunduru$hpo.term == tmp.terms[j]]
      
      hpo.prevalence <- final.freqs$Frequency[final.freqs$synd == bunduru.final.key[i] & final.freqs$term == tmp.terms.id[j]] ####need to solve this correspondence
      # draw the frequency from a uniform distribution if the prevalence is a range
      if(hpo.prevalence == "HP:0040281") hpo.prevalence <- runif(1, 80, 99)/100
      if(hpo.prevalence == "HP:0040282") hpo.prevalence <- runif(1, 30, 79)/100
      if(hpo.prevalence == "HP:0040283") hpo.prevalence <- runif(1, 5, 29)/100
      if(hpo.prevalence == "HP:0040284") hpo.prevalence <- runif(1, 1, 4)/100
      hpo.prevalence <- as.numeric(hpo.prevalence)
      
      tmp.prevalence <- sample(1:length(tmp.hpo.pred), length(tmp.hpo.pred) * hpo.prevalence)
      tmp.hpo.pred[tmp.prevalence] <- loocv.pred[hdrda.df$synd == unique(ultimate.bunduru$synd)[i]][tmp.prevalence]
      sim.synd.hpo.results <- c(sim.synd.hpo.results, tmp.hpo.pred)
    }
  }
  permuted.results[,k] <- sim.synd.hpo.results
  save(permuted.results, file = "80PC_permuted_results_rank1.Rdata")
  
  print(k)
}
save(permuted.results, file = "80PC_permuted_results_rank1.Rdata")


#top 3####
# define equivalent loocv result
#rank2 for original model####
rank2 <- 3
rank2.result <- loocv.pred
for(i in 1:length(unique(hdrda.df$synd))){
  tmp.hpo.post <- loocv.post[hdrda.df$synd == levels(hdrda.df$synd)[i],]
  tmp.hpo.pred <-loocv.pred[hdrda.df$synd == levels(hdrda.df$synd)[i]]
  rank2.check <- rep(NA, nrow(tmp.hpo.post))
  for(l in 1: length(rank2.check)){
    if(sum(levels(hdrda.df$synd)[order(tmp.hpo.post[l,], decreasing = T)][1:rank2] == levels(hdrda.df$synd)[i]) > 0){
      rank2.check[l] <- as.character(levels(hdrda.df$synd)[i])
    } else{rank2.check[l] <- tmp.hpo.pred[l]}
  }
  rank2.result[hdrda.df$synd == levels(hdrda.df$synd)[i]] <- rank2.check
}
face.only.rank2 <- rank2.result
permuted.results.rank2 <- matrix(NA, nrow = nrow(ultimate.bunduru), ncol = 1000)

#top 3 simulation####
for(k in 1:1000){
  sim.synd.hpo.results <- NULL
  for(i in 1:length(unique(ultimate.bunduru$synd))){
    tmp.terms <- unique(ultimate.bunduru$hpo.term[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    tmp.terms.id <- unique(ultimate.bunduru$hpo.id[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    
    for(j in 1:length(unique(tmp.terms))){
      #grab synd i, term j posteriors and check top 3 posteriors. Then simulate .5 hpo term prevalence
      tmp.hpo.post <- ultimate.bunduru[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i] & ultimate.bunduru$hpo.term == tmp.terms[j], -1:-4]
      rank2 <- 3
      rank2.check <- rep(NA, nrow(tmp.hpo.post))
      for(l in 1: length(rank2.check)){
        if(sum(levels(hdrda.df$synd)[order(tmp.hpo.post[l,], decreasing = T)][1:rank2] == unique(ultimate.bunduru$synd)[i]) > 0){
          rank2.check[l] <- as.character(unique(ultimate.bunduru$synd)[i])
        } else{rank2.check[l] <- names(sort(tmp.hpo.post[l,], decreasing = T)[1])}
      }
      
      hpo.prevalence <- final.freqs$Frequency[final.freqs$synd == bunduru.final.key[i] & final.freqs$term == tmp.terms.id[j]] ####need to solve this correspondence
      # draw the frequency from a uniform distribution if the prevalence is a range
      if(hpo.prevalence == "HP:0040281") hpo.prevalence <- runif(1, 80, 99)/100
      if(hpo.prevalence == "HP:0040282") hpo.prevalence <- runif(1, 30, 79)/100
      if(hpo.prevalence == "HP:0040283") hpo.prevalence <- runif(1, 5, 29)/100
      if(hpo.prevalence == "HP:0040284") hpo.prevalence <- runif(1, 1, 4)/100
      hpo.prevalence <- as.numeric(hpo.prevalence)
      
      tmp.prevalence <- sample(1:length(rank2.check), length(rank2.check) * hpo.prevalence)
      rank2.check[tmp.prevalence] <- face.only.rank2[hdrda.df$synd == unique(ultimate.bunduru$synd)[i]][tmp.prevalence]
      sim.synd.hpo.results <- c(sim.synd.hpo.results, rank2.check)
    }
  }
  permuted.results.rank2[,k] <- sim.synd.hpo.results
  save(permuted.results.rank2, file = "80PC_permuted_results_rank2.Rdata")
  print(k)
}

save(permuted.results.rank2, file = "80PC_permuted_results_rank2.Rdata")

#top 10####
# define equivalent loocv result
#rank3 for original model####
rank3 <- 10
rank3.result <- loocv.pred

colnames(loocv.post) <- levels(hdrda.df$synd)

for(i in 1:length(unique(hdrda.df$synd))){
  tmp.hpo.post <- loocv.post[hdrda.df$synd == levels(hdrda.df$synd)[i],]
  tmp.hpo.pred <-loocv.pred[hdrda.df$synd == levels(hdrda.df$synd)[i]]
  rank3.check <- rep(NA, nrow(tmp.hpo.post))
  for(l in 1: length(rank3.check)){
    if(sum(levels(hdrda.df$synd)[order(tmp.hpo.post[l,], decreasing = T)][1:rank3] == levels(hdrda.df$synd)[i]) > 0){
      rank3.check[l] <- as.character(levels(hdrda.df$synd)[i])
    } else{rank3.check[l] <- tmp.hpo.pred[l]}
    # print(l)
    # print(tmp.hpo.post[l,order(tmp.hpo.post[l,], decreasing = T)][1:rank3])
  }
  rank3.result[hdrda.df$synd == levels(hdrda.df$synd)[i]] <- rank3.check
}
face.only.rank3 <- rank3.result


permuted.results.rank3 <- matrix(NA, nrow = nrow(ultimate.bunduru), ncol = 1000)

for(k in 1:1000){
  sim.synd.hpo.results <- NULL
  for(i in 1:length(unique(ultimate.bunduru$synd))){
    tmp.terms <- unique(ultimate.bunduru$hpo.term[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    tmp.terms.id <- unique(ultimate.bunduru$hpo.id[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    
    for(j in 1:length(unique(tmp.terms))){
      #grab synd i, term j posteriors and check top 3 posteriors. Then simulate .5 hpo term prevalence
      tmp.hpo.post <- ultimate.bunduru[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i] & ultimate.bunduru$hpo.term == tmp.terms[j], -1:-4]
      rank3 <- 10
      rank3.check <- rep(NA, nrow(tmp.hpo.post))
      for(l in 1: length(rank3.check)){
        if(sum(levels(hdrda.df$synd)[order(tmp.hpo.post[l,], decreasing = T)][1:rank3] == unique(ultimate.bunduru$synd)[i]) > 0){
          rank3.check[l] <- as.character(unique(ultimate.bunduru$synd)[i])
        } else{rank3.check[l] <- names(sort(tmp.hpo.post[l,], decreasing = T)[1])}
        # print(l)
        # print(unique(ultimate.bunduru$synd)[i])
        # print(tmp.hpo.post[l,order(tmp.hpo.post[l,], decreasing = T)][1:rank3])
      }
      
      
      hpo.prevalence <- final.freqs$Frequency[final.freqs$synd == bunduru.final.key[i] & final.freqs$term == tmp.terms.id[j]] ####need to solve this correspondence
      # draw the frequency from a uniform distribution if the prevalence is a range
      if(hpo.prevalence == "HP:0040281") hpo.prevalence <- runif(1, 80, 99)/100
      if(hpo.prevalence == "HP:0040282") hpo.prevalence <- runif(1, 30, 79)/100
      if(hpo.prevalence == "HP:0040283") hpo.prevalence <- runif(1, 5, 29)/100
      if(hpo.prevalence == "HP:0040284") hpo.prevalence <- runif(1, 1, 4)/100
      hpo.prevalence <- as.numeric(hpo.prevalence)
      
      tmp.prevalence <- sample(1:length(rank3.check), length(rank3.check) * hpo.prevalence)
      rank3.check[tmp.prevalence] <- face.only.rank3[hdrda.df$synd == unique(ultimate.bunduru$synd)[i]][tmp.prevalence]
      sim.synd.hpo.results <- c(sim.synd.hpo.results, rank3.check)
    }
    
  }
  permuted.results.rank3[,k] <- sim.synd.hpo.results
  save(permuted.results.rank3, file = "80PC_permuted_results_rank3.Rdata")
  print(k)
}

save(permuted.results.rank3, file = "80PC_permuted_results_rank3.Rdata")














