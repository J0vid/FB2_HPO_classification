#HPO analysis on dense PCs####
library(ggplot2)
library(Morpho)
library(geomorph)
library(rgl)
library(dplyr)
library(Rvcg)
library(sparsediscrim)
library(caret)
library(readr)
library(plotly)

# load("Classification_demo_legacy/demo_objects.Rdata")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# save(atlas, d.meta, front.face, PC.eigenvectors, PC.scores, synd.mshape, phenotype.df.synd, hpo, hpo.pos, hdrda.mod, hdrda.df, official.names, file = "data_combined.Rdata")
# save(age.sex.lm, atlas, d.meta, front.face, PC.eigenvectors, PC.scores, synd.mshape, phenotype.df.synd, hpo, hpo.pos, hdrda.mod, hdrda.df, official.names, file = "adjusted_data_combined.Rdata")

# load("data.Rdata")
load("adjusted_data_combined.Rdata")
#load these in sequence to update database then save final combined.rdata file
#load("combined_PCs.Rdata")
#load("adjusted_PCs.Rdata")
#PC.scores <- adjusted.PC.scores

# orig.preds <- (predict(hdrda.mod, hdrda.df[,-1]))
# View(confusionMatrix(factor(loocv.pred, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass)

#start with a one-off example of how changing priors changes an individual's posterior distribution####
#when someone selects specific terms, find the omims/syndrome name in hpo.pos
print(hpo.pos[hpo.pos[,5] == names(hpo$name)[hpo$name == "Strabismus"],3])

in.hpo <- rep(NA, length(unique(hdrda.df$synd)))

for(i in 1:length(in.hpo)) in.hpo[i] <- length(grep(official.names[i], x = hpo.pos[hpo.pos[,5] == "HP:0009601",3])) > 0

in.hpo[official.names == "Non-syndromic"] <- TRUE

#priors are adjusted to uniformly sharing x% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list

updated.priors <- rep(NA, length(official.names))
names(updated.priors) <- official.names

updated.priors[in.hpo == F] <- .05 / length(official.names[in.hpo == F])
updated.priors[in.hpo == T] <- .95 / length(official.names[in.hpo == T])

hdrda.updated <- hdrda(synd ~ ., data = hdrda.df, prior = updated.priors)

posterior.distribution <- predict(hdrda.updated, newdata = hdrda.df[which(d.meta$Syndrome == "Nager Syndrome"),-1], type = "prob")$post[3,]

View(predict(hdrda.updated, newdata = hdrda.df[which(d.meta$Syndrome == "Nager Syndrome"),-1], type = "prob")$post)
predict(hdrda.updated, newdata = hdrda.df[which(d.meta$Syndrome == "Nager Syndrome"),-1], type = "prob")$class

posterior.distribution <- sort(posterior.distribution, decreasing = T)

#used to be part of plot.df: ID = as.factor(1:10), 
plot.df <- data.frame(Probs = round(as.numeric(posterior.distribution[1:10]), digits = 4), Syndrome = as.factor(gsub("_", " ", names(posterior.distribution[1:10]))))
plot.df$Syndrome <- as.character(plot.df$Syndrome)
plot.df$Syndrome[plot.df$Syndrome == "Unaffected Unrelated"] <- "Non-syndromic"


plot_ly(data = plot.df, x = ~Syndrome, y = ~Probs, type = "bar", color = I("grey"), hoverinfo = paste0("Syndrome: ", "x", "<br>", "Probability: ", "y")) %>%
  layout(xaxis = list(tickvals = gsub("_", " ", plot.df$Syndrome), tickangle = 45, ticktext = c(Syndrome = plot.df$Syndrome, Probability = plot.df$Probs), title = "<b>Syndrome</b>"),
         yaxis = list(title = "<b>Class probability</b>"),
         paper_bgcolor='rgba(245, 245, 245, .9)',
         margin = list(b = 125, l = 50, r = 100)
  )

#can we incorporate frequencies?####
for(i in 1 : length(unique(hdrda.df$synd))){
  for(j in 1 : length(unique(hpo.pos[hpo.pos$V3 == official.names[i],5]))){
  hpo.term <- unique(hpo.pos[hpo.pos$V3 == official.names[i],5])[j]

#get frequency of term with current syndrome####
tmp.hpo.frequency <- phenotype.df$Frequency[phenotype.df$HPO_ID == hpo.term][grep(official.names[i], phenotype.df[phenotype.df$HPO_ID == hpo.term,2], ignore.case = T)]
print(paste0(levels(hdrda.df$synd)[i], ": [", hpo$name[names(hpo$name) == hpo.term], "] Freq:", tmp.hpo.frequency,"."))
}
}

#most frequencies are blank--what should we do? We could limit the analysis to only the syndromes and terms with frequency stats or we could simulate frequency...####
# phenotype.df.synd <- phenotype.df[phenotype.df$DiseaseName %in% official.names,]
View(phenotype.df.synd[is.na(phenotype.df.synd$Frequency) == F,])




#plot result priors with simulated term prevalences####
#outside of one-off examples, I have moved the all of the heavy computation out of this script and into job scripts. Visualization continues below.
load("hpo_results_NA_54_4.Rdata")
load("full_hpo_results_NA_545.Rdata") #pre-loocv results

hpo.perf <- synd.hpo.result %>%
  group_by(synd) %>%
  summarise(top1.min = min(sensitivity), top1.max = max(sensitivity), top1.mean = mean(sensitivity))

hdrda.orig.preds <- predict(hdrda.mod, newdata = hdrda.df[,-1])$class
original.sens <- confusionMatrix(hdrda.orig.preds, hdrda.df$synd)$byClass[,1]

hpo.df <- data.frame(orig.sens = original.sens[match(hpo.perf$synd, levels(hdrda.df$synd))], hpo.perf)

ggplot(aes(x = reorder(synd, -orig.sens), y = top1.mean), data = hpo.df) +
  geom_bar(stat = "identity",  fill = "slategrey") + 
  geom_errorbar(aes(ymin = top1.min, ymax = top1.max), width=.2, position=position_dodge(.9)) +
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")

#we should mark non-informative hpo terms
load("C:/Users/David A/Downloads/FB2_HPO_classification/hpo_results_NA_54_4.Rdata")
View(data.frame(unique(synd.hpo.result$hpo.name)))
bad.hpos <- c("Autosomal dominant inheritance", "X-linked dominant inheritance", "Autosomal recessive inheritance", "X-linked inheritance", "X-linked recessive inheritance", "Variable expressivity", "Stillbirth")
synd.hpo.result_5 <- synd.hpo.result
#plot result priors with varying simulated term prevalences####
hpo.perf_1 <- synd.hpo.result_1 %>%
  group_by(synd) %>%
  summarise(top1.min = min(sensitivity), top1.max = max(sensitivity), top1.mean = mean(sensitivity))

hpo.perf_5 <- synd.hpo.result_5 %>%
  group_by(synd) %>%
  summarise(top1.min = min(sensitivity), top1.max = max(sensitivity), top1.mean = mean(sensitivity))

hpo.perf_25 <- synd.hpo.result_25 %>%
  group_by(synd) %>%
  summarise(top1.min = min(sensitivity), top1.max = max(sensitivity), top1.mean = mean(sensitivity))


hdrda.orig.preds <- predict(hdrda.mod, newdata = hdrda.df[,-1])$class
levels(hdrda.orig.preds) <- levels(hdrda.df$synd)
original.sens <- confusionMatrix(factor(loocv.pred, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1]

hpo.df <- data.frame(orig.sens = original.sens[match(hpo.perf_1$synd, levels(hdrda.df$synd))], synd = hpo.perf_1$synd, mean1 = hpo.perf_1$top1.mean, mean5 = hpo.perf_5$top1.mean, mean25 = hpo.perf_25$top1.mean)
hpo.df <- data.frame(orig.sens = original.sens[match(hpo.perf_5$synd, levels(hdrda.df$synd))], synd = hpo.perf_5$synd, mean1 = hpo.perf_5$top1.mean)
View(hpo.df)
ggplot(aes(x = reorder(synd, -orig.sens), y = mean1), data = hpo.df) +
  geom_bar(stat = "identity",  fill = "#0F084B") + 
  # geom_bar(stat = "identity", aes(y = mean5), fill = "#6066AC", alpha = 1) +
  # geom_bar(stat = "identity", aes(y = mean25), fill = "#B0B5EF", alpha = 1) +
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#DFE1F9") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  # scale_colour_manual(name = "Classification /n approach", values=c("#0F084B", "slategrey"), labels = c("Shape only", "With HPO")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")

#what insights can we glean about the utility of hpo terms?####
#relationship between sensitivity and number of syndromes associated?
#rank the differences in hpo sensitivity (top 10 most useful)

#visualizations for how posteriors change with term info####
load("C:/Users/David A/Downloads/FB2_HPO_classification/hpo_results_loocv_full.Rdata")

#demonstrative syndromes: Marfan Syndrome, Kabuki Syndrome
#hists facet wrapped by term don't show how inds change
selected.synd <- "Stickler Syndrome"
  
hist.df <- data.frame(term = hpo.meta[hpo.meta[,1] == selected.synd,3], X1 =  hpo.distribution[hpo.meta[,1] == selected.synd, which(levels(hdrda.df$synd) == selected.synd)])
hist.df <- rbind.data.frame(data.frame(term = "No term", X1 = loocv.post[hdrda.df$synd == selected.synd, which(levels(hdrda.df$synd) == selected.synd)]),
hist.df)

hist.df$term <- factor(hist.df$term, levels = c("No term", unique(hist.df$term)[unique(hist.df$term) != "No term"]))

ggplot(aes(x = X1, y = ..density..), data = hist.df) +
  geom_histogram(bins = 10) +
  xlab(paste0("Posterior values for ", selected.synd)) +
  ylab("Frequency") +
  theme_bw() +
  facet_wrap(term ~ ., strip.position = "top")

#per-term faceted syndromic hist comparison to non-syndromic hist####
#the trouble with this visualization is that it doesn't convey the increase in confidence for the true syndrome
selected.synd <- "Stickler Syndrome"

hist.df <- data.frame(term = hpo.meta[hpo.meta[,1] == selected.synd,3], X1 =  hpo.distribution[hpo.meta[,1] == selected.synd, c(1,which(levels(hdrda.df$synd) == selected.synd))])
hist.df <- rbind.data.frame(data.frame(term = "No term", X1 = loocv.post[hdrda.df$synd == selected.synd, c(1,which(levels(hdrda.df$synd) == selected.synd))]),
                            hist.df)

hist.df$term <- factor(hist.df$term, levels = c("No term", unique(hist.df$term)[unique(hist.df$term) != "No term"]))

colnames(hist.df)[2:3] <- c("Non-syndromic", selected.synd)

melt.df <- reshape2::melt(hist.df, id = "term")

ggplot(aes(x = value, y = ..density..), data = melt.df) +
  geom_histogram(bins = 10) +
  xlab(paste0("Posterior values for ", selected.synd)) +
  ylab("Frequency") +
  theme_bw() +
  facet_wrap(~ term + variable, strip.position = "top")

#violin plot showing per-term change in posteriors####
selected.synd <- "Stickler Syndrome"

hist.df <- data.frame(term = hpo.meta[hpo.meta[,1] == selected.synd,3], X1 =  hpo.distribution[hpo.meta[,1] == selected.synd, which(levels(hdrda.df$synd) == selected.synd)])
hist.df <- rbind.data.frame(data.frame(term = "No term", X1 = loocv.post[hdrda.df$synd == selected.synd, which(levels(hdrda.df$synd) == selected.synd)]),
                            hist.df)

hist.df$term <- factor(hist.df$term, levels = c("No term", unique(hist.df$term)[unique(hist.df$term) != "No term"]))

ggplot(aes(x = term, y = X1), data = hist.df) +
  geom_violin() +
  ylab(paste0("Posterior values for ", selected.synd)) +
  xlab("HPO term") +
  theme_bw()

# delta plot for each term####
#shows change for each term for each individual AND if they were correctly classified
#loop across all syndromes and save plot
for(j in 1:length(unique(hpo.meta[,1]))){
  #subtract no term from every term
  selected.synd <- unique(hpo.meta[,1])[j]

  hist.df <- data.frame(term = hpo.meta[hpo.meta[,1] == selected.synd,3], X1 =  hpo.distribution[hpo.meta[,1] == selected.synd,])
  no.term.df <- data.frame(X1 = loocv.post[hdrda.df$synd == selected.synd, ])
  
  delta.df <- NULL
  for(i in 1:length(unique(hist.df$term))) delta.df <- rbind(delta.df, hist.df[hist.df$term == unique(hist.df$term)[i], -1] - no.term.df[,])
  
  pred.colors <- c(paste0("Not ", selected.synd), selected.synd)[1 + (hpo.pred[hpo.meta[,1] == selected.synd] == selected.synd)]
  delta.df <- data.frame(term = hist.df$term, posterior = delta.df[, which(levels(hdrda.df$synd) == selected.synd)], correct.pred = pred.colors)
  
  pdf(paste0("results/ind_diffs/individual_term_changes_", gsub(pattern = "/", replacement = "_", selected.synd), ".pdf"), width = 10, height = 7)
  p <- ggplot(aes(x = term, y = posterior), data = delta.df) +
    geom_jitter(aes(colour = correct.pred), alpha = 1, cex = 2, width = 0.25) +
    scale_colour_manual(name = "Prediction", values = c("red", "black")) +
    ylab(paste0("Change in posterior values for ", selected.synd)) +
    xlab("HPO term") +
    ylim(-.5,1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  print(p)
  dev.off()
  
  print(j)
}

#for a given syndrome, how much does the term change top 1,5,10 preds####
#compare top 1 sensitivity####
#error bars####
hpo.perf <- synd.hpo.result %>%
  group_by(synd) %>%
  summarise(top1.min = min(sensitivity), top1.max = max(sensitivity), top1.mean = mean(sensitivity))


original.sens <- confusionMatrix(factor(loocv.pred, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1]
hpo.perf.full <- hpo.perf[match(levels(hdrda.df$synd), as.data.frame(hpo.perf[,1])$synd), ]
hpo.perf.full[,1] <- levels(hdrda.df$synd)


hpo.df <- data.frame(orig.sens = original.sens, hpo.perf.full)

fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/top1comparison_full.pdf", width = 12, height = 8)
ggplot(aes(x = reorder(synd, -orig.sens), y = top1.mean), data = hpo.df) +
  geom_bar(stat = "identity",  fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = top1.min, ymax = top1.max)) +
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  scale_colour_manual(name = "Classification /n approach", values=c("#0F084B", "slategrey"), labels = c("Shape only", "With HPO")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")

dev.off()

#no error bars and mixed prevalences####
#simulate prevalences by randomly mixing in hpo and non-hpo preds####
#for each syndrome, for each term 
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
    tmp.prevalence <- sample(1:length(tmp.hpo.pred), length(tmp.hpo.pred) * .5)
    tmp.hpo.pred[tmp.prevalence] <- loocv.pred[hdrda.df$synd == unique(ultimate.bunduru$synd)[i]][tmp.prevalence]
    sim.synd.hpo.results <- c(sim.synd.hpo.results, tmp.hpo.pred)
  }
}
  permuted.results[,k] <- sim.synd.hpo.results
  print(k)
}

#get mean max min sd of sensitivities from 1000 simulated analysis at 50% term prevalence
permuted.sens <- apply(permuted.results, 2, function(x) confusionMatrix(factor(x, levels = levels(hdrda.df$synd)), factor(ultimate.bunduru$synd, levels = levels(hdrda.df$synd)))$byClass[,1])
permuted.mean <- apply(permuted.sens, 1, mean)
permuted.sd  <- apply(permuted.sens, 1, sd)
permuted.max  <- apply(permuted.sens, 1, max)
permuted.min  <- apply(permuted.sens, 1, min)

#bind permuted.results with true synd
hpo.df <- data.frame(synd = levels(hdrda.df$synd), orig.sens = original.sens, mean = permuted.mean, sd = permuted.sd, max = permuted.max, min = permuted.min)

fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/top1comparison_prevalence.pdf", width = 12, height = 8)
ggplot(aes(x = reorder(synd, -orig.sens), y = mean), data = hpo.df) +
  geom_bar(stat = "identity", fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = mean - (2*sd), ymax = mean + (2*sd)),  fill = 'black') + 
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")

dev.off()

#redo sim checking top 3 posterior probabilities####
#simulate prevalences by randomly mixing in hpo and non-hpo preds####
#for each syndrome, for each term 
colnames(hpo.distribution) <- levels(hdrda.df$synd)
ultimate.bunduru <- data.frame(hpo.meta, pred = hpo.pred, hpo.distribution)
colnames(ultimate.bunduru)[1:3] <- c("synd", "hpo.id", "hpo.term")

permuted.results <- matrix(NA, nrow = nrow(ultimate.bunduru), ncol = 1000)

for(k in 1:1000){
  sim.synd.hpo.results <- NULL
  for(i in 1:length(unique(ultimate.bunduru$synd))){
    tmp.terms <- unique(ultimate.bunduru$hpo.term[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
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
      
      tmp.prevalence <- sample(1:length(rank2.check), length(rank2.check) * .5)
      rank2.check[tmp.prevalence] <- loocv.pred[hdrda.df$synd == unique(ultimate.bunduru$synd)[i]][tmp.prevalence]
      sim.synd.hpo.results <- c(sim.synd.hpo.results, rank2.check)
    }
  }
  permuted.results[,k] <- sim.synd.hpo.results
  print(k)
}

#bind permuted.results with true synd
permuted.results.rank2 <- permuted.results
permuted.sens3 <- apply(permuted.results.rank2, 2, function(x) confusionMatrix(factor(x, levels = levels(hdrda.df$synd)), factor(ultimate.bunduru$synd, levels = levels(hdrda.df$synd)))$byClass[,1])
permuted.mean3 <- apply(permuted.sens, 1, mean)
permuted.sd3  <- apply(permuted.sens, 1, sd)
permuted.max3  <- apply(permuted.sens, 1, max)
permuted.min3  <- apply(permuted.sens, 1, min)

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

original.sens.rank2 <-  confusionMatrix(factor(rank2.result, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1]
hpo.df3 <- data.frame(synd = levels(hdrda.df$synd), orig.sens = original.sens.rank2, mean = permuted.mean3, sd = permuted.sd3, max = permuted.max3, min = permuted.min3)

fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/top3comparison_prevalence.pdf", width = 12, height = 8)
ggplot(aes(x = reorder(synd, -orig.sens), y = mean), data = hpo.df3) +
  geom_bar(stat = "identity", fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = mean - (2*sd), ymax = mean + (2*sd)),  fill = 'black') + 
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")

dev.off()


#redo sim checking top 10 posterior probabilities####
#simulate prevalences by randomly mixing in hpo and non-hpo preds####
#for each syndrome, for each term 
colnames(hpo.distribution) <- levels(hdrda.df$synd)
ultimate.bunduru <- data.frame(hpo.meta, pred = hpo.pred, hpo.distribution)
colnames(ultimate.bunduru)[1:3] <- c("synd", "hpo.id", "hpo.term")

permuted.results <- matrix(NA, nrow = nrow(ultimate.bunduru), ncol = 1000)

for(k in 1:1000){
  sim.synd.hpo.results <- NULL
  for(i in 1:length(unique(ultimate.bunduru$synd))){
    tmp.terms <- unique(ultimate.bunduru$hpo.term[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i]])
    for(j in 1:length(unique(tmp.terms))){
      #grab synd i, term j posteriors and check top 3 posteriors. Then simulate .5 hpo term prevalence
      tmp.hpo.post <- ultimate.bunduru[ultimate.bunduru$synd == unique(ultimate.bunduru$synd)[i] & ultimate.bunduru$hpo.term == tmp.terms[j], -1:-4]
      rank3 <- 10
      rank3.check <- rep(NA, nrow(tmp.hpo.post))
      for(l in 1: length(rank3.check)){
        if(sum(levels(hdrda.df$synd)[order(tmp.hpo.post[l,], decreasing = T)][1:rank3] == unique(ultimate.bunduru$synd)[i]) > 0){
          rank3.check[l] <- as.character(unique(ultimate.bunduru$synd)[i])
        } else{rank3.check[l] <- names(sort(tmp.hpo.post[l,], decreasing = T)[1])}
      }
      
      tmp.prevalence <- sample(1:length(rank3.check), length(rank3.check) * .5)
      rank3.check[tmp.prevalence] <- loocv.pred[hdrda.df$synd == unique(ultimate.bunduru$synd)[i]][tmp.prevalence]
      sim.synd.hpo.results <- c(sim.synd.hpo.results, rank3.check)
    }
    
  }
  permuted.results[,k] <- sim.synd.hpo.results
  print(k)
}

#bind permuted.results with true synd
permuted.results.rank3 <- permuted.results
permuted.sens10 <- apply(permuted.results.rank3, 2, function(x) confusionMatrix(factor(x, levels = levels(hdrda.df$synd)), factor(ultimate.bunduru$synd, levels = levels(hdrda.df$synd)))$byClass[,1])
permuted.mean10 <- apply(permuted.sens, 1, mean)
permuted.sd10  <- apply(permuted.sens, 1, sd)
permuted.max10  <- apply(permuted.sens, 1, max)
permuted.min10  <- apply(permuted.sens, 1, min)

#rank3 for original model####
rank3 <- 3
rank3.result <- loocv.pred
for(i in 1:length(unique(hdrda.df$synd))){
  tmp.hpo.post <- loocv.post[hdrda.df$synd == levels(hdrda.df$synd)[i],]
  tmp.hpo.pred <-loocv.pred[hdrda.df$synd == levels(hdrda.df$synd)[i]]
  rank3.check <- rep(NA, nrow(tmp.hpo.post))
  for(l in 1: length(rank3.check)){
    if(sum(levels(hdrda.df$synd)[order(tmp.hpo.post[l,], decreasing = T)][1:rank3] == levels(hdrda.df$synd)[i]) > 0){
      rank3.check[l] <- as.character(levels(hdrda.df$synd)[i])
    } else{rank3.check[l] <- tmp.hpo.pred[l]}
  }
  rank3.result[hdrda.df$synd == levels(hdrda.df$synd)[i]] <- rank3.check
}

original.sens.rank3 <-  confusionMatrix(factor(rank3.result, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1]
hpo.df10 <- data.frame(synd = levels(hdrda.df$synd), orig.sens = original.sens.rank3, mean = permuted.mean10, sd = permuted.sd10, max = permuted.max10, min = permuted.min10)

fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/top10comparison_prevalence.pdf", width = 12, height = 8)
ggplot(aes(x = reorder(synd, -orig.sens), y = mean), data = hpo.df10) +
  geom_bar(stat = "identity", fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = mean - (2*sd), ymax = mean + (2*sd)),  fill = 'black') + 
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")

dev.off()





#how bad is it to supply an incorrect term?####








