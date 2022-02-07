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
#delta plot main figs####
for(j in 1:length(unique(hpo.meta[,1]))){
  #subtract no term from every term
  selected.synd <- unique(hpo.meta[,1])[j]

  hist.df <- data.frame(term = hpo.meta[hpo.meta[,1] == selected.synd,3], X1 =  hpo.distribution[hpo.meta[,1] == selected.synd,])
  no.term.df <- data.frame(X1 = loocv.post[hdrda.df$synd == selected.synd, ])
  
  delta.df <- NULL
  for(i in 1:length(unique(hist.df$term))) delta.df <- rbind(delta.df, hist.df[hist.df$term == unique(hist.df$term)[i], -1] - no.term.df[,])
  
  pred.colors <- c(paste0("Not ", selected.synd), selected.synd)[1 + (hpo.pred[hpo.meta[,1] == selected.synd] == selected.synd)]
  pred.colors <- factor(pred.colors, levels = c(selected.synd, paste0("Not ", selected.synd)))
  delta.df <- data.frame(term = hist.df$term, posterior = delta.df[, which(levels(hdrda.df$synd) == selected.synd)], correct.pred = pred.colors)
  
  pdf(paste0("results/ind_diffs/individual_term_changes_", gsub(pattern = "/", replacement = "_", selected.synd), ".pdf"), width = 8.5, height = 4.5)
  p <- ggplot(aes(x = term, y = posterior), data = delta.df) +
    geom_jitter(aes(colour = correct.pred), alpha = 1, cex = 4, width = 0.25) +
    scale_colour_manual(name = "Prediction", values = c("red", "black")) +
    ylab(paste0("")) +
    xlab("") +
    ylim(-.5,1) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 25, hjust = .5, vjust = .55, size = 17),
          axis.text.y = element_text(size = 17),
          axis.title.x = element_text(size = 17, face = "bold"),
          axis.title.y = element_text(size = 17, face = "bold"),
          legend.position = "none")
  print(p)
  dev.off()
  
  #x/ylabs for supplemental:
  # ggtitle(paste0("Change in posterior values for ", selected.synd, "\n using HPO")) +
  # ylab(paste0("Change in posterior")) +
  #   xlab("HPO term") +
  # theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 17),
  #       axis.text.y = element_text(size = 17),
  #       axis.title.x = element_text(size = 17, face = "bold"),
  #       axis.title.y = element_text(size = 17, face = "bold"),
  #       legend.position = "top")
  
  print(j)
}

#delta plot supplemental fig####
for(j in 1:length(unique(hpo.meta[,1]))){
  #subtract no term from every term
  selected.synd <- unique(hpo.meta[,1])[j]
  
  hist.df <- data.frame(term = hpo.meta[hpo.meta[,1] == selected.synd,3], X1 =  hpo.distribution[hpo.meta[,1] == selected.synd,])
  no.term.df <- data.frame(X1 = loocv.post[hdrda.df$synd == selected.synd, ])
  
  delta.df <- NULL
  for(i in 1:length(unique(hist.df$term))) delta.df <- rbind(delta.df, hist.df[hist.df$term == unique(hist.df$term)[i], -1] - no.term.df[,])
  
  pred.colors <- c(paste0("Not ", selected.synd), selected.synd)[1 + (hpo.pred[hpo.meta[,1] == selected.synd] == selected.synd)]
  pred.colors <- factor(pred.colors, levels = c(selected.synd, paste0("Not ", selected.synd)))
  delta.df <- data.frame(term = hist.df$term, posterior = delta.df[, which(levels(hdrda.df$synd) == selected.synd)], correct.pred = pred.colors)
  
  # pdf(paste0("results/ind_diffs/individual_term_changes_", gsub(pattern = "/", replacement = "_", selected.synd), ".pdf"), width = 8.5, height = 4.5)
  assign(paste0("p", j), 
         ggplotGrob(ggplot(aes(x = term, y = posterior), data = delta.df) +
          geom_jitter(aes(colour = correct.pred), alpha = 1, cex = 2.5, width = 0.25) +
          scale_colour_manual(name = "Prediction", values = c("red", "black")) +
          ylab(paste0("")) +
          xlab("") +
          ylim(-.5,1) +
          ggtitle(selected.synd) +
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 25, hjust = .5, vjust = .55, size = 10.5),
                axis.text.y = element_text(size = 14),
                axis.title.x = element_text(size = 15, face = "bold"),
                axis.title.y = element_text(size = 15, face = "bold"),
                legend.position = "none")
        )
  )
  # print(p)
  # dev.off()
  
  #x/ylabs for supplemental:
  # ggtitle(paste0("Change in posterior values for ", selected.synd, "\n using HPO")) +
  # ylab(paste0("Change in posterior")) +
  #   xlab("HPO term") +
  # theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 17),
  #       axis.text.y = element_text(size = 17),
  #       axis.title.x = element_text(size = 17, face = "bold"),
  #       axis.title.y = element_text(size = 17, face = "bold"),
  #       legend.position = "top")
  
  print(j)
}

#divide delta posterior results into 7 pages####
library(grid)
pdf("results/supp1_pg1.pdf", height = 40, width = 8)
p <- get(paste0("p", 1))
for(i in 2:10) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()

pdf("results/supp1_pg2.pdf", height = 40, width = 8)
p <- get(paste0("p", 11))
for(i in 12:20) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()

pdf("results/supp1_pg3.pdf", height = 40, width = 8)
p <- get(paste0("p", 21))
for(i in 22:30) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()

pdf("results/supp1_pg4.pdf", height = 40, width = 8)
p <- get(paste0("p", 31))
for(i in 32:40) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()


pdf("results/supp1_pg5.pdf", height = 40, width = 8)
p <- get(paste0("p", 41))
for(i in 42:50) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()

pdf("results/supp1_pg6.pdf", height = 40, width = 8)
p <- get(paste0("p", 51))
for(i in 52:60) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()

pdf("results/supp1_pg7.pdf", height = 40, width = 8)
p <- get(paste0("p", 61))
for(i in 62:70) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()

pdf("results/supp1_pg8.pdf", height = 40, width = 8)
p <- get(paste0("p", 71))
for(i in 72:79) p <- rbind(p, get(paste0("p", i)), size = "first")
p$widths <- unit.pmax(p1$widths)
grid.newpage()
grid.draw(p)
dev.off()

#for a given syndrome, how much does the term change top 1,5,10 preds####
#compare top 1 sensitivity####
#error bars####
hpo.perf <- synd.hpo.result %>%
  group_by(synd) %>%
  summarise(top1.min = min(sensitivity), top1.max = max(sensitivity), top1.mean = mean(sensitivity))


original.sens <- confusionMatrix(factor(loocv.pred, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1]
hpo.perf.full <- hpo.perf[match(levels(hdrda.df$synd), as.data.frame(hpo.perf[,1])$synd), ]
hpo.perf.full[,1] <- publication.synd.names


hpo.df <- data.frame(orig.sens = original.sens, hpo.perf.full)

fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/80PC_top1comparison_full.pdf", width = 13, height = 7)
ggplot(aes(x = reorder(synd, -orig.sens), y = top1.mean), data = hpo.df) +
  geom_bar(stat = "identity",  fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = top1.min, ymax = top1.max)) +
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8.5),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  publication.theme

dev.off()

#error bars and mixed prevalences####
#get mean max min sd of sensitivities from 1000 simulated analysis at 50% term prevalence
load("80PC_permuted_results_rank1.Rdata")
permuted.sens <- apply(permuted.results, 2, function(x) confusionMatrix(factor(x, levels = levels(hdrda.df$synd)), factor(ultimate.bunduru$synd, levels = levels(hdrda.df$synd)))$byClass[,1])
permuted.mean <- apply(permuted.sens, 1, mean)
permuted.sd  <- apply(permuted.sens, 1, sd)
permuted.max  <- apply(permuted.sens, 1, max)
permuted.min  <- apply(permuted.sens, 1, min)

#bind permuted.results with true synd
publication.synd.names <- levels(hdrda.df$synd)
publication.synd.names <- gsub(" Syndrome", "", publication.synd.names)
publication.synd.names <- gsub("deletion", "del", publication.synd.names)
publication.synd.names <- gsub("Deletion", "del", publication.synd.names)
publication.synd.names[publication.synd.names == "X-Linked Hypohidrotic Ectodermal Dysplasia"] <- "XLHED"
publication.synd.names[publication.synd.names == "Ectrodactyly-Ectodermal Dysplasia-Cleft Lip/Palate"] <- "EEC"
publication.synd.names[publication.synd.names == "Epileptic Encephalopathy Early Infantile Type 2"] <- "EIEE2"
publication.synd.names[publication.synd.names == "Rhizomelic Chondrodysplasia Punctata"] <- "Rhizo Chond Punct"
publication.synd.names <- gsub("Dysplasia", "dysp", publication.synd.names)

publication.theme <- theme(axis.title = element_text(size = 16, face = "bold"), axis.text.y = element_text(size = 14), axis.text.x =  element_text(size = 12))

hpo.df <- data.frame(synd = publication.synd.names, orig.sens = original.sens, mean = permuted.mean, sd = permuted.sd, max = permuted.max, min = permuted.min)

fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/80_top1comparison_prevalence.pdf", width = 12, height = 6)
ggplot(aes(x = reorder(synd, -orig.sens), y = mean), data = hpo.df) +
  geom_bar(stat = "identity", fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = mean - (2*sd), ymax = mean + (2*sd)),  fill = 'black') + 
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  publication.theme

dev.off()

#redo sim checking top 3 posterior probabilities####
#simulate prevalences by randomly mixing in hpo and non-hpo preds####

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

original.sens.rank2 <-  confusionMatrix(factor(rank2.result, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1]

#bind permuted.results with true synd
permuted.results.rank2 <- permuted.results
load("C:/Users/David A/Downloads/FB2_HPO_classification/80PC_permuted_results_rank2.Rdata")
permuted.sens3 <- apply(permuted.results.rank2, 2, function(x) confusionMatrix(factor(x, levels = levels(hdrda.df$synd)), factor(ultimate.bunduru$synd, levels = levels(hdrda.df$synd)))$byClass[,1])
permuted.mean3 <- apply(permuted.sens3, 1, mean)
permuted.sd3  <- apply(permuted.sens3, 1, sd)
permuted.max3  <- apply(permuted.sens3, 1, max)
permuted.min3  <- apply(permuted.sens3, 1, min)

hpo.df3 <- data.frame(synd = publication.synd.names, orig.sens = original.sens.rank2, mean = permuted.mean3, sd = permuted.sd3, max = permuted.max3, min = permuted.min3)
#for checkpointing: hpo.df3 <- data.frame(synd = levels(hdrda.df$synd), orig.sens = original.sens.rank2, mean = top10_half_prev[,1], sd = permuted.sd3, max = permuted.max3, min = permuted.min3)


fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/80PC_top3comparison_prevalence.pdf", width = 12, height = 6)
ggplot(aes(x = reorder(synd, -orig.sens), y = mean), data = hpo.df3) +
  geom_bar(stat = "identity", fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = mean - (2*sd), ymax = mean + (2*sd)),  fill = 'black') + 
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  publication.theme

dev.off()


#redo sim checking top 10 posterior probabilities####
#simulate prevalences by randomly mixing in hpo and non-hpo preds####
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

original.sens.rank3 <-  confusionMatrix(factor(rank3.result, levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1]

#bind permuted.results with true synd
load("C:/Users/David A/Downloads/FB2_HPO_classification/80PC_permuted_results_rank3.Rdata")
permuted.sens10 <- apply(permuted.results.rank3, 2, function(x) confusionMatrix(factor(x, levels = levels(hdrda.df$synd)), factor(ultimate.bunduru$synd, levels = levels(hdrda.df$synd)))$byClass[,1])
permuted.mean10 <- apply(permuted.sens10, 1, mean)
permuted.sd10  <- apply(permuted.sens10, 1, sd)
permuted.max10  <- apply(permuted.sens10, 1, max)
permuted.min10  <- apply(permuted.sens10, 1, min)

#compare to rank2: View(data.frame(permuted.mean3, permuted.mean10))
hpo.df10 <- data.frame(synd = publication.synd.names, orig.sens = original.sens.rank3, mean = permuted.mean10, sd = permuted.sd10, max = permuted.max10, min = permuted.min10)

#for checkpointing: hpo.df10 <- data.frame(synd = levels(hdrda.df$synd), orig.sens = original.sens.rank3, mean = top10_half_prev[,1], sd = permuted.sd10, max = permuted.max10, min = permuted.min10)

fill <- c("#0F084B", "#3D60A7", "#A0D2E7")

pdf("results/80PC_top10comparison_prevalence.pdf", width = 12, height = 6)
ggplot(aes(x = reorder(synd, -orig.sens), y = mean), data = hpo.df10) +
  geom_bar(stat = "identity", fill = "#3D60A7") + 
  geom_errorbar(aes(ymin = mean - (2*sd), ymax = mean + (2*sd)),  fill = 'black') + 
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none") + 
  publication.theme

dev.off()

#save out all permuted results####
# save(permuted.results, permuted.results.rank2, permuted.results.rank3, file = "updated_permuted_results_ranked_plots.Rdata")

#put it all together####
# hpo.df, hpo.df3, hpo.df10
hpo.df.combo <- data.frame(synd = publication.synd.names, rank1 = hpo.df$mean, rank2 = hpo.df3$mean, rank3 = hpo.df10$mean)
# well need to fuse in the original values for synds without terms
hpo.df.combo$rank1[which(is.na(hpo.df.combo$rank1))] <- hpo.df$orig.sens[which(is.na(hpo.df.combo$rank1))]
hpo.df.combo$rank2[which(is.na(hpo.df.combo$rank2))] <- hpo.df3$orig.sens[which(is.na(hpo.df.combo$rank2))]
hpo.df.combo$rank3[which(is.na(hpo.df.combo$rank3))] <- hpo.df10$orig.sens[which(is.na(hpo.df.combo$rank3))]
View(hpo.df.combo)

pdf("results/80PC_top1_3_10_prevalence.pdf", width = 12, height = 6)
ggplot(aes(x = reorder(synd, -rank1), y = rank3), data = hpo.df.combo) +
  geom_bar(fill =  "#A0D2E7", stat = "identity") +
  geom_bar(aes(x = synd, y = rank2), fill = "#3D60A7", stat = "identity") +
  geom_bar(aes(x = synd, y = rank1), fill = "#0F084B", stat = "identity") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none") +
  publication.theme
dev.off()

#by syndrome means
colMeans(hpo.df.combo[,-1])


#how bad is it to supply an incorrect term?####








