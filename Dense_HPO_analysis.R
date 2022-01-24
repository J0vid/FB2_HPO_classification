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


#when someone selects specific terms, find the omims/syndrome name in hpo.pos
print(hpo.pos[hpo.pos[,5] == names(hpo$name)[hpo$name == "Strabismus"],3])

in.hpo <- rep(NA, length(unique(hdrda.df$synd)))

View(data.frame(official.names))

for(i in 1:length(in.hpo)) in.hpo[i] <- length(grep(official.names[i], x = hpo.pos[hpo.pos[,5] == "HP:0009601",3])) > 0

in.hpo[official.names == "Non-syndromic"] <- TRUE

#priors are adjusted to uniformly sharing 20% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list

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

#most frequencies are blank--what should we do? We could limit the analysis to only the syndromes and terms with frequency stats or we could simulate frequency...
# phenotype.df.synd <- phenotype.df[phenotype.df$DiseaseName %in% official.names,]
View(phenotype.df.synd[is.na(phenotype.df.synd$Frequency) == F,])

#outside of one-off examples, I have moved the all of the heavy computation out of this script and into job scripts. Visualization continues below.
load("hpo_results_NA_54_4.Rdata")
load("full_hpo_results_NA_545.Rdata")


#plot result priors with simulated term prevalances####
hpo.perf <- synd.hpo.result %>%
  group_by(synd) %>%
  summarise(top1.min = min(sensitivity), top1.max = max(sensitivity), top1.mean = mean(sensitivity))

hdrda.orig.preds <- predict(hdrda.mod, newdata = hdrda.df[,-1])$class
original.sens <- confusionMatrix(hdrda.orig.preds, hdrda.df$synd)$byClass[,1]

hpo.df <- data.frame(orig.sens = original.sens[match(hpo.perf$synd, levels(hdrda.df$synd))], hpo.perf)

ggplot(aes(x = reorder(synd, -orig.sens), y = top1.mean), data = hpo.df) +
  geom_bar(stat = "identity",  fill = "slategrey") + 
  geom_errorbar(aes(ymin = top1.min, ymax = top1.max), width=.2,
                position=position_dodge(.9)) +
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#0F084B") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  # scale_colour_manual(name = "Classification /n approach", values=c("#0F084B", "slategrey"), labels = c("Shape only", "With HPO")) +
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



#how bad is it to supply an incorrect term?####








