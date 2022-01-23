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

# orig.preds <- (predict(hdrda.mod, hdrda.df[,-1]))
# View(confusionMatrix(orig.preds$class, as.factor(hdrda.df[,1]))$byClass)

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

View(standardized.freqs)

#NAs are defined as obligate####
standardized.freqs$Frequency[is.na(standardized.freqs$Frequency)] <- .545

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
    if(length(tmp.hpo.frequency) == 0) tmp.hpo.frequency <- .545
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
}

# save(synd.hpo.result, file = "hpo_results_NA_54_4.Rdata")

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
original.sens <- confusionMatrix(hdrda.orig.preds, as.factor(hdrda.df$synd))$byClass[,1]

hpo.df <- data.frame(orig.sens = original.sens[match(hpo.perf_1$synd, levels(hdrda.df$synd))], synd = hpo.perf_1$synd, mean1 = hpo.perf_1$top1.mean, mean5 = hpo.perf_5$top1.mean, mean25 = hpo.perf_25$top1.mean)

ggplot(aes(x = reorder(synd, -orig.sens), y = mean1), data = hpo.df) +
  geom_bar(stat = "identity",  fill = "#0F084B") + 
  geom_bar(stat = "identity", aes(y = mean5), fill = "#6066AC", alpha = 1) +
  geom_bar(stat = "identity", aes(y = mean25), fill = "#B0B5EF", alpha = 1) +
  geom_bar(stat = "identity", aes(y = orig.sens), fill = "#DFE1F9") +
  ylab("Sensitivity") +
  xlab("Syndrome") +
  # scale_colour_manual(name = "Classification /n approach", values=c("#0F084B", "slategrey"), labels = c("Shape only", "With HPO")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 9),
        plot.background = element_rect(fill = "transparent"),
        legend.position = "none")











