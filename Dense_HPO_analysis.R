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
# save(atlas, d.meta, front.face, PC.eigenvectors, PC.scores, synd.mshape, phenotype.df, hpo, hpo.pos, hdrda.mod, hdrda.df, official.names, file = "data.Rdata")
load("data.Rdata")

#when someone selects specific terms, find the omims/syndrome name in hpo.pos
print(hpo.pos[hpo.pos[,5] == names(hpo$name)[hpo$name == "Strabismus"],3])

in.hpo <- rep(NA, length(unique(hdrda.df$synd)))

View(data.frame(official.names))

for(i in 1:length(in.hpo)) in.hpo[i] <- length(grep(official.names[i], x = hpo.pos[hpo.pos[,5] == "HP:0009601",3])) > 0

in.hpo[official.names == "Unaffected Unrelated"] <- TRUE

#priors are adjusted to uniformly sharing 20% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list

updated.priors <- rep(NA, length(official.names))
names(updated.priors) <- official.names

updated.priors[in.hpo == F] <- .05 / length(official.names[in.hpo == F])
updated.priors[in.hpo == T] <- .95 / length(official.names[in.hpo == T])

hdrda.updated <- hdrda(synd ~ ., data = hdrda.df, prior = updated.priors)

posterior.distribution <- predict(hdrda.updated, newdata = hdrda.df[which(d.meta$Syndrome == "Nager Syndrome"),-1], type = "prob")$post[1,]

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


#testing the method with one HPO term at a time, assuming equal prevelance####
#for all the people with syndrome i, let's look at the sensitivity with HPO term j
#how to deal with differing number of terms for each synd? Save only the mean top1,3,10 sens and the min/max, and the number of HPO terms associated

hpo.perf <- matrix(NA, nrow = length(unique(hdrda.df$synd)), ncol = 10)
colnames(hpo.perf) <- c("top1.mean", "top3.mean","top10.mean", "top1.min", "top3.min", "top10.min", "top1.max", "top3.max", "top10.max", "N.hpo")

for(i in 1 : length(unique(hdrda.df$synd))){
  
  hpo.perf[i,10] <- length(unique(hpo.pos[hpo.pos$V3 == official.names[i],5]))
  tmp.means <- rep(NA,  hpo.perf[i,10])
  
  for(j in 1 : hpo.perf[i,10]){
    
    #make vector to store perf for each term
    
    # hdrda.df[hdrda.df$synd == unique(hdrda.df$synd)[i], -1]
    
    for(k in 1:length(in.hpo)) in.hpo[k] <- length(grep(official.names[k], x = hpo.pos[hpo.pos[,5] == unique(hpo.pos[hpo.pos$V3 == official.names[i],5])[j],3])) > 0
    
    in.hpo[official.names == "Non-syndromic"] <- TRUE
    
    #priors are adjusted to uniformly sharing 10% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list
    
    updated.priors <- rep(NA, length(official.names))
    names(updated.priors) <- official.names
    
    updated.priors[in.hpo == F] <- 0.1 / length(official.names[in.hpo == F])
    updated.priors[in.hpo == T] <- 0.9 / length(official.names[in.hpo == T])
    
    updated.priors <- updated.priors/sum(updated.priors)
    
    hdrda.updated <- hdrda(synd ~ ., data = hdrda.df, prior = updated.priors)
    
    posterior.distribution <- predict(hdrda.updated, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$post
    
    posterior.class <- predict(hdrda.updated, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$class
    
    tmp.gs <- hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i], 1]
    levels(tmp.gs) <- levels(posterior.class)
    
    print(paste0(unique(official.names)[i], ": term ", j, " out of ", hpo.perf[i,10], "-- mean: ", round(confusionMatrix(posterior.class, factor(tmp.gs, levels = levels(posterior.class)))$byClass[i,1], digits = 3)))
    
    tmp.means[j] <- confusionMatrix(posterior.class, factor(tmp.gs, levels = levels(posterior.class)))$byClass[i,1]
    
  }
  
  hpo.perf[i,1] <- mean(tmp.means)
  hpo.perf[i,4] <- min(tmp.means)
  hpo.perf[i,7] <- max(tmp.means)
  
  
}

#plot result####
original.sens <- c(.666, .5, .16666, .47, .75, .482, .772, .477, .59, .43, .8, .5, .55, .58, .68, .61, .82, .056, .138, .48, .28, .31, .55,.124,.28,.34,.21,.519,.64,.08,.96,.45,.35,.25,.5,.56,.1, .17, .88, .47, .48, .4,.22,.68,.333,.22,.61,.53,.54,.47,.64,.59)
hpo.df <- data.frame(Syndrome = levels(hdrda.df$synd), orig.sens = original.sens, hpo.perf)

ggplot(aes(x = reorder(Syndrome, -orig.sens), y = top1.mean), data = hpo.df) +
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



#testing the method with one HPO term at a time with simulated term prevalances####
#for all the people with syndrome i, let's look at the sensitivity with HPO term j
#how to deal with differing number of terms for each synd? Save only the mean top1,3,10 sens and the min/max, and the number of HPO terms associated

hpo.perf <- matrix(NA, nrow = length(unique(hdrda.df$synd)), ncol = 10)
colnames(hpo.perf) <- c("top1.mean", "top3.mean","top10.mean", "top1.min", "top3.min", "top10.min", "top1.max", "top3.max", "top10.max", "N.hpo")
hdrda.orig <- hdrda(Syndrome ~ ., data = hdrda.df)

for(i in 1 : length(unique(hdrda.df$synd))){
  
  hpo.perf[i,10] <- length(unique(hpo.pos[hpo.pos$V3 == official.names[i],5]))
  tmp.means <- rep(NA,  hpo.perf[i,10])
  
  for(j in 1 : hpo.perf[i,10]){
    
    #make vector to store perf for each term
    
    hdrda.df[hdrda.df$synd == unique(hdrda.df$synd)[i], -1]
    
    hpo.term <- unique(hpo.pos[hpo.pos$V3 == official.names[i],5])[j]
    
    for(k in 1:length(in.hpo)){
      in.hpo[k] <- length(grep(official.names[k], x = hpo.pos[hpo.pos[,5] == hpo.term,3])) > 0
      
    }
    
    in.hpo[official.names == "Non-syndromic"] <- TRUE

    #get frequency of term with current syndrome####
    
    tmp.hpo.frequency <- phenotype.df$Frequency[phenotype.df$HPO_ID == hpo.term][grep(official.names[i], phenotype.df[phenotype.df$HPO_ID == hpo.term,2], ignore.case = T)]
    
    if(length(tmp.hpo.frequency) == 0){ tmp.hpo.frequency <- 0
    } else if(is.na(tmp.hpo.frequency)){ tmp.hpo.frequency <- 0
    } else if(tmp.hpo.frequency == "HP:0040280"){ tmp.hpo.frequency <- 1
    } else if(tmp.hpo.frequency == "HP:0040281"){ tmp.hpo.frequency <- (80+99/2)/100
    } else if(tmp.hpo.frequency == "HP:0040282"){ tmp.hpo.frequency <- (30+79/2)/100
    } else if(tmp.hpo.frequency == "HP:0040283"){ tmp.hpo.frequency <- (5+29/2)/100
    } else if(tmp.hpo.frequency == "HP:0040284"){ tmp.hpo.frequency <- (1+4/2)/100
    } else if(tmp.hpo.frequency == "HP:0040285"){ tmp.hpo.frequency <- 0
    } else if(length(grep("%", tmp.hpo.frequency)) > 0){tmp.hpo.frequency <- as.numeric(substr(tmp.hpo.frequency, 1, nchar(tmp.hpo.frequency)-1))/100
    } else if(length(grep("/", tmp.hpo.frequency)) > 0){tmp.hpo.frequency <- eval(parse(text = tmp.hpo.frequency))
    }
    
    if(tmp.hpo.frequency > 0){
      #priors are adjusted to uniformly sharing 20% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list
      
      updated.priors <- rep(NA, length(official.names))
      names(updated.priors) <- official.names
      
      updated.priors[in.hpo == F] <- 0.1 / length(official.names[in.hpo == F])
      updated.priors[in.hpo == T] <- 0.9 / length(official.names[in.hpo == T])
      
      updated.priors <- updated.priors/sum(updated.priors)
      
      hdrda.updated <- hdrda(Syndrome ~ ., data = hdrda.df, prior = updated.priors)
      hdrda.df$synd[hdrda.df$synd == "Control"] <- "Non-syndromic"
      
      
      #calculate an index of which observations get the HPO term boost
      synd.count <- 1:length(hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],2])
      updated.posterior.sample <- sample(synd.count, length(synd.count)*tmp.hpo.frequency)
      
      
      posterior.distribution <- predict(hdrda.updated, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$post
      orig.posterior.distribution <- predict(hdrda.cva, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")
      
      posterior.class <- predict(hdrda.updated, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$class
      levels.to.keep <- levels(posterior.class)
      old.posterior.class <- predict(hdrda.orig, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")$class
      
      posterior.class <- as.character(posterior.class)
      old.posterior.class <- as.character(old.posterior.class)
      # orig.posterior.class <- as.character(predict(hdrda.cva, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1]))
      # orig.posterior.class[orig.posterior.class == "Control"] <- "Non-syndromic"
      # orig.posterior.class <- factor(orig.posterior.class, levels = levels(posterior.class))
      
      #merge simulated hpo data w/original data
      posterior.distribution[-updated.posterior.sample,] <- orig.posterior.distribution[-updated.posterior.sample,]
      posterior.class[-updated.posterior.sample] <- old.posterior.class[-updated.posterior.sample]
      
      posterior.class <- factor(posterior.class, levels = levels.to.keep)
      tmp.gs <- factor((hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i], 2]), levels = levels(posterior.class))
    } else {
      
      posterior.distribution <- predict(hdrda.cva, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1], type = "prob")
      posterior.class <- predict(hdrda.cva, newdata = hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i],-1])
      tmp.gs <- factor(make.names(hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i], 2]), levels = levels(posterior.class))
    }
    
    # tmp.gs <- as.factor(hdrda.df[hdrda.df$synd == levels(hdrda.df$synd)[i], 2])
    # levels(tmp.gs) <- levels(hdrda.df$synd)
    
    
    tmp.means[j] <- confusionMatrix(posterior.class, tmp.gs)$byClass[which(levels(tmp.gs) == tmp.gs[1]),1]
    
    print(paste0(unique(official.names)[i], ": term ", j, " out of ", hpo.perf[i,10], "-- mean: ", tmp.means[j]))
    
    
    
  }
  
  hpo.perf[i,1] <- mean(tmp.means)
  hpo.perf[i,4] <- min(tmp.means)
  hpo.perf[i,7] <- max(tmp.means)
  
  
}

#plot result priors with simulated term prevalances####
original.sens <- c(.666, .5, .16666, .47, .75, .482, .772, .477, .59, .43, .8, .5, .55, .58, .68, .61, .82, .056, .138, .48, .28, .31, .55,.124,.28,.34,.21,.519,.64,.08,.96,.45,.35,.25,.5,.56,.1, .17, .88, .47, .48, .4,.22,.68,.333,.22,.61,.53,.54,.47,.64,.59)
hpo.df <- data.frame(Syndrome = levels(hdrda.df$synd), orig.sens = original.sens, hpo.perf)

ggplot(aes(x = reorder(Syndrome, -orig.sens), y = top1.mean), data = hpo.df) +
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














