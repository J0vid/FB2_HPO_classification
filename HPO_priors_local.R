#hpo tests
library(ggplot2)
library(Morpho)
library(geomorph)
library(rgl)
library(dplyr)
library(Rvcg)
library(sparsediscrim)
library(caret)


load("~/shiny/shinyapps/Classification_demo/demo_objects.Rdata")
library(readr)
phenotype_2 <- read_delim("Downloads/phenotype-2.hpoa", 
                          "\t", escape_double = FALSE, trim_ws = TRUE, 
                          skip = 4)

phenotype.df <- as.data.frame(phenotype_2)

max.sample <- read.table("/Volumes/Storage/Hallgrimsson/Collaborations/FaceBase/FaceBase_2/Analysis/2019/Analysis Data/03 Max Synd Sample/PR data sym corrected 20yo.txt")
max.class <- read.table("/Volumes/Storage/Hallgrimsson/Collaborations/FaceBase/FaceBase_2/Analysis/2019/Analysis Data/03 Max Synd Sample/classifiers.txt")

FB.lms2 <- data.frame(Syndrome = max.class$Syndrome, lms = max.sample)

filtered.lms.max <- as.data.frame(FB.lms2 %>% group_by(Syndrome) %>% filter(n()> 2))
filtered.lms.max$Syndrome <- as.character(filtered.lms.max$Syndrome)
filtered.lms.max$Syndrome <- as.factor(filtered.lms.max$Syndrome)
View(table(filtered.lms.max$Syndrome))

updated.priors <- table(filtered.lms$Syndrome)

#zero out known impossibilities based on HPO
#when someone selects specific terms, find the omims/syndrome name in hpo_pos
#here's an example with facial asymmetry selected: hpo_pos[hpo_pos[,5] == "HP:0000322",]
print(hpo_pos[hpo_pos[,5] == names(hpo$name)[hpo$name == "Strabismus"],3])

in.hpo <- rep(NA, length(unique(filtered.lms$Syndrome)))

cva.data <- data.frame(ID = filtered.lms$ID, Syndrome = as.character(filtered.lms$Syndrome), scores = bad.cva$CVscores)
mod.data <- data.frame(Syndrome = as.character(filtered.lms$Syndrome), scores = bad.cva$CVscores)

official.names <- levels(mod.data$Syndrome)

# filtered.lms$Syndrome[filtered.lms$Syndrome == "Cleft_Lip_Palate"] <- "CLEFT LIP/PALATE"
#no 18p tetrasomy, goldenhar in hpo
official.names[2] <- "CHROMOSOME 1P36 DELETION SYNDROME"
official.names[3] <- "CHROMOSOME 22q11.2 DELETION SYNDROME, DISTAL"
official.names[4] <- "CRI-DU-CHAT SYNDROME"
official.names[5] <- "#100800 ACHONDROPLASIA; ACH"
official.names[6] <- "ANGELMAN SYNDROME"
official.names[7] <- "#101200 APERT SYNDROME;;ACROCEPHALOSYNDACTYLY, TYPE I; ACS1;;ACS IAPERT-CROUZON DISEASE, INCLUDED;;ACROCEPHALOSYNDACTYLY, TYPE II, INCLUDED;;ACS II, INCLUDED;;VOGT CEPHALODACTYLY, INCLUDED"
official.names[8] <- "CARDIOFACIOCUTANEOUS SYNDROME 1; CFC1"
official.names[9] <- "CHARGE SYNDROME"
official.names[10] <- "CLEFT LIP/PALATE WITH ABNORMAL THUMBS AND MICROCEPHALY"
official.names[11] <- "COCKAYNE SYNDROME A; CSA"
official.names[12] <- "COFFIN-SIRIS SYNDROME 1; CSS1"
official.names[13] <- "COHEN SYNDROME"
official.names[14] <- "CORNELIA DE LANGE SYNDROME 1"
official.names[15] <- "COSTELLO SYNDROME; CSTLO"
official.names[16] <- "CROUZON SYNDROME"
official.names[17] <- "DOWN SYNDROMETRISOMY 21"
official.names[18] <- "EHLERS-DANLOS SYNDROME, TYPE I"
official.names[19] <- "FRAGILE X MENTAL RETARDATION SYNDROME"
official.names[21] <- "JACOBSEN SYNDROME"
official.names[22] <- "JOUBERT SYNDROME 1"
official.names[23] <- "KABUKI SYNDROME 1"
official.names[25] <- "LOEYS-DIETZ SYNDROME, TYPE 1A LOEYS-DIETZ AORTIC ANEURYSM SYNDROME"
official.names[26] <- "MARFAN SYNDROME; MFS"
official.names[27] <- "MOEBIUS SYNDROME; MBS"
official.names[28] <- "MUCOPOLYSACCHARIDOSIS TYPE IIIA"
official.names[29] <- "ACROFACIAL DYSOSTOSIS 1, NAGER TYPE; AFD1"
official.names[30] <- "NEUROFIBROMATOSIS, TYPE I"
official.names[32] <- "#163950 NOONAN SYNDROME 1; NS1;;NOONAN SYNDROME;;MALE TURNER SYNDROME;;FEMALE PSEUDO-TURNER SYNDROME;;TURNER PHENOTYPE WITH NORMAL KARYOTYPEPTERYGIUM COLLI SYNDROME, INCLUDED"
official.names[33] <- "OSTEOGENESIS IMPERFECTA, TYPE I"
official.names[34] <- "PHELAN-MCDERMID SYNDROME; PHMDS"
official.names[35] <- "311895 PIERRE ROBIN SEQUENCE WITH FACIAL AND DIGITAL ANOMALIES"
official.names[36] <- "PITT-HOPKINS SYNDROME"
official.names[37] <- "PSEUDOACHONDROPLASIA"
official.names[38] <- "RETT SYNDROME; RTT"
official.names[39] <- "#215100 RHIZOMELIC CHONDRODYSPLASIA PUNCTATA, TYPE 1; RCDP1;;PEROXISOME BIOGENESIS DISORDER 9; PBD9;;CHONDRODYSPLASIA PUNCTATA, RHIZOMELIC FORM; CDPR;;CHONDRODYSTROPHIA CALCIFICANS PUNCTATA"
official.names[40] <- "RUBINSTEIN-TAYBI SYNDROME 1; RSTS1"
official.names[41] <- "SILVER-RUSSELL SYNDROME; SRS"
official.names[42] <- "SMITH-LEMLI-OPITZ SYNDROME; SLOS"
official.names[43] <- "SMITH-MAGENIS SYNDROME; SMS"
official.names[44] <- "#117550 SOTOS SYNDROME 1; SOTOS1;;SOTOS SYNDROME;;CEREBRAL GIGANTISM;;CHROMOSOME 5q35 DELETION SYNDROME"
official.names[45] <- "SPONDYLOEPIPHYSEAL DYSPLASIA WITH CONGENITAL JOINT DISLOCATIONS"
official.names[46] <- "STICKLER SYNDROME, TYPE I"
official.names[47] <- "TREACHER COLLINS-FRANCESCHETTI SYNDROME"
official.names[48] <- "TRISOMY 18-LIKE SYNDROME"
official.names[49] <- "MENTAL RETARDATION, X-LINKED, SYNDROMIC, TURNER TYPE; MRXST"
official.names[50] <- "VAN DER WOUDE SYNDROME"
official.names[51] <- "WILLIAMS-BEUREN SYNDROME; WBS"
official.names[52] <- "ECTODERMAL DYSPLASIA 1, HYPOHIDROTIC, X-LINKED; XHED"



for(i in 1:length(in.hpo)) in.hpo[i] <- length(grep(official.names[i], x = hpo_pos[hpo_pos[,5] == "HP:0009601",3])) > 0

in.hpo[official.names == "Non-syndromic"] <- TRUE

# #old code that zeroed out priors
# mod.data <- cva.data[filtered.lms$Syndrome %in% official.names[in.hpo == T],-1]
# mod.data$Syndrome <- factor(mod.data$Syndrome, levels = official.names[in.hpo == T])

#priors are adjusted to uniformly sharing 20% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list

updated.priors <- rep(NA, length(official.names))
names(updated.priors) <- official.names

updated.priors[in.hpo == F] <- .05 / length(official.names[in.hpo == F])
updated.priors[in.hpo == T] <- .95 / length(official.names[in.hpo == T])

hdrda.updated <- hdrda(Syndrome ~ ., data = mod.data, prior = updated.priors)

posterior.distribution <- predict(hdrda.updated, newdata = rbind(cva.data[cva.data$ID == "160728105214",-1:-2],cva.data[cva.data$ID == "160728105214",-1:-2]), type = "prob")$post[1,]

posterior.distribution <- sort(posterior.distribution, decreasing = T)

#used to be part of plot.df: ID = as.factor(1:10), 
plot.df <- data.frame(Probs = round(as.numeric(posterior.distribution[1:10]), digits = 4), Syndrome = as.factor(gsub("_", " ", names(posterior.distribution[1:10]))))
plot.df$Syndrome <- as.character(plot.df$Syndrome)
plot.df$Syndrome[plot.df$Syndrome == "Control"] <- "Non-syndromic"


plot_ly(data = plot.df, x = ~Syndrome, y = ~Probs, type = "bar", color = I("grey"), hoverinfo = paste0("Syndrome: ", "x", "<br>", "Probability: ", "y")) %>%
  layout(xaxis = list(tickvals = gsub("_", " ", plot.df$Syndrome), tickangle = 45, ticktext = c(Syndrome = plot.df$Syndrome, Probability = plot.df$Probs), title = "<b>Syndrome</b>"),
         yaxis = list(title = "<b>Class probability</b>"),
         paper_bgcolor='rgba(245, 245, 245, .9)',
         margin = list(b = 125, l = 50, r = 100)
  )


#testing the method with one HPO term at a time####
#for all the people with syndrome i, let's look at the sensitivity with HPO term j
#how to deal with differing number of terms for each synd? Save only the mean top1,3,10 sens and the min/max, and the number of HPO terms associated

hpo.perf <- matrix(NA, nrow = length(unique(filtered.lms$Syndrome)), ncol = 10)
colnames(hpo.perf) <- c("top1.mean", "top3.mean","top10.mean", "top1.min", "top3.min", "top10.min", "top1.max", "top3.max", "top10.max", "N.hpo")

for(i in 1 : length(unique(filtered.lms$Syndrome))){
  
 
  
  hpo.perf[i,10] <- length(unique(hpo_pos[hpo_pos$V3 == official.names[i],5]))
  tmp.means <- rep(NA,  hpo.perf[i,10])
  
for(j in 1 : hpo.perf[i,10]){
  
  #make vector to store perf for each term

filtered.lms[filtered.lms$Syndrome == unique(filtered.lms$Syndrome)[i], -1:-2]

for(k in 1:length(in.hpo)) in.hpo[k] <- length(grep(official.names[k], x = hpo_pos[hpo_pos[,5] == unique(hpo_pos[hpo_pos$V3 == official.names[i],5])[j],3])) > 0

in.hpo[official.names == "Non-syndromic"] <- TRUE

# #old code that zeroed out priors
# mod.data <- cva.data[filtered.lms$Syndrome %in% official.names[in.hpo == T],-1]
# mod.data$Syndrome <- factor(mod.data$Syndrome, levels = official.names[in.hpo == T])

#priors are adjusted to uniformly sharing 20% for each syndrome not in the selected HPO, the rest of the weight is uniform with the remaining syndromes in the HPO list

updated.priors <- rep(NA, length(official.names))
names(updated.priors) <- official.names

updated.priors[in.hpo == F] <- 0.1 / length(official.names[in.hpo == F])
updated.priors[in.hpo == T] <- 0.9 / length(official.names[in.hpo == T])

updated.priors <- updated.priors/sum(updated.priors)

hdrda.updated <- hdrda(Syndrome ~ ., data = mod.data, prior = updated.priors)

posterior.distribution <- predict(hdrda.updated, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2], type = "prob")$post

posterior.class <- predict(hdrda.updated, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2], type = "prob")$class

tmp.gs <- filtered.lms[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i], 2]
levels(tmp.gs) <- levels(posterior.class)

print(paste0(unique(official.names)[i], ": term ", j, " out of ", hpo.perf[i,10], "-- mean: ", confusionMatrix(posterior.class, factor(tmp.gs, levels = levels(posterior.class)))$byClass[i,1]))

tmp.means[j] <- confusionMatrix(posterior.class, factor(tmp.gs, levels = levels(posterior.class)))$byClass[i,1]

}
  
  hpo.perf[i,1] <- mean(tmp.means)
  hpo.perf[i,4] <- min(tmp.means)
  hpo.perf[i,7] <- max(tmp.means)
  
  
}

#plot result####
original.sens <- c(.666, .5, .16666, .47, .75, .482, .772, .477, .59, .43, .8, .5, .55, .58, .68, .61, .82, .056, .138, .48, .28, .31, .55,.124,.28,.34,.21,.519,.64,.08,.96,.45,.35,.25,.5,.56,.1, .17, .88, .47, .48, .4,.22,.68,.333,.22,.61,.53,.54,.47,.64,.59)
hpo.df <- data.frame(Syndrome = levels(mod.data$Syndrome), orig.sens = original.sens, hpo.perf)

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

hpo.perf <- matrix(NA, nrow = length(unique(filtered.lms$Syndrome)), ncol = 10)
colnames(hpo.perf) <- c("top1.mean", "top3.mean","top10.mean", "top1.min", "top3.min", "top10.min", "top1.max", "top3.max", "top10.max", "N.hpo")
hdrda.orig <- hdrda(Syndrome ~ ., data = mod.data)

for(i in 1 : length(unique(filtered.lms$Syndrome))){
  
  hpo.perf[i,10] <- length(unique(hpo_pos[hpo_pos$V3 == official.names[i],5]))
  tmp.means <- rep(NA,  hpo.perf[i,10])
  
  for(j in 1 : hpo.perf[i,10]){
    
    #make vector to store perf for each term
    
    filtered.lms[filtered.lms$Syndrome == unique(filtered.lms$Syndrome)[i], -1:-2]
    
    hpo.term <- unique(hpo_pos[hpo_pos$V3 == official.names[i],5])[j]
    
    for(k in 1:length(in.hpo)){
      in.hpo[k] <- length(grep(official.names[k], x = hpo_pos[hpo_pos[,5] == hpo.term,3])) > 0
      
    }
    
    in.hpo[official.names == "Non-syndromic"] <- TRUE
    
    # #old code that zeroed out priors
    # mod.data <- cva.data[filtered.lms$Syndrome %in% official.names[in.hpo == T],-1]
    # mod.data$Syndrome <- factor(mod.data$Syndrome, levels = official.names[in.hpo == T])
    
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
    
    hdrda.updated <- hdrda(Syndrome ~ ., data = mod.data, prior = updated.priors)
    mod.data$Syndrome[mod.data$Syndrome == "Control"] <- "Non-syndromic"
    
    
    #calculate an index of which observations get the HPO term boost
    synd.count <- 1:length(cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],2])
    updated.posterior.sample <- sample(synd.count, length(synd.count)*tmp.hpo.frequency)
    
    
    posterior.distribution <- predict(hdrda.updated, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2], type = "prob")$post
    orig.posterior.distribution <- predict(hdrda.cva, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2], type = "prob")
    
    posterior.class <- predict(hdrda.updated, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2], type = "prob")$class
    levels.to.keep <- levels(posterior.class)
    old.posterior.class <- predict(hdrda.orig, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2], type = "prob")$class
    
    posterior.class <- as.character(posterior.class)
    old.posterior.class <- as.character(old.posterior.class)
    # orig.posterior.class <- as.character(predict(hdrda.cva, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2]))
    # orig.posterior.class[orig.posterior.class == "Control"] <- "Non-syndromic"
    # orig.posterior.class <- factor(orig.posterior.class, levels = levels(posterior.class))
    
    #merge simulated hpo data w/original data
    posterior.distribution[-updated.posterior.sample,] <- orig.posterior.distribution[-updated.posterior.sample,]
    posterior.class[-updated.posterior.sample] <- old.posterior.class[-updated.posterior.sample]
    
    posterior.class <- factor(posterior.class, levels = levels.to.keep)
    tmp.gs <- factor((filtered.lms[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i], 2]), levels = levels(posterior.class))
    } else {
      
      posterior.distribution <- predict(hdrda.cva, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2], type = "prob")
      posterior.class <- predict(hdrda.cva, newdata = cva.data[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i],-1:-2])
      tmp.gs <- factor(make.names(filtered.lms[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i], 2]), levels = levels(posterior.class))
    }
    
    # tmp.gs <- as.factor(filtered.lms[filtered.lms$Syndrome == levels(mod.data$Syndrome)[i], 2])
    # levels(tmp.gs) <- levels(mod.data$Syndrome)
    
    
    tmp.means[j] <- confusionMatrix(posterior.class, tmp.gs)$byClass[which(levels(tmp.gs) == tmp.gs[1]),1]
    
    print(paste0(unique(official.names)[i], ": term ", j, " out of ", hpo.perf[i,10], "-- mean: ", tmp.means[j]))
    
    
    
  }
  
  hpo.perf[i,1] <- mean(tmp.means)
  hpo.perf[i,4] <- min(tmp.means)
  hpo.perf[i,7] <- max(tmp.means)
  
  
}

#plot result priors with simulated term prevalances####
original.sens <- c(.666, .5, .16666, .47, .75, .482, .772, .477, .59, .43, .8, .5, .55, .58, .68, .61, .82, .056, .138, .48, .28, .31, .55,.124,.28,.34,.21,.519,.64,.08,.96,.45,.35,.25,.5,.56,.1, .17, .88, .47, .48, .4,.22,.68,.333,.22,.61,.53,.54,.47,.64,.59)
hpo.df <- data.frame(Syndrome = levels(mod.data$Syndrome), orig.sens = original.sens, hpo.perf)

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









