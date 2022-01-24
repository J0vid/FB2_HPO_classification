library(Morpho)

load("D:/hpo_data/30k.Rdata")
d.meta.ns <- d.meta

tmp.mesh <- file2mesh("whoami3.ply")

sample1k <- sample(1:27903, 1000)
rot.ns <- array(NA, dim = c(27903, 3, dim(d.lms.combined)[3]))
tmp.mesh <- atlas
for(i in 1:dim(rot.ns)[3]){
  
  tmp.mesh$vb[-4,] <- t(d.lms.combined[,,i])  
  
  trafo1k <- rotmesh.onto(tmp.mesh, d.lms.combined[sample1k,,i], d.registered$mshape[sample1k,], scale = T)
  rot.ns[,,i] <- t(trafo1k$mesh$vb[-4,])
  print(i)
}

#combine datasets
d.lms <- abind::abind(rot.ns, d.registered$rotated, along = 3)

colnames(d.meta.ns)

ns.meta <- cbind(d.meta.ns[-1913,c(1,3,6)], rep("Non-syndromic", 2272), d.meta.ns[-1913,2])
colnames(ns.meta)[4:5] <- c("Syndrome", "scan") 

#match sex cov
ns.meta$Sex <- ns.meta$Sex - 1 
d.meta.combined$Sex <- as.numeric(d.meta.combined$Sex == 'F')
d.meta.combined$Syndrome <- as.character(d.meta.combined$Syndrome)
d.meta.combined$Syndrome[d.meta.combined$Syndrome == "Unaffected Unrelated"] <- "Non-syndromic"
d.meta.combined$Syndrome[d.meta.combined$Syndrome == "Osteogenesis imperfecta"] <- "Osteogenesis Imperfecta"
#remove synds
bad.synds <- c("Rett Syndrome_CDKL5", "15q26.3 Deletion", "Mowat-Wilson Syndrome", "XXX", "Facial Dysmorphia, possibly genetic", "Connective Tissue Disorder")
d.lms <- d.lms[,, -which(d.meta.combined$Syndrome %in% bad.synds)] #get rid of bad lms

d.meta.combined <- d.meta.combined[ -which(d.meta.combined$Syndrome %in% bad.synds),]
d.meta <- rbind(d.meta.combined, ns.meta)
#remove bad ns ind reg
tmp.mesh$vb[-4,] <- t(d.lms[,,(nrow(d.meta.combined) + 1913)])  
plot3d(vcgSmooth(tmp.mesh), aspect = "iso", col = "lightgrey", specular = 1)
d.lms <- d.lms[,,-(nrow(d.meta.combined) + 1913)]

View(d.meta)
#cleanup before big mem needed
rm(d.registered)
rm(d.lms.combined)
gc()
#takes forever! d.pca <- prcompfast(vecx(d.lms))

plot(d.pca$x[,1], d.pca$x[,2], col = (d.meta$Syndrome == "Non-syndromic") + 1)

#construct eigenvectors
PC.eigenvectors <- d.pca$rotation[,1:200]
PC.scores <- d.pca$x[,1:200]
synd.mshape <- geomorph::arrayspecs(rbind(d.pca$center,d.pca$center), p = 27903, k = 3)[,,1]

library(sparsediscrim)
#relevel d.meta
d.meta$Syndrome <- factor(d.meta$Syndrome, levels = c("Non-syndromic", unique(d.meta$Syndrome)[unique(d.meta$Syndrome) != "Non-syndromic"]))

hdrda.df <- data.frame(synd = d.meta$Syndrome, PC.scores)
hdrda.mod <- hdrda(synd ~ ., data = hdrda.df)

save(hdrda.df, hdrda.mod, synd.mshape, PC.eigenvectors, PC.scores, d.meta, file = "combined_PCs.Rdata")

#adjust pc scores for effects of age and sex, save resids and coefficients separately
age.sex.lm <- lm(PC.scores ~ d.meta$Sex + poly(d.meta$Age,3))

adjusted.PC.coefs <- age.sex.lm$coefficients

adjusted.PC.scores <- age.sex.lm$residuals

hdrda.df <- data.frame(synd = d.meta$Syndrome, adjusted.PC.scores)
hdrda.mod <- hdrda(synd ~ ., data = hdrda.df)

save(hdrda.df, hdrda.mod, age.sex.lm, adjusted.PC.scores, file = "adjusted_PCs.Rdata")

plot(d.pca$x[,1], d.pca$x[,2])
points(adjusted.PC.scores[,1], adjusted.PC.scores[,2], col = 2)




