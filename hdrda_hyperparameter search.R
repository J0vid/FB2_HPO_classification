#PC_hyperparam search for HDRDA model

library(sparsediscrim)
library(caret)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# setup 20 folds
kfold.index <- createFolds(hdrda.df$synd, k = 20)
num.pcs <- seq(21, 201, 20)

loocv.pred <- matrix(NA, nrow = nrow(hdrda.df), ncol = length(num.pcs))
loocv.post <- array(NA,  dim = c(nrow(hdrda.df), length(unique(hdrda.df$synd)), length(num.pcs)))
hdrda.df$synd <- as.factor(hdrda.df$synd)

for(j in 1:length(num.pcs)){
  for(i in 1:length(kfold.index)){
    loocv.mod <- hdrda(synd ~ ., data = hdrda.df[-kfold.index[[i]],1:num.pcs[j]])
    tmp.pred <- predict(loocv.mod, hdrda.df[kfold.index[[i]],2:num.pcs[j]])
    loocv.pred[kfold.index[[i]],j] <- as.character(tmp.pred$class)
    loocv.post[kfold.index[[i]],,j] <- tmp.pred$post

  }
  print(j)
  save(loocv.pred, loocv.post, file = "20fold_hdrda_hyperparam_search.Rdata")
}

#what's the best model?
#confusion matrix for each model
kfold.cmat <- matrix(NA, nrow = length(levels(hdrda.df$synd)), ncol = length(num.pcs))

for(i in 1:length(num.pcs)) kfold.cmat[,i] <- (confusionMatrix(factor(loocv.pred[,i], levels = levels(hdrda.df$synd)), hdrda.df[,1])$byClass[,1])

View(kfold.cmat)   

#80 PCs looks best
