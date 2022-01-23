library(sparsediscrim)
library(caret)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

loocv.pred <- rep(NA, nrow(hdrda.df))
loocv.post <- matrix(NA, nrow = nrow(hdrda.df), ncol = length(unique(hdrda.df$synd)))
hdrda.df$synd <- as.factor(hdrda.df$synd)

for(i in 1:nrow(hdrda.df)){
loocv.mod <- hdrda(synd ~ ., data = hdrda.df[-i,])
tmp.pred <- predict(loocv.mod, hdrda.df[i,-1])
loocv.pred[i] <- as.character(tmp.pred$class)
loocv.post[i,] <- tmp.pred$post
print(paste0(i, "th pred: ",loocv.pred[i], ". True synd: ", hdrda.df$synd[i]))
save(loocv.pred, loocv.post, file = "loocv_hdrda.Rdata")
}

