#discrimination fig

dist1 <- rnorm(40)

dist1 <- c(rnorm(20, -.5, .75), rnorm(20, .5, .75))
#combination of features
dist2 <- rnorm(40, 0, 1)

#was dist2
# dist1 <- c(rnorm(20, -1, .5), rnorm(20, 1, .5))
png("~/Downloads/lda_points.png", height = 600, width = 600)
plot(dist2 ~ dist1, col = c(rep(rgb(1,0,0,.5), 20), rep(rgb(0,0,0,.5), 20)), pch = 19, xlab = "", ylab = "", cex = 2, xlim = c(-2.3, 2.3), ylim = c(-2,2))
es <- 1.75
segments(mean(dist1[1:20]) - es, mean(dist2[1:20]) - es, mean(dist1[21:40]) + es, mean(dist2[21:40]) + es, lty = 2, lwd = 2.5)
points(c(mean(dist1[1:20]) - es, mean(dist1[21:40]) + es), c(mean(dist2[1:20])-es, mean(dist2[21:40])+es), pch = 19, cex = 2)
dev.off()

#bad discriminator
png("~/Downloads/marginal_bad.png", height = 300, width = 600)
plot(density(dist2), type = "n", ylim = c(0,.7), axes = F, xlab= "", ylab = "", main = "")
polygon(density(dist2[1:20]), col = rgb(1,0,0,.5))
polygon(density(dist2[21:40]), col = rgb(0,0,0,.5))
dev.off()
#good discriminator
png("~/Downloads/marginal_good.png", height = 300, width = 600)
plot(density(dist1), type = "n", ylim = c(0,.7), axes = F, xlab= "", ylab = "", main = "")
polygon(density(dist1[1:20]), col = rgb(1,0,0,.5))
polygon(density(dist1[21:40]), col = rgb(0,0,0,.5))
dev.off()

#best discriminator
projection <- sparsediscrim::hdrda(cbind(dist1, dist2), c(rep(1, 20), rep(2, 20)), lambda = 1)

# plot(density(proj), type = "n", ylim = c(0,1.3), axes = F, xlab= "", ylab = "", main = "")
# polygon(density(proj[1:20]), col = rgb(1,0,0,.3))
# polygon(density(proj[21:40]), col = rgb(0,0,1,.3))


decision.points <- expand.grid(seq(min(dist1), max(dist1), length.out = 50), seq(min(dist2), max(dist2), length.out = 50))

posteriors <- predict(projection, newdata = decision.points)$post

reds <- rep(NA, 2500)

r2b <- colorRampPalette(c("black", "white", "red"))

for(i in 1:2500) reds[i] <- r2b(50)[which.min(abs(posteriors[i,1] - seq(0,1, length.out = 50)))]

png("~/Downloads/lda2.png", height = 600, width = 600)
plot(decision.points[,2] ~ decision.points[,1], col = reds, pch = 19, cex = .3, xlab = "", ylab = "")
# points(dist2 ~ dist1, col = c(rep("red", 20), rep("black", 20)), pch = 19, bg = "black")
# points(decision.points[,2] ~ decision.points[,1], col = blues, pch = 19)
dev.off()


#knn####

train <- cbind(dist1, dist2)
cl <- c(rep(1, 20), rep(2, 20))
knn.ex <- knn3Train(train, decision.points, cl, k = 5, prob = TRUE)

knn.ex

reds <- rep(NA, 2500)

r2b <- colorRampPalette(c("black", "white", "red"))

for(i in 1:2500) reds[i] <- r2b(50)[which.min(abs(attr(knn.ex, "prob")[i,1] - seq(0,1, length.out = 50)))]

pdf("~/Downloads/knn.pdf", height = 4, width = 7)
plot(decision.points[,2] ~ decision.points[,1], col = reds, pch = 19, cex = .3, xlab = "", ylab = "")

dev.off()

pdf("~/Downloads/knn_discrete.pdf", height = 4, width = 7)
plot(decision.points[,2] ~ decision.points[,1], col = knn.ex, pch = 19, cex = .3, xlab = "", ylab = "")

dev.off()

#plsr####
pls.df <- data.frame(cl = as.factor(cl), train)
model.pls <- train(cl ~ ., data = pls.df, method = 'pls')

model.pls

colnames(decision.points) <- colnames(pls.df[-1])

pls.preds <- predict(model.pls, newdata =  decision.points, type = "prob")

r2b <- colorRampPalette(c("black", "white", "red"))

for(i in 1:2500) reds[i] <- r2b(50)[which.min(abs(pls.preds[i,1] - seq(0,1, length.out = 50)))]

pdf("~/Downloads/pls.pdf", height = 4, width = 7)
plot(decision.points[,2] ~ decision.points[,1], col = reds, pch = 19, cex = .3, xlab = "", ylab = "")
dev.off()

#rf####
rf.df <- data.frame(cl = as.factor(cl), train)
model.rf <- train(cl ~ ., data = rf.df, method = 'rf')

model.rf

colnames(decision.points) <- colnames(rf.df[-1])

rf.preds <- predict(model.rf, newdata =  decision.points, type = "prob")

r2b <- colorRampPalette(c("black", "white", "red"))

for(i in 1:2500) reds[i] <- r2b(50)[which.min(abs(rf.preds[i,1] - seq(0,1, length.out = 50)))]

pdf("~/Downloads/rf.pdf", height = 4, width = 7)
plot(decision.points[,2] ~ decision.points[,1], col = reds, pch = 19, cex = .3, xlab = "", ylab = "")
dev.off()

#nn####
nn.df <- data.frame(cl = as.factor(cl), train)
model.nn <- train(cl ~ ., data = nn.df, method = 'mlpML', tuneGrid = expand.grid(layer1 = 1:3, layer2 = 0:2, layer3 = 0:2))

model.nn

colnames(decision.points) <- colnames(nn.df[-1])

nn.preds <- predict(model.nn, newdata =  decision.points, type = "prob")

r2b <- colorRampPalette(c("black", "white", "red"))

for(i in 1:2500) reds[i] <- r2b(50)[which.min(abs(nn.preds[i,1] - seq(0,1, length.out = 50)))]

pdf("~/Downloads/nn.pdf", height = 4, width = 7)
plot(decision.points[,2] ~ decision.points[,1], col = reds, pch = 19, cex = .3, xlab = "", ylab = "")
dev.off()








