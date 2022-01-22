#register non-synds with to atlas one at a time

load("/Users/jovid/Downloads/Segmented_PLS/30k.Rdata")
d.meta.ns <- d.meta

sample1k <- sample(1:27903, 1000)
rot.ns <- array(NA, dim = c(27903, 3, 2273))
tmp.mesh <- atlas
for(i in 1:dim(rot.ns)[3]){

tmp.mesh$vb[-4,] <- t(d.registered$rotated[,,i])  
  
trafo1k <- rotmesh.onto(tmp.mesh, d.registered$rotated[sample1k,,i], synd.mshape[sample1k,])
rot.ns[,,i] <- t(trafo1k$mesh$vb[-4,])
print(i)
}

ns.PC.scores <- getPCscores(rot.ns, PC.eigenvectors, synd.mshape)

#put data together
load("data.Rdata")

View(d.meta)

PC.scores[,1:2]

rownames(ns.PC.scores) <- d.meta.ns$FBID

PC.scores <- rbind(PC.scores, ns.PC.scores)

colnames(d.meta.ns)

colnames(d.meta)

ns.meta <- cbind(d.meta.ns[,c(1,3,6)], rep("Non-syndromic", 2273), d.meta.ns[,2])
colnames(ns.meta)[4:5] <- c("Syndrome", "scan") 

#match sex cov
# tmp.mesh$vb[-4,] <- t(d.registered$rotated[,,6])  
# plot3d(vcgSmooth(tmp.mesh), aspect = "iso", col = "lightgrey", specular = 1)
ns.meta$Sex <- ns.meta$Sex - 1 

d.meta <- rbind(d.meta, ns.meta)
