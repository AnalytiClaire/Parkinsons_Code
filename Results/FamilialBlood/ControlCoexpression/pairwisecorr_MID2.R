#####MID2####
uniqueresult <- read.csv("MID2_DEGfilt.csv", row.names = 1)

##For loop for generating regression values and p values
CorExprMat <- t(uniqueresult)

reg <- matrix(0, ncol(CorExprMat), ncol(CorExprMat))
p.value <- matrix(0, ncol(CorExprMat), ncol(CorExprMat))

for (i in 1:ncol(CorExprMat)){
  for (j in 1:ncol(CorExprMat)){
    reg[i,j] <- cor.test(CorExprMat[,i], CorExprMat[,j], method = "spearman")$estimate
  }}

rownames(reg) <- colnames(reg) <- colnames(CorExprMat)

for (i in 1:ncol(CorExprMat)){
  for (j in 1:ncol(CorExprMat)){
    p.value[i,j] <- cor.test(CorExprMat[,i], CorExprMat[,j], method = "pearson")$p.value
  }}

rownames(p.value) <- colnames(p.value) <- colnames(CorExprMat)

##Only take upper triangle without diagonal (all comparisons are currently doubled)
ptri <- p.value
ptri[lower.tri(ptri, diag = TRUE)] <- NA

#Turn into vector
library(gdata)
p.vec <- unmatrix(ptri)
#Remove NA values
p.vec <- na.omit(p.vec)
#Multiple hypothesis testing correction 
p.adj <- p.adjust(p.vec, method = "fdr", n = length(p.vec))

#Create results table
reg.mat <- unmatrix(reg)
reg.mat <- as.data.frame(reg.mat)
p.adj <- as.data.frame(p.adj)
p.mat <- as.data.frame(p.vec)

pvals <- merge(p.adj, p.mat, by.x = "row.names", by.y = "row.names")
rownames(pvals)<- pvals$Row.names
pvals[,1] <- NULL
results <- merge(pvals, reg.mat, by.x = "row.names", by.y = "row.names")
rownames(results)<- results$Row.names
results[,1] <- NULL
results <- results[order(results$p.vec),]

write.csv(results, "MID2corresult.csv")