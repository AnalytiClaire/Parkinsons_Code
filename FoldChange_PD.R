#############################################################
################### FOLD CHANGE DEG #########################
#############################################################

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

LEW <- read.csv("LEWfilteredresult.csv")
LEW <- LEW[order(LEW$P.Value),]
MID3 <- read.csv("MID3filteredresult.csv")
MID3 <- MID3[order(MID3$P.Value),]
MID4 <- read.csv("MID4filteredresult.csv")
MID4 <- MID4[order(MID4$P.Value),]
MOR.FC <- read.csv("MOR.FCfilteredresult.csv")
MOR.FC <- MOR.FC[order(MOR.FC$P.Value),]
DIJ <- read.csv("DIJfilteredresult.csv")
DIJ <- DIJ[order(DIJ$P.Value),]
FFR <- read.csv("FFRfilteredresult.csv")
FFR <- FFR[order(FFR$P.Value),]
MID1 <- read.csv("MID1filteredresult.csv")
MID1 <- MID1[order(MID1$P.Value),]
MID2 <- read.csv("MID2filteredresult.csv")
MID2 <- MID2[order(MID2$P.Value),]
MOR.SN <- read.csv("MOR.SNfilteredresult.csv")
MOR.SN <- MOR.SN[order(MOR.SN$P.Value),]


thresh <- 1

upLEW <- subset(LEW, LEW$Fold.Change >= thresh)
upLEWgene <- upLEW$Gene.Symbol

upMID3<- subset(MID3, MID3$Fold.Change >= thresh)
upMID3gene <- upMID3$Gene.Symbol

upMID4 <- subset(MID4, MID4$Fold.Change >= thresh)
upMID4gene <- upMID4$Gene.Symbol

upMOR.FC <- subset(MOR.FC, MOR.FC$Fold.Change >= thresh)
upMOR.FCgene <- upMOR.FC$Gene.Symbol

# upDUM <- subset(DUM, DUM$FoldChange >= thresh)
# upDUMgene <- upDUM$hgnc_symbol

upDIJ <- subset(DIJ, DIJ$Fold.Change >= thresh)
upDIJgene <- upDIJ$Gene.Symbol

upFFR<- subset(FFR, FFR$Fold.Change >= thresh)
upFFRgene <- upFFR$Gene.Symbol

upMID1<- subset(MID1, MID1$Fold.Change >= thresh)
upMID1gene <- upMID1$Gene.Symbol

upMID2 <- subset(MID2, MID2$Fold.Change >= thresh)
upMID2gene <- upMID2$Gene.Symbol

upMOR.SN <- subset(MOR.SN, MOR.SN$Fold.Change >= thresh)
upMOR.SNgene <- upMOR.SN$Gene.Symbol


INTUP <- Reduce(intersect, list(upLEWgene, upMID3gene, upMID4gene, upMOR.FCgene,
                                upDIJgene, upFFRgene, upMID1gene, upMID2gene, upMOR.SNgene))



#### DOWN ####
thresh <- -1

downLEW <- subset(LEW, LEW$Fold.Change <= thresh)
downLEWgene <- downLEW$Gene.Symbol

downMID3 <- subset(MID3, MID3$Fold.Change <= thresh)
downMID3gene <- downMID3$Gene.Symbol

downMID4 <- subset(MID4, MID4$Fold.Change <= thresh)
downMID4gene <- downMID4$Gene.Symbol

downMOR.FC <- subset(MOR.FC, MOR.FC$Fold.Change <= thresh)
downMOR.FCgene <- downMOR.FC$Gene.Symbol

# downDUM <- subset(DUM, DUM$FoldChange <= thresh)
# downDUMgene <- downDUM$hgnc_symbol

downDIJ <- subset(DIJ, DIJ$Fold.Change <= thresh)
downDIJgene <- downDIJ$Gene.Symbol

downFFR<- subset(FFR, FFR$Fold.Change <= thresh)
downFFRgene <- downFFR$Gene.Symbol

downMID1<- subset(MID1, MID1$Fold.Change <= thresh)
downMID1gene <- downMID1$Gene.Symbol

downMID2 <- subset(MID2, MID2$Fold.Change <= thresh)
downMID2gene <- downMID2$Gene.Symbol

downMOR.SN <- subset(MOR.SN, MOR.SN$Fold.Change <= thresh)
downMOR.SNgene <- downMOR.SN$Gene.Symbol

INTDOWN <- Reduce(intersect, list(downLEWgene, downMID3gene, downMID4gene, downMOR.FCgene,
                                  downDIJgene, downFFRgene, downMID1gene, downMID2gene, downMOR.SNgene))




########################### COMMON GENES ##############################
all <- c(INTUP, INTDOWN)
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/FoldChange")
write.table(INTUP,"upDEGs.txt", row.names = F, col.names = F, quote = F) 
write.table(INTDOWN,"downDEGs.txt", row.names = F, col.names = F, quote = F) 
write.table(all, "allDEGs.txt", row.names = F, col.names = F, quote = F)

PDgenes <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/ParkinsonsDiseaseMalacards.txt")

intersect(INTUP, PDgenes)
intersect(INTDOWN, PDgenes)
