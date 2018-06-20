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
DUM <- read.csv("DUM_UniqueGene_DESeq2.csv")
DUM <- DUM[order(DUM$pvalue),]
BOT <- read.csv("BOTrankeduniqueresult.csv")
BOT <- BOT[order(BOT$P.Value),]
BOT2 <- read.csv("BOT2rankeduniqueresult.csv")
BOT2 <- BOT2[order(BOT2$P.Value),]




thresh <- 1

upLEW <- subset(LEW, LEW$Fold.Change >= thresh)
upLEWgene <- upLEW$Gene.Symbol

upMID3<- subset(MID3, MID3$Fold.Change >= thresh)
upMID3gene <- upMID3$Gene.Symbol

upMID4 <- subset(MID4, MID4$Fold.Change >= thresh)
upMID4gene <- upMID4$Gene.Symbol

upMOR.FC <- subset(MOR.FC, MOR.FC$Fold.Change >= thresh)
upMOR.FCgene <- upMOR.FC$Gene.Symbol

upDUM <- subset(DUM, DUM$log2FoldChange >= 0)
upDUMgene <- upDUM$hgnc_symbol

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

upBOT <- subset(BOT, BOT$Fold.Change >= thresh)
upBOTgene <- upBOT$Gene.Symbol

upBOT2 <- subset(BOT2, BOT2$Fold.Change >= thresh)
upBOT2gene <- upBOT2$Gene.Symbol


INTUP <- Reduce(intersect, list(upLEWgene, upMID3gene, upMID4gene, upMOR.FCgene,
                                upDIJgene, upFFRgene, upMID1gene, upMID2gene, 
                                upMOR.SNgene, upDUMgene, upBOTgene, upBOT2gene))



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

downDUM <- subset(DUM, DUM$log2FoldChange <= 0)
downDUMgene <- downDUM$hgnc_symbol

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


downBOT <- subset(BOT, BOT$Fold.Change <= thresh)
downBOTgene <- downBOT$Gene.Symbol

downBOT2 <- subset(BOT2, BOT2$Fold.Change <= thresh)
downBOT2gene <- downBOT2$Gene.Symbol

INTDOWN <- Reduce(intersect, list(downLEWgene, downMID3gene, downMID4gene, downMOR.FCgene,
                                  downDIJgene, downFFRgene, downMID1gene, downMID2gene, 
                                  downMOR.SNgene, downDUMgene, downBOTgene, downBOT2gene))




########################### COMMON GENES ##############################
all <- c(INTUP, INTDOWN)
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/FoldChange")
# write.table(INTUP,"Sp_LRRK2_upDEGs.txt", row.names = F, col.names = F, quote = F)
# write.table(INTDOWN,"Sp_LRRK2_downDEGs.txt", row.names = F, col.names = F, quote = F)
# write.table(all, "Sp_LRRK2_allDEGs.txt", row.names = F, col.names = F, quote = F)

PDgenes <- readLines("/Users/clairegreen/Documents/PhD/Parkinsons/ParkinsonsDiseaseMalacards.txt")

intersect(INTUP, PDgenes)
intersect(INTDOWN, PDgenes)



####### ALS signature ##########

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/noMedian/")

C9 <- read.csv("C9_unique.csv")
C9 <- C9[order(C9$P.Value),]
sals <- read.csv("sals_unique.csv")
sals <- sals[order(sals$P.Value),]

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/TDP-43_DEseq2/")

pet <- read.csv("PET_results_keepfiltering.csv")
rav <- read.csv("RAV_results_keepfiltering.csv")

setwd("/users/clairegreen/Documents/PhD/TDP-43/TDP-43_Code/Results/GeneExpression/non-TDP/")

FUS <- read.csv("FUSrankeduniqueresult.csv")
FUS <- FUS[order(FUS$P.Value),]
SOD1 <- read.csv("SOD1rankeduniqueresult.csv")
SOD1 <- SOD1[order(SOD1$P.Value),]


thresh <- 1

upC9 <- subset(C9, C9$Fold.Change >= thresh)
upC9gene <- upC9$Gene.Symbol

upSALS <- subset(sals, sals$Fold.Change >= thresh)
upSALSgene <- upSALS$Gene.Symbol

upPET <- subset(pet, pet$FoldChange >= thresh)
upPETgene <- upPET$hgnc_symbol

upRAV <- subset(rav, rav$FoldChange >= thresh)
upRAVgene <- upRAV$hgnc_symbol

upFUS <- subset(FUS, FUS$Fold.Change >= thresh)
upFUSgene <- upFUS$Gene.Symbol

upSOD1 <- subset(SOD1, SOD1$Fold.Change >= thresh)
upSOD1gene <- upSOD1$Gene.Symbol

INTUP_ALS <- Reduce(intersect, list(upC9gene, upSALSgene, upPETgene, upRAVgene, upFUSgene, upSOD1gene))


#### DOWN 
thresh <- -1

downC9 <- subset(C9, C9$Fold.Change <= thresh)
downC9gene <- downC9$Gene.Symbol

downSALS <- subset(sals, sals$Fold.Change <= thresh)
downSALSgene <- downSALS$Gene.Symbol

downPET <- subset(pet, pet$FoldChange <= thresh)
downPETgene <- downPET$hgnc_symbol

downRAV <- subset(rav, rav$FoldChange <= thresh)
downRAVgene <- downRAV$hgnc_symbol

downFUS <- subset(FUS, FUS$Fold.Change <= thresh)
downFUSgene <- downFUS$Gene.Symbol

downSOD1 <- subset(SOD1, SOD1$Fold.Change <= thresh)
downSOD1gene <- downSOD1$Gene.Symbol

INTDOWN_ALS <- Reduce(intersect, list(downC9gene, downSALSgene, downPETgene, downRAVgene, downFUSgene, downSOD1gene))



##### COMMON GENES ###
upremove <- Reduce(intersect, list (INTUP, INTUP_ALS))
downremove <- Reduce(intersect, list(INTDOWN, INTDOWN_ALS))



###### REMOVE COMMON GENES ###
resultsup <- subset(INTUP, !(INTUP %in% upremove))
resultsdown <- subset(INTDOWN, !(INTDOWN %in% downremove))
results <- c(resultsup, resultsdown)

####### sPD Blood signature ##########

setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")

AMA <- read.csv("AMAfilteredresult.csv")
AMA <- AMA[order(AMA$P.Value),]
RON <- read.csv("RONfilteredresult.csv")
RON <- RON[order(RON$P.Value),]


thresh <- 1

upAMA <- subset(AMA, AMA$Fold.Change >= thresh)
upAMAgene <- upAMA$Gene.Symbol

upRON <- subset(RON, RON$Fold.Change >= thresh)
upRONgene <- upRON$Gene.Symbol


INTUP_blood <- Reduce(intersect, list(upAMAgene, upRONgene))


#### DOWN ###
thresh <- -1

downAMA <- subset(AMA, AMA$Fold.Change <= thresh)
downAMAgene <- downAMA$Gene.Symbol

downRON <- subset(RON, RON$Fold.Change <= thresh)
downRONgene <- downRON$Gene.Symbol


INTDOWN_blood <- Reduce(intersect, list(downAMAgene, downRONgene))



##### COMMON GENES ###
upremove2 <- Reduce(intersect, list(INTUP, INTUP_blood))
downremove2 <- Reduce(intersect, list(INTDOWN, INTDOWN_blood))






##### REMOVE COMMON GENES ###
resultsup <- subset(resultsup, !(resultsup %in% upremove3))
resultsdown <- subset(resultsdown, !(resultsdown %in% downremove3))
results <- c(resultsup, resultsdown)




setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/FamilialBlood/")
write.table(resultsup, "ALS_sfblood_UPgenes.txt", quote = F, row.names = F, col.names = F)
write.table(resultsdown, "ALS_sfblood_DOWNgenes.txt", quote = F, row.names = F, col.names = F)
write.table(results, "ALS_sfblood_ALLgenes.txt", quote = F, row.names = F, col.names = F)
cat(resultsup, sep="\n")

intersect(resultsup, PDgenes)
intersect(resultsdown, PDgenes)





# ####### LRRK2/Parkin Fibroblast Signature ##########
# 
# setwd("/users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
# 
# LRRK2F <- read.csv("LRRK2Fibrankeduniqueresult.csv")
# LRRK2F <- LRRK2F[order(LRRK2F$P.Value),]
# PARK2F <- read.csv("PARK2Fibrankeduniqueresult.csv")
# PARK2F <- PARK2F[order(PARK2F$P.Value),]
# 
# 
# thresh <- 1
# 
# upLRRK2F <- subset(LRRK2F, LRRK2F$FoldChange >= thresh)
# upLRRK2Fgene <- as.character(upLRRK2F$Gene.Symbol)
# 
# upPARK2F <- subset(PARK2F, PARK2F$FoldChange >= thresh)
# upPARK2Fgene <- as.character(upPARK2F$Gene.Symbol)
# 
# 
# INTUP_LRRK2F <- Reduce(intersect, list(resultsup, upLRRK2Fgene))
# INTUP_PARK2F <- Reduce(intersect, list(resultsup, upPARK2Fgene))
# 
# 
# #### DOWN ####
# thresh <- -1
# 
# downLRRK2F <- subset(LRRK2F, LRRK2F$FoldChange <= thresh)
# downLRRK2Fgene <- downLRRK2F$Gene.Symbol
# 
# downPARK2F <- subset(PARK2F, PARK2F$FoldChange <= thresh)
# downPARK2Fgene <- downPARK2F$Gene.Symbol
# 
# 
# INTDOWN_LRRK2F <- Reduce(intersect, list(resultsdown, downLRRK2Fgene))
# INTDOWN_PARK2F <- Reduce(intersect, list(resultsdown, downPARK2Fgene))
# 
# 
# 
# ################ REMOVE COMMON GENES ######################
# resultsup <- subset(resultsup, !(resultsup %in% INTUP_LRRK2F))
# resultsup <- subset(resultsup, !(resultsup %in% INTUP_PARK2F))
# 
# resultsdown <- subset(resultsdown, !(resultsdown %in% INTDOWN_LRRK2F))
# resultsdown <- subset(resultsdown, !(resultsdown %in% INTDOWN_PARK2F))
# 
# results <- c(resultsup, resultsdown)
# 



# setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/ALSPARK2LRRK2/")
# write.table(resultsup, "removeallUPgenes.txt", quote = F, row.names = F, col.names = F)
# write.table(resultsdown, "removeallDOWNgenes.txt", quote = F, row.names = F, col.names = F)
# write.table(results, "removeallALLgenes.txt", quote = F, row.names = F, col.names = F)
# cat(resultsup, sep="\n")

intersect(resultsup, PDgenes)
intersect(resultsdown, PDgenes)
