##Differential Expression of Genes##

library(affy)
library(tcltk)
library(widgetTools)
library(DynDoc)
library(tools)
library(Biobase)
library(tkWidgets)
library(plyr)
###########################
##### ALPHA SYNUCLEIN #####
###########################

##### LEW GSE19587 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE19587_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"LEW" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(8,14)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)

##### MOR.FC GSE8397 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE8397_Cortex/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MOR.FC" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(3,5)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)



##### MID3 GSE20168 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20168_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MID3" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(15,14)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)





##### MID4 GSE20291 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20291_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MID4" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(20,15)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/AlphaSynuclein/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)





###########################
##### MITOCHONDRIA ########
###########################

##### MID1 GSE20141 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20141_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"DAV" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(8,10)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)


##### FFR GSE7621 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE7621_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"FFR" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A_2.na35.annot.csv/HG-U133A_2.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(9,16)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), sep="\t", row.names=FALSE, quote = FALSE)
##### MID2 GSE20292 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE20292_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MID" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(18,11)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)


##### DIJ GSE49036 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE49036_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"DIJ" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133_Plus_2.na35.annot.csv/HG-U133_Plus_2.na35_SHORT.annot.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat<-factor(rep(c("Control", "Patient"),c(8,15)), levels=c("Control", "Patient"))
design<-model.matrix(~Treat)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)


##### MOR.SN GSE8397 #####
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Data/GSE8397_RAW/")

#run program to choose .CEL files from directory
celfiles <- fileBrowser(textToShow = "Choose CEL files", testFun = hasSuffix("[cC][eE][lL]"))
Data<-ReadAffy(filenames=celfiles) #read in files
rmaEset<-rma(Data) #normalise using RMA
analysis.name<-"MOR" #Label analysis
dataMatrixAll<-exprs(rmaEset) #takes expression from normalised expression set


#mas5call generates presence/absence calls for each probeset
mas5call<-mas5calls(Data)
callMatrixAll<-exprs(mas5call)
colnames(callMatrixAll)<-sub(".CEL", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
colnames(callMatrixAll)<-sub(".cel", ".mas5-Detection", colnames(callMatrixAll),fixed=TRUE)
callMatrixAll<-as.data.frame(callMatrixAll)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPf<-function(x){
  sum(x=="P")
}

#count how many samples have presence calls
countPl<-apply(callMatrixAll, 1, countPf)
callMatrixAll$ProbeSetID<-rownames(callMatrixAll)
countPdf<-data.frame(ProbeSetID=names(countPl), countP=countPl) 

#read annotation
# USING ANNOTATION FILE (if .csv, convert to .txt using excel)
annotation.file<-"/Users/clairegreen/Documents/PhD/TDP-43/TDP-43_Data/HG-U133A.na35.annot.csv/HG-U133A.na35.annot_SHORT.txt"
annotation<-read.table(annotation.file, header = TRUE, row.names=NULL, sep="\t", skip=0, stringsAsFactors=F, quote = "", comment.char="!", fill = TRUE, as.is = TRUE)
dim(annotation)
nrow(annotation)
annotation<-subset( annotation, subset=(Gene.Symbol !="---")) #if no gene symbol, discount

# Remove rows in which genes are noted to have negative strand matching probes
idxNegativeStrand<-grep("Negative Strand Matching Probes", annotation$Annotation.Notes)
if(length(idxNegativeStrand)>0)
{
  annotation<-annotation[-idxNegativeStrand,]
}


expressionMatrix<-exprs(rmaEset)
colnames(expressionMatrix)

#this is for matched samples
Treat <- factor(rep(c("Control", "Patient"),c(15,24)), levels=c("Control", "Patient"))
Tissue <- factor(c(rep(c("SN", "LN"),c(7,9)),rep(c("SN", "LN"),c(8,15))), levels=c("SN", "LN"))

design<-model.matrix(~Treat+Tissue)
rownames(design)<-colnames(expressionMatrix)
design

#Conduct statistical analysis of expression
library(limma)
fit<-lmFit(expressionMatrix, design) #linear model fit
fit<-eBayes(fit) 
result<-topTable(fit, coef="TreatPatient", adjust="BH", number=nrow(expressionMatrix)) #"BH" adjust for multiple hypothesis testing
#toptable normally takes top number but this takes all
setwd("/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression")
write.csv(result, file = paste(analysis.name, "_result.csv", sep = ""))

result$"ProbeSetID"<-rownames(result) #make probeset IDs the row names
head(result$"ProbeSetID") 
result$"Fold Change"<-2^result$logFC 
result$"Fold Change"[result$"Fold Change"<1]<-(-1)/result$"Fold Change"[result$"Fold Change"<1] #converts log fold change into a linear value above or below 0
expressionLinear<-as.data.frame(2^expressionMatrix)
expressionLinear$ProbeSetID<-rownames(expressionLinear)
result<-merge(result, expressionLinear, by.x="ProbeSetID", by.y="ProbeSetID") #merge values into one array
result<-merge(annotation, result, by.x="Probe.Set.ID", by.y="ProbeSetID")
result<-merge(result, countPdf, by.x="Probe.Set.ID", by.y="ProbeSetID")
result$Gene.Symbol <- sapply(strsplit(result$Gene.Symbol,"///"), `[`, 1)
result$Ensembl <- sapply(strsplit(result$Ensembl,"///"), `[`, 1)

setwd(dir = "/Users/clairegreen/Documents/PhD/Parkinsons/Parkinsons_Code/Results/GeneExpression/")
write.csv(result, file=paste(analysis.name, "_finalresult.csv", sep=""), row.names=FALSE, quote = FALSE)

genesort <- result[order(result$P.Value),]
uniqueresult <- genesort[!duplicated(genesort[,5]),]
write.csv(uniqueresult, file=paste(analysis.name, "rankeduniqueresult.csv", sep=""), row.names=FALSE, quote = FALSE)

result<-subset(result, Gene.Symbol!="") #removes any probes for which there are no gene symbols
result<-subset(result, subset=(countP>2)) #only takes results that have at least 2 samples with a presence call for a probe
write.csv(uniqueresult, file=paste(analysis.name, "filteredresult.csv", sep=""), row.names=FALSE, quote = FALSE)

