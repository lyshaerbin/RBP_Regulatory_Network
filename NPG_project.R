setwd("C:/projects/RBP/newidea")
RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
CGC<-read.csv("Census_allSun Nov 13 17-37-13 2016.tsv",stringsAsFactors=F,sep="\t",skip=0,header =T)
library(vioplot)
setwd("C:/projects/TCGAmutation")
Allfreq<-read.csv("LIHC_freq.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
RBP_frq<-merge(RBP,Allfreq,by.x = "HGNC.symbol",by.y = "V1")
CGC_frq<-merge(CGC,Allfreq,by.x = "Gene.Symbol",by.y = "V1")
RC_gene<-union(RBP$HGNC.symbol,CGC$Gene.Symbol)
other_g<-setdiff(Allfreq$V1,RC_gene)
other_g<-t(other_g)
other_g<-t(other_g)
other_frq<-merge(other_g,Allfreq,by.x="V1",by.y="V1")
##############random selected genes
othergene=setdiff(Allfreq$V1,RBP$HGNC.symbol)
othergene=t(t(othergene))
other_frq<-merge(othergene,Allfreq,by.x="V1",by.y="V1")
X2=other_frq$V2
Randfreq<-c()
for (i in 1:100) {
  SS=sample(15137,1350)
  Randfreq=rbind(Randfreq,X2[SS])
}
boxplot(as.numeric(RBP_frq$V2),as.numeric(Randfreq),ylab="Frequency",ylim=c(0,0.05))
wilcox.test(RBP_frq$V2,Randfreq)

#######RBP-sample mutation matrix
setwd("C:/projects/RBP/newidea")
RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
setwd("C:/projects/TCGAmutation")
mutation<-read.csv("LIHC.maf",stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
LIHC=unique(mutation$Tumor_Sample_Barcode)
Mutmat<-c()
for (i in 1:dim(RBP)[1]) {
  RBP_Sa<-c()
  for (j in 1:length(LIHC)) {
    x<-which(mutation$Hugo_Symbol==RBP$HGNC.symbol[i]&mutation$Tumor_Sample_Barcode==LIHC[j])
    if(length(x)>0){
      RBP_Sa<-rbind(RBP_Sa,1)
    } else{
      RBP_Sa<-rbind(RBP_Sa,0)
    }
  }
  Mutmat<-cbind(Mutmat,RBP_Sa)
}
KK=apply(Mutmat,2,sum)
RBP_num=c()
for (i in 1:10) {
  x=which(KK==i)
  RBP_num<-rbind(RBP_num,length(x))
}
x=which(KK>10)
RBP_num<-rbind(RBP_num,length(x))
par(mfrow=c(3,2))
barplot(t(RBP_num),ylab = "Number of samples")
library(pheatmap)
pheatmap(Mutmat)
WW=apply(Mutmat, 1, sum)
x=which(KK>0)
y=which(WW>0)
Mutmat_1=Mutmat[y,x]
hist(WW)

############################download the LIHC gene expression datasets
library(TCGAbiolinks)
# Downloading and prepare 
query <- GDCquery(project = "TCGA-LIHC", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
LIHCdata <- GDCprepare(query, save = TRUE, 
                       save.filename = "LIHCExpression.rda",
                       remove.files.prepared = TRUE)
library(SummarizedExperiment)
LIHCmatrix<-assay(LIHCdata,1,"FPKM")############gene exprssion
###################filter the gene expression profile
filter_line<-c()
for (i in 1:dim(LIHCmatrix)[1]) {
  n=which(LIHCmatrix[i,]<0.1)
  if(length(n)>424*0.3){
    filter_line<-rbind(filter_line,i)
  }
}
LIHCmatrix<-LIHCmatrix[-filter_line,]
LIHCclinal<-LIHCdata@colData@listData$definition#############sample information tumor vs normal
LIHCclinal=t(t(LIHCclinal))
RBP<-read.csv("EnsemRBP.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
LIHCmatrix_nor=log2(LIHCmatrix+0.05)
M_exp=c()
V_exp=c()
for (i in 1:dim(LIHCmatrix_nor)[1]) {
  print(i)
  M_exp=rbind(M_exp,mean(LIHCmatrix_nor[i,]))
  V_exp=rbind(V_exp,var(LIHCmatrix_nor[i,]))
}
GeneMV=cbind(rownames(LIHCmatrix),M_exp,V_exp)
RBPMV=merge(RBP,GeneMV,by.x="V1",by.y="V1")
othergene=setdiff(GeneMV[,1],RBPMV[,1])
othergene=t(t(othergene))
OtherMV=merge(othergene,GeneMV,by.x="V1",by.y="V1")
X1=as.numeric(as.character(RBPMV[,6]))
X2=as.numeric(as.character(OtherMV$V2))
X3=as.numeric(as.character(RBPMV[,7]))
X4=as.numeric(as.character(OtherMV$V3))
library(vioplot)
par(mfrow=c(3,3))
vioplot(X1,X2)
vioplot(X3,X4)
wilcox.test(X1,X2)
wilcox.test(X3,X4)#########p<2.2e-16 RBPs have higher expression and lower varation
write.table(RBPMV,file = "RBP_MV.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)

##############Identify the differentially expressed RBPs 
FC<-c()
Wp<-c()
Normal_s=which(LIHCclinal=="Solid Tissue Normal")
Cancer_s=which(LIHCclinal=="Primary solid Tumor"|LIHCclinal=="Recurrent Solid Tumor")
for (i in 1:dim(LIHCmatrix_nor)[1]) {
  print(i)
  FC=rbind(FC,mean(LIHCmatrix[i,Cancer_s])/mean(LIHCmatrix[i,Normal_s]))
  p1=wilcox.test(LIHCmatrix_nor[i,Cancer_s],LIHCmatrix_nor[i,Normal_s])$p.value
  Wp<-rbind(Wp,p1)
}
FDR=p.adjust(Wp,method = "fdr")
Diffexp<-cbind(rownames(LIHCmatrix),FC,Wp,FDR)
Up=which(Diffexp[,2]>2&Diffexp[,4]<0.01)
Down=which(Diffexp[,2]<0.5&Diffexp[,4]<0.01)
Sigdiff<-rbind(Diffexp[Up,],Diffexp[Down,])

#############rbp differentially expressed plot
RBPdexp<-merge(RBP,Diffexp,by.x="V1",by.y="V1")
write.table(RBPdexp,file = "RBP_Differnt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

lfc <- 1
pval <- 0.01
y <- -log10(as.numeric(as.character(RBPdexp$FDR)))
x <- log2(as.numeric(as.character(RBPdexp$V2.y)))
c <- -log10(as.numeric(as.character(RBPdexp$V3.y)))
bb <- data.frame(x,y,c)
upGenes <- which((x>1) & (y>2) )
downGenes <- which((x<(-1)) & (y>2))

pdf("volcano.pdf",width = 10,height = 6)
par(mar = c(5,4,4,4))
plot(bb$x,bb$y,pch=16, col="lightgray", cex=0.6, xlab = "log2fold-change",ylab="-log10(FDR)")
points(bb[upGenes,]$x,bb[upGenes,]$y,pch=16,cex=0.8,col="firebrick")
points(bb[downGenes,]$x,bb[downGenes,]$y,pch=16,cex=0.8,col="dodgerblue3")
abline(h=-log10(pval),col="lightgray",lty=2)
abline(v=c(-1,1),col="lightgray",lty=2)
mtext(paste("pval = ", pval),side=4, at=-log10(pval),cex=0.7,line=0.5,las=1)
mtext(c(paste("-",lfc,"fold"),paste("+",lfc,"fold")),side=3,at=c(-lfc,lfc),cex = 0.8,line=0.5)
dev.off()

##############differentially expressed heatmap
library(gplots)
dexin=union(upGenes,downGenes)
DexRBP<-RBPdexp$V1[dexin]
######DF RBP expression profile
dfRBPexp=c()
for (i in 1:length(DexRBP)) {
  x=which(rownames(LIHCmatrix_nor)==DexRBP[i])
  dfRBPexp<-rbind(dfRBPexp,LIHCmatrix_nor[x,])
}
color.map <- function(samp) { if (samp=="Solid Tissue Normal") "springgreen4" else "orchid4" }
patientcolors <- sapply(LIHCclinal, color.map)
heatmap.2(dfRBPexp, labRow = DexRBP, labCol=NULL,Rowv = FALSE, col=bluered(100), scale="row",ColSideColors=patientcolors,symkey=FALSE, density.info="density", trace="none", cexRow=0.5)

library(pheatmap)
pheatmap(dfRBPexp,scale = "row",cluster_rows=FALSE,color = patientcolors)


######################RBP network analysis
setwd("C:/projects/RBP/newidea")
RBPnet<-read.csv("RBP_gene_inter_HepG2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBP=unique(RBPnet$V1)
RBP_degree<-c()
for (i in 1:length(RBP)) {
  x=which(RBPnet$V1==RBP[i])
  RBP_degree<-rbind(RBP_degree,length(x))
}
RBP_degree<-cbind(RBP,RBP_degree)
RBPMV=read.csv("RBP_MV.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
RBP_dmv=merge(RBP_degree,RBPMV,by.x="RBP",by.y="V2.x")
X=as.numeric(as.character(RBP_dmv$V2))
Y=as.numeric(as.character(RBP_dmv$V2.y))
cor(log10(X),Y,method = "spearman")
plot(log10(X),Y)

######################RBP mutation analysis
setwd("C:/projects/RBP/newidea")
RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
setwd("C:/projects/TCGAmutation")
mutation<-read.csv("LIHC.maf",stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
setwd("C:/projects/RBP/newidea")
mutation_score<-read.csv("TCGA.LIHC.maf.hg38_multianno.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
RBP_mut<-merge(RBP,mutation,by.x = "HGNC.symbol",by.y = "Hugo_Symbol")
RBPmut_score<-merge(RBP_mut,mutation_score,by.x = c("Chromosome","Start_Position","End_Position","Tumor_Seq_Allele1","Tumor_Seq_Allele2"),
                    by.y = c("Chr","Start","End","Ref","Alt"))
Real_score=RBPmut_score$integrated_fitCons_score
Real_score=as.numeric(as.character(Real_score))
X=which(is.na(Real_score)==FALSE)
Real_score=Real_score[X]
n=length(Real_score)
N=dim(mutation_score)[1]
Rand_score=c()
for (i in 1:1000) {
  print(i)
  xx=sample(N,n)
  rr=mutation_score$integrated_fitCons_score[xx]
  rr=as.numeric(as.character(rr))
  y=which(is.na(rr)==FALSE)
  rr=rr[y]
  rr=t(t(rr))
  Rand_score=rbind(Rand_score,rr)
}
wilcox.test(Real_score,Rand_score,alternative = "greater")
library(vioplot)
par(mfrow=c(2,2))
vioplot(Real_score,Rand_score)

XXX=matrix(c(897,0,1350-897,38,151,1350-38-150,34,0,1316),nrow = 3)
XXX[,1]=XXX[,1]/sum(XXX[,1])
XXX[,2]=XXX[,2]/sum(XXX[,2])
XXX[,3]=XXX[,3]/sum(XXX[,3])
barplot(XXX, horiz =T)

#####################surface analysis
setwd("C:/projects/ORFchar")
subdir <- list.files(pattern="result",recursive=T)
write.table(subdir,file = "file.txt",sep = "\t",quote = F,row.names = F,col.names = F)
Surfscore=read.csv("Surf_score.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
TCGA=read.csv("TCGA-hg38-orfeome.tsv",stringsAsFactors = FALSE,sep = "\t",header = T)
xx=which(TCGA$cancer=="LIHC")
LIHC_orf=TCGA[xx,]
LIHC_surf=merge(LIHC_orf,Surfscore,by.x=c("ORF_ID","residue_."),by.y = c("V3","V4"))
setwd("C:/projects/TCGAmutation")
mutation<-read.csv("LIHC.maf",stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
mutation_surf<-merge(mutation,LIHC_surf,by.x = c("Chromosome","Start_Position","HGVSc"),by.y = c("chromosome","position","HGVSc"))
setwd("C:/projects/RBP/newidea")
RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
RBPmut_score<-merge(RBP,mutation_surf,by.x = "HGNC.symbol",
                    by.y = "Hugo_Symbol")
Real_score=RBPmut_score$V5
Real_score=as.numeric(as.character(Real_score))
X=which(is.na(Real_score)==FALSE)
Real_score=Real_score[X]
n=length(Real_score)

gene=unique(mutation_surf$Hugo_Symbol)
othergene=setdiff(gene,RBP$HGNC.symbol)
othergene=t(t(othergene))
otherscore=merge(mutation_surf,othergene,by.x = "Hugo_Symbol",by.y = "V1")
N=dim(otherscore)[1]
Rand_score=c()
KK=c()
for (i in 1:1000) {
  print(i)
  xx=sample(N,n)
  rr=otherscore$V5[xx]
  rr=as.numeric(as.character(rr))
  y=which(is.na(rr)==FALSE)
  rr=rr[y]
  rr=t(t(rr))
  Rand_score=rbind(Rand_score,rr)
  KK=rbind(KK,sum(rr>0.25))
}

wilcox.test(Real_score,Rand_score,alternative = "greater")
hist(KK,20, prob=TRUE,xlim = c(400,500)) 
lines(density(KK))
points(498,0)

library(vioplot)
par(mfrow=c(3,2))
vioplot(Real_score,Rand_score)

#####################analysis of RBP-RBP network
setwd("C:/projects/RBP/newidea")
RBPnet<-read.csv("Net.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBP=unique(RBPnet$V2)
RBP_RBP=c()
for (i in 1:dim(RBPnet)[1]) {
  print(i)
  x=which(RBP==RBPnet$V4[i])
  if(length(x)>0){
    RBP_RBP=rbind(RBP_RBP,RBPnet[i,])
  }
}

write.table(RBP_RBP,file = "RBP_RBP_net.txt",sep = "\t",quote = F,row.names = F,col.names = F)
RBP_RBP_id=c()
for (i in 1:dim(RBP_RBP)) {
  x1=which(RBP==RBP_RBP$V2[i])
  x2=which(RBP==RBP_RBP$V4[i])
  RBP_RBP_id=rbind(RBP_RBP_id,cbind(x1,x2))
}
write.table(RBP_RBP_id,file = "RBP_RBP_netid.txt",sep = "\t",quote = F,row.names = F,col.names = F)

RBP_matrix=c()
for (i in 1:65) {
  RBP_col=c()
  for (j in 1:65) {
    x=which(RBPnet$V2==RBP[i]&RBPnet$V4==RBP[j])
    if(length(x)>0){
      RBP_col=cbind(RBP_col,1)
    } else {
      RBP_col=cbind(RBP_col,0)
    }
  }
  RBP_matrix=rbind(RBP_matrix,RBP_col)
}
colnames(RBP_matrix)=RBP
rownames(RBP_matrix)=RBP
write.table(RBP_matrix,file = "RBP_RBP_MATRIX.txt",sep = "\t",quote = F,row.names = T,col.names = T)
####gene degree
gene=unique(RBPnet$V4)
genedegree=c()
for (i in 1:length(gene)) {
  x=which(RBPnet$V4==gene[i])
  RR=unique(RBPnet$V2[x])
  genedegree=rbind(genedegree,length(RR))
}
Mygenedegree=cbind(gene,genedegree)
par(mfrow=c(3,3))
colfunc <- colorRampPalette(c("white", "skyblue"))
XX=colfunc(33)
hist(genedegree,30,col =XX)
n1=length(which(genedegree>2))
n2=length(which(genedegree==1))
pie(c(n1,n2))
###############co-regulation vs interaction
HI3=read.csv("C:/projects/RBP/Pancancer/RBP_network/SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBP_coregulation=c()
for (i in 1:length(RBP)) {
  for (j in 1:length(RBP)) {
    x1=which(RBPnet$V2==RBP[i])
    x2=which(RBPnet$V2==RBP[j])
    target1=unique(RBPnet$V4[x1])
    target2=unique(RBPnet$V4[x2])
    nn=length(intersect(target1,target2))
    RBP_coregulation=rbind(RBP_coregulation,cbind(RBP[i],RBP[j],nn))
  }
}
Interornot=c()
for (i in 1:dim(RBP_coregulation)[1]) {
  x1=which(HI3$Symbol.for.A==RBP_coregulation[i,1]&HI3$Symbol.for.B==RBP_coregulation[i,2])
  x2=which(HI3$Symbol.for.A==RBP_coregulation[i,2]&HI3$Symbol.for.B==RBP_coregulation[i,1])
  if(length(x1)>0|length(x2)>0){
    Interornot=rbind(Interornot,1)
  } else{
    Interornot=rbind(Interornot,0)
  }
}
RBP_coregulation=cbind(RBP_coregulation,Interornot)
X1=as.numeric(as.character(RBP_coregulation[,3]))
X2=as.numeric(as.character(RBP_coregulation[,4]))
boxplot(X1~X2,boxwex=c(0.5,0.5),col=c("Gray","SkyBlue"))
x1=which(X2==0)
x2=which(X2==1)
wilcox.test(X1[x1],X1[x2])
RBP_inter_coreg=RBP_coregulation[x2,]

RBPnet1=RBPnet
xx=sample(dim(RBPnet1)[1],dim(RBPnet1)[1])
RBPnet$V4=RBPnet$V4[xx]
RBP_coregulation1=c()
for (i in 1:length(RBP)) {
  for (j in 1:length(RBP)) {
    x1=which(RBPnet$V2==RBP[i])
    x2=which(RBPnet$V2==RBP[j])
    target1=unique(RBPnet$V4[x1])
    target2=unique(RBPnet$V4[x2])
    nn=length(intersect(target1,target2))
    RBP_coregulation1=rbind(RBP_coregulation1,cbind(RBP[i],RBP[j],nn))
  }
}
Y1=as.numeric(as.character(RBP_coregulation1[,3]))
par(mfrow=c(3,3))
hist(Y1,20,probability = T,col = "gray")
lines(density(Y1),col="gray")
hist(X1,20,probability = T,col = "skyblue",add=T)
lines(density(X1),col="skyblue")
write.table(RBP_inter_coreg,file = "Intercoregu.txt",quote =FALSE,row.names = F,col.names = F,sep = "\t")
######################whether mutations are likely to locate in the binding sites of RBPs
setwd("C:/projects/RBP/newidea/random_regions")
subdir <- list.files(pattern="TCGA_rand",recursive=T)
Rand_N=c()
for (i in 1:100) {
  randm=read.csv(subdir[i],sep = "\t",header = FALSE,stringsAsFactors = FALSE)
  Rand_N=rbind(Rand_N,dim(randm)[1])
}
par(mfrow=c(3,3))
hist(Rand_N,30,probability = T,col = "gray")
lines(density(Rand_N),col="gray")
points(2164,0,col="red")
########################################
##whether the mutation located in the binding sites are more likely to be detelerous mutations
######################RBP mutation analysis
setwd("C:/projects/RBP/newidea")
mymut<-read.csv("mymut.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
setwd("C:/projects/TCGAmutation")
mutation<-read.csv("LIHC.maf",stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
setwd("C:/projects/RBP/newidea")
mutation_score<-read.csv("TCGA.LIHC.maf.hg38_multianno.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)

RBPmut_score<-merge(mymut,mutation_score,by.x = c("V6","V7","V8","V13","V14"),
                    by.y = c("Chr","Start","End","Ref","Alt"))
Real_score=RBPmut_score$MutationAssessor_score
Real_score=as.numeric(as.character(Real_score))
X=which(is.na(Real_score)==FALSE)
Real_score=Real_score[X]
n=length(Real_score)
N=dim(mutation_score)[1]
Rand_score=c()
for (i in 1:1000) {
  print(i)
  xx=sample(N,n)
  rr=mutation_score$MutationAssessor_score[xx]
  rr=as.numeric(as.character(rr))
  y=which(is.na(rr)==FALSE)
  rr=rr[y]
  rr=t(t(rr))
  Rand_score=rbind(Rand_score,rr)
}
wilcox.test(Real_score,Rand_score,alternative = "greater")
library(vioplot)
par(mfrow=c(3,3))
vioplot(Real_score,Rand_score)

#################whether mutations select hub genes in RBP-gene network
setwd("C:/projects/RBP/newidea")
mymut<-read.csv("mymut.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
mutgene=unique(mymut$V2)
ens2sym=read.csv("ensmble2symbol.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
mutgene=t(t(mutgene))
mutens=merge(mutgene,ens2sym,by.x="V1",by.y="HGNC.symbol")

RBPnet<-read.csv("RBP_gene_inter_HepG2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
gene=unique(RBPnet$V2)
genedegree=c()
for (i in 1:length(gene)) {
  x=which(RBPnet$V2==gene[i])
  rr=unique(RBPnet$V1[x])
  genedegree=rbind(genedegree,length(rr))
}
genedegree=cbind(gene,genedegree)
mutdegree=merge(mutens,genedegree,by.x = "Ensembl.Gene.ID",by.y = "gene")
other_g=setdiff(gene,mutens$Ensembl.Gene.ID)
other_g=t(t(other_g))
other_degree=merge(other_g,genedegree,by.x="V1",by.y = "gene")
X1=as.numeric(as.character(mutdegree$V2))
Y1=as.numeric(as.character(other_degree$V2))
wilcox.test(X1,Y1,alternative = "greater")
library(vioplot)
par(mfrow=c(3,3))
vioplot(X1,Y1)

#############KEGG pathway analysis
XX=c(0.00125,0.00341,0.00979,0.00982,0.0114)
par(mfrow=c(3,3))
barplot(-log2(XX[5:1]),horiz = T)

###identify the mutation mediated RBP-target interaction in liver cancer
#####load NormallizeExp.RData

setwd("C:/projects/RBP/newidea")
RBP_target<-read.csv("RBP_target_gene.txt",sep = "\t",header = FALSE,stringsAsFactors = FALSE)
RBP_target_U=RBP_target[!duplicated(RBP_target),]############RBP-protein coding gene interactions
Mutgene<-read.csv("mutation_gene.txt",sep = "\t",header = FALSE,stringsAsFactors = FALSE)
U_mutation<-unique(Mutgene$V2) ###unique mutations in liver cancer
allmutRBP<-c()
expsample=colnames(LIHCmatrix_nor)
expsample=t(t(expsample))
xx=strsplit(expsample[,1],"-")
Nor_Exp_samp=c()
for (i in 1:dim(expsample)[1]) {
  aaa=paste(xx[[i]][1],xx[[i]][2],xx[[i]][3],sep = "-")
  Nor_Exp_samp=rbind(Nor_Exp_samp,aaa)
}
en2symbol=read.csv("ensmble2symbol.txt",sep = "\t",header = T,stringsAsFactors = FALSE)
###for each mutation, find the RBP-taget pairs and then calculated the fc for the target genes in mut and non-mut sample
for (i in 1:length(U_mutation)) {
  
  mut_RBP_target<-c()
  x<-which(Mutgene$V2==U_mutation[i]) ####mutation loc in mutationgene file
  mutsample<-unique(Mutgene$V1[x])
  if (length(intersect(Nor_Exp_samp[,1],mutsample))>0){ #####the mut sample with expression profiles
    print(i)
    ###The column of the mutation samples in expression profile
    mutloc<-match(mutsample,Nor_Exp_samp[,1])
    mut_target=Mutgene$V3[x]###target gene symbol
    #####transfer to ensg
    x1=which(en2symbol$HGNC.symbol==mut_target)
    if(length(x1)>0){
      mutensg=en2symbol$Ensembl.Gene.ID[x1]
      yy1=which(rownames(LIHCmatrix_nor)==mutensg)#########mut gene loc in the expression profile
      if(length(yy1)>0){  #mutated genes are in expression profile
        ###find the target genes amd RBP-binding sites
        x2_1<-which(RBP_target$V4==Mutgene$V3[x[1]]
                    &RBP_target$V10==Mutgene$V7[x[1]]
                    &RBP_target$V8<=Mutgene$V8[x[1]]
                    &RBP_target$V9>=Mutgene$V8[x[1]])
        x2_2<-which(RBP_target$V4==Mutgene$V3[x[1]]
                    &RBP_target$V10==Mutgene$V7[x[1]]
                    &RBP_target$V8<=Mutgene$V9[x[1]]
                    &RBP_target$V9>=Mutgene$V9[x[1]])
        x2<-union(x2_1,x2_2)
        if(length(x2)>0){
          mutexp=LIHCmatrix_nor[yy1,mutloc]
          mutexp=2^mutexp
          nonmutsample<-setdiff(Cancer_s,mutloc)
          nonmutexp<-LIHCmatrix_nor[yy1,nonmutsample]
          nonmutexp=2^nonmutexp
          if (mean(as.numeric(as.character(mutexp)))>1&mean(as.numeric(as.character(nonmutexp)))>1){
            FC<-mean(as.numeric(as.character(mutexp)))/mean(as.numeric(as.character(nonmutexp)))
            for (mm in 1:length(x2)) {
              print(mm)
              mut_RBP_target<-cbind(Mutgene[x,c(1,2,3,11)],RBP_target[x2[mm],],FC)
            }
          }
        }
     
      }
    }
    
  } 
  allmutRBP<-rbind(allmutRBP,mut_RBP_target)
}

write.table(allmutRBP,file = "mutation_RBP_target.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

X<-which(allmutRBP$FC>2)
Y<-which(allmutRBP$FC<0.5)
Z<-union(X,Y)
write.table(allmutRBP[Z,],file = "diff_mutation_RBP_target.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)

diffmut<-allmutRBP[Z,]
bindingIA<-c()
muteffect<-c()
for (i in 1:467) {
  if (diffmut$V6[i]>2){
    bindingIA<-rbind(bindingIA,"inhibit")
    if(diffmut$FC[i]>2){
      muteffect<-rbind(muteffect,"Prom")
    } else {
      muteffect<-rbind(muteffect,"Inbit")
    }
  } else {
    bindingIA<-rbind(bindingIA,"active")
    if(diffmut$FC[i]>2){
      muteffect<-rbind(muteffect,"Prom")
    } else {
      muteffect<-rbind(muteffect,"Inbit")
    }
  }
  
  
}
diffmut=cbind(diffmut,bindingIA,muteffect)
colnames(diffmut)=c("sample","mutloc","targetgene","mutation_type","RBP_ensg","RBP_symbol","Target_ensg","Target_symbol",
                    "target_type","shRNA_FC","peak_loc","peak_start","peak_end","peak_chr","mut_FC","shRNA_effect","mut_effect")
write.table(diffmut,file = "mutation_RBP_target_IA.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)

x1=which(diffmut$shRNA_FC<=0.5&diffmut$mut_FC<=0.5)
x2=which(diffmut$shRNA_FC<0.5&diffmut$mut_FC>=2)
x3=which(diffmut$shRNA_FC>=2&diffmut$mut_FC<=0.5)
x4=which(diffmut$shRNA_FC>=2&diffmut$mut_FC>=2)

xx1=diffmut$FC[x1]-diffmut$V6[x1]
par(mfrow=c(3,3))
pie(c(length(x1),length(x2)))
pie(c(length(x3),length(x4)))

barplot(c(133,76))

xx=union(x1,x4)
satified_net=diffmut[xx,]
write.table(satified_net,file = "mutation_RBP_target_satisfied.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)

##################extract the altered RBP-gene network
RBP=unique(satified_net$RBP_symbol)
tt=unique(satified_net$Target_symbol)
Geneall=union(RBP,tt)

write.table(Geneall,file = "mutation_RBP_target_geneall.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

setwd("C:/projects/RBP/newidea/network")
RBPgene=read.csv("RBPgene.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBP=read.csv("RBP.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBPloc=merge(RBP,RBPgene,by.x = "V1",by.y = "V4")
write.table(RBPloc,file = "rbploc.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
target=setdiff(unique(RBPgene$V4),RBP$V1)
target=t(t(target))
targetloc=merge(target,RBPgene,by.x = "V1",by.y = "V4")
write.table(targetloc,file = "targetloc.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)



#########################clinical analysis
Muteff_sample=unique(satified_net$sample)
Muteff_sample=unique(allmutRBP$V1)
setwd("C:/projects/AS/tcgaCLI")
library(survival)
Mysample=read.csv("LIHC_Mysample.txt",stringsAsFactors=F,sep="\t",skip=0,header = TRUE)
myclass=c()
for (i in 1:dim(Mysample)[1]) {
  xx=which(Muteff_sample==Mysample$V3[i])
  if(length(xx)>0){
    myclass=rbind(myclass,1)
  } else{
    myclass=rbind(myclass,0)
  }
}
Mysample=cbind(Mysample,myclass)
ddd<-survfit(Surv(Mysample$FinalTime/365,Mysample$vital_status=="dead")~Mysample$myclass)
plot(ddd,col = c("green","red"))
fit<-survdiff(Surv(Mysample$FinalTime/365,Mysample$vital_status=="dead")~Mysample$myclass)
fit
install.packages("survminer")
library(survminer)
ggsurvplot(ddd,  size = 1,  # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palettes
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           legend.labs = c("Without", "With"), # Change legend labels
           risk.table.height = 0.15, # Useful to change when you have multiple groups
           ggtheme = theme_bw() # Change ggplot2 theme
)
#############################RBP-lincRNA analysis
setwd("C:/projects/RBP/newidea")
RBPgene=read.csv("RBP_gene_inter_HepG2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
x=which(RBPgene$V3=="lincRNA")
RBPlinc=RBPgene[x,]
RBP=unique(RBPlinc$V1)
lincRNA=unique(RBPlinc$V2)
write.table(RBPlinc,file = "RBP_lincRNA_net.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
cancer_linc=read.csv("cancerlinc.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
lincRNA_degree=c()
for (i in 1:length(lincRNA)) {
  x=which(RBPlinc$V2==lincRNA[i])
  y=which(cancer_linc==lincRNA[i])
  RBP_r=unique(RBPlinc$V1[x])
  lincRNA_degree=rbind(lincRNA_degree,cbind(lincRNA[i],length(RBP_r),length(y)))
}
x1=which(lincRNA_degree[,3]==1)
X1=as.numeric(as.character(lincRNA_degree[x1,2]))
x2=which(lincRNA_degree[,3]==0)
X2=as.numeric(as.character(lincRNA_degree[x2,2]))
boxplot(X1,X2)
wilcox.test(X1,X2)##################no difference in degree of cancer and other lincRNAs
################################################noncoding mutations
noncoding=read.csv("Liver_clean_somatic_mutations_for_signature_analysis.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
noncodingsample=unique(noncoding$V1)################88 samples from liver cancer
noncodingsample=t(noncodingsample)
write.table(noncodingsample,file = "liver_noncoding_sample.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
Biofunction <- file("CosmicNCV.tsv", "r")
line=readLines(Biofunction,n=1)
Result <- c()
i=1;
while( length(line) != 0 ) {
  print(i)
  cosline <- strsplit(line,"\t")
  cossam <- cosline[[1]][1]
  x=which(noncodingsample==cossam)
  if(length(x)>0){
    Result=rbind(Result,cosline[[1]][1:18])
  }
  i=1+i
  line=readLines(Biofunction,n=1);
}
close(Biofunction)
write.table(Result,file = "liver_noncoding_mutation.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
#################revised noncoding mutation
noncoding=read.csv("Liver_clean_somatic_mutations_for_signature_analysis.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
noncodingsample=unique(noncoding$V1)################88 samples from liver cancer
noncodingsample=t(noncodingsample)
noncodingsample=t(noncodingsample)
noncoding=read.csv("CosmicNCV.tsv",stringsAsFactors=F,sep="\t",skip=0,header = T)
ncmutation=merge(noncodingsample,noncoding,by.x="V1",by.y="Sample.name")
write.table(ncmutation,file = "liver_noncoding_mutation.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
###########################ananlysis of RBP-lincRNA based on expression
setwd("C:/projects/RBP/newidea/RBP_lincRNA")
RBP_lincRNA=read.csv("RBP_lincRNA.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
lincRNA=unique(RBP_lincRNA$V3)
allgene=rownames(LIHCmatrix_nor)
XX=intersect(allgene,lincRNA)
#############identified the differentially expressed lincRNAs
FC=c()
P_value=c()
for (i in 1:length(XX)) {
  x=which(allgene==XX[i])
  ca_exp=LIHCmatrix_nor[x,Cancer_s]
  na_exp=LIHCmatrix_nor[x,Normal_s]
  ca_exp=as.numeric(as.character(ca_exp))
  na_exp=as.numeric(as.character(na_exp))
  FC=rbind(FC,mean(2^ca_exp)/mean(2^na_exp))
  P_value=rbind(P_value,wilcox.test(ca_exp,na_exp)$p.value)
}
FDR=p.adjust(P_value,method = "fdr")
lincRNA_DF=cbind(XX,FC,P_value,FDR)
xx=which(FDR<0.01&FC>2)
yy=which(FDR<0.01&FC<0.5)
zz=union(xx,yy)
diffLncRNA=lincRNA_DF[zz,]
Diff_RBP_linc=merge(RBP_lincRNA,diffLncRNA,by.x="V3",by.y="XX")
write.table(Diff_RBP_linc,file = "Diff_interaction.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

##################ENSG00000248323 expression  ENSG00000225680
x=which(allgene=="ENSG00000212694")
ca_exp=LIHCmatrix_nor[x,Cancer_s]
na_exp=LIHCmatrix_nor[x,Normal_s]
ca_exp=as.numeric(as.character(ca_exp))
na_exp=as.numeric(as.character(na_exp))
par(mfrow=c(3,3))
library(vioplot)
vioplot(na_exp,ca_exp)
wilcox.test(na_exp,ca_exp)
###############linc01089
x=which(allgene=="ENSG00000212694")
ca_exp=LIHCmatrix_nor[x,Cancer_s]
na_exp=LIHCmatrix_nor[x,Normal_s]
ca_exp=as.numeric(as.character(ca_exp))
Allsample=colnames(LIHCmatrix_nor)
Allcancersam=Allsample[Cancer_s]
KK=cbind(Allcancersam,ca_exp)
Threesampe=c()
for (i in 1:dim(KK)[1]) {
  aa=strsplit(KK[i,1],"-")
  bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
  Threesampe=rbind(Threesampe,bb)
}
KK=cbind(KK,Threesampe)
Sa<-read.csv("LIHC_Mysample.txt",stringsAsFactors=F,sep="\t",skip=0,header = TRUE)
Lincsample=merge(KK,Sa,by.x="V3",by.y="V3")
Lincsample$XX=1
aa=as.numeric(as.character(Lincsample$ca_exp))
bb=which(aa<median(aa))
# Q75=quantile(aa,0.70)[[1]]
# Q25=quantile(aa,0.30)[[1]]
Lincsample$XX[bb]=0
# bb=which(aa>Q75)
# Lincsample$XX[bb]=2
# bb=which(aa<Q25)
# Lincsample$XX[bb]=0
# cc=which(Lincsample$XX==0|Lincsample$XX==2)
library(survival)
ddd<-survfit(Surv(Lincsample$FinalTime,Lincsample$vital_status=="dead")~Lincsample$XX)
plot(ddd,col = c("green","red"))
fit<-survdiff(Surv(Lincsample$FinalTime,Lincsample$vital_status=="dead")~Lincsample$XX)
fit

library(survminer)
ggsurvplot(ddd,  size = 1,  # change line size
           palette = c("#E7B800", "#2E9FDF"), # custom color palettes
           conf.int = TRUE, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = TRUE, # Add risk table
           risk.table.col = "strata", # Risk table color by groups
           legend.labs = c("Without", "With"), # Change legend labels
           risk.table.height = 0.15, # Useful to change when you have multiple groups
           ggtheme = theme_bw() # Change ggplot2 theme
)



xx1=which(Lincsample$tumor_stage=="stage i")
xx2=which(Lincsample$tumor_stage=="stage ii")
xx3=which(Lincsample$tumor_stage=="stage iii"|Lincsample$tumor_stage=="stage iiia"|
            Lincsample$tumor_stage=="stage iiib"|Lincsample$tumor_stage=="stage iiic")
xx4=which(Lincsample$tumor_stage=="stage iv"|Lincsample$tumor_stage=="stage iva"|
            Lincsample$tumor_stage=="stage ivb")
boxplot(as.numeric(Lincsample$ca_exp[xx1]),as.numeric(Lincsample$ca_exp[xx2]),as.numeric(Lincsample$ca_exp[xx3]),as.numeric(Lincsample$ca_exp[xx4]))


wilcox.test(as.numeric(Lincsample$ca_exp[xx1]),as.numeric(Lincsample$ca_exp[xx4]))

xx1=which(Lincsample$XX==0)
xx2=which(Lincsample$XX==1)
X1=Lincsample$bmi[xx1]
X2=Lincsample$bmi[xx2]
X1=as.numeric(as.character(X1))
X2=as.numeric(as.character(X2))
wilcox.test(X1,X2)

########################RBP-lincRNA-mutation analysis

setwd("C:/projects/RBP/newidea/RBP_lincRNA")
RRBP_bindsite=read.csv("Target_bed_ensg.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
BP_lincRNA=read.csv("Diff_interaction.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBP_lincRNA_site=merge(RBP_lincRNA,RBP_bindsite,by.x = c("V3","V1"),by.y = c("V2","V1"))
write.table(RBP_lincRNA_site,file = "Diff_lincRNA_site.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

#####################Mutated target gene pathway enrichment
KEGG=read.csv("mutgeneKEGG.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
x1=which(KEGG$PValue<0.01)
x2=which(KEGG$PValue>0.01)
par(mfrow=c(3,3))
plot(KEGG$Fold.Enrichment[x2],-log10(as.numeric(KEGG$PValue[x2])),ylim = c(0.5,3.5),col="gray")
points(KEGG$Fold.Enrichment[x1],-log10(as.numeric(KEGG$PValue[x1])),col="red")
abline(h=-log10(0.01),col="lightgray",lty=2)

##################################################
#########rbp network and PPI network integration analysis
setwd("C:/projects/RBP/newidea")
RBP_net=read.csv("Net.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
PPI=read.csv("HI3_Y2H_102416.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
PPI_gene1=unique(PPI$Symbol.for.A)
PPI_gene2=unique(PPI$Symbol.for.B)
PPI_gene=union(PPI_gene1,PPI_gene2)
#######PPI_degree
PPI_degree=c()
for (i in 1:length(PPI_gene)) {
  x1=which(PPI$Symbol.for.A==PPI_gene[i])
  x2=which(PPI$Symbol.for.B==PPI_gene[i])
  x=union(x1,x2)
  PPI_degree=rbind(PPI_degree,cbind(PPI_gene[i],length(x)))
}
#############rbp network indegree and outdegree
RBP=unique(RBP_net$V2)
target_gene=unique(RBP_net$V4)
outdegree=c()
indegree=c()
for (i in 1:length(RBP)) {
  x=which(RBP_net$V2==RBP[i])
  gg=unique(RBP_net$V4[x])
  outdegree=rbind(outdegree,cbind(RBP[i],length(gg)))
}
for (i in 1:length(target_gene)) {
  x=which(RBP_net$V4==target_gene[i])
  gg=unique(RBP_net$V2[x])
  indegree=rbind(indegree,cbind(target_gene[i],length(gg)))
}
in_ppi=merge(indegree,PPI_degree,by.x="V1",by.y="V1")
out_ppi=merge(outdegree,PPI_degree,by.x="V1",by.y="V1")
x1=as.numeric(as.character(in_ppi$V2.x))
x2=as.numeric(as.character(in_ppi$V2.y))
RR=cor(x1,x2)
cor.test(x1,x2)###########R=0.042,P=0.06
x3=as.numeric(as.character(out_ppi$V2.x))
x4=as.numeric(as.character(out_ppi$V2.y))
RR=cor(x3,x4)
cor.test(x3,x4) ###########R=0.26,p=0.17
###########RBP-RBP-coregulation
RBP_coregulation=c()
for (i in 1:length(RBP)) {
  for (j in 1:length(RBP)) {
    x1=which(RBP_net$V2==RBP[i])
    x2=which(RBP_net$V2==RBP[j])
    target1=unique(RBP_net$V4[x1])
    target2=unique(RBP_net$V4[x2])
    nn=length(intersect(target1,target2))
    sco=nn/min(length(target1),length(target2))
    RBP_coregulation=rbind(RBP_coregulation,cbind(RBP[i],RBP[j],nn,sco))
  }
}
KK_sim=as.numeric(as.character(RBP_coregulation[,4]))
KK_matrix=matrix(KK_sim,nrow = 65)
rownames(KK_matrix)=RBP
colnames(KK_matrix)=RBP
library(pheatmap)
pheatmap(KK_matrix)
#############divided into low-coregulation, high-coregulation
highco=read.csv("high_co.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
lowco=read.csv("low_co.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
High_PPI=c()
for (i in 1:45) {
  for (j in 1:45) {
    x1=which(PPI$Symbol.for.A==highco$V1[i])
    x2=which(PPI$Symbol.for.B==highco$V1[i])
    target1=unique(PPI$Symbol.for.B[x1])
    target2=unique(PPI$Symbol.for.A[x2])
    target_A=union(target1,target2)
    
    x1=which(PPI$Symbol.for.A==highco$V1[j])
    x2=which(PPI$Symbol.for.B==highco$V1[j])
    target1=unique(PPI$Symbol.for.B[x1])
    target2=unique(PPI$Symbol.for.A[x2])
    target_B=union(target1,target2)
    
    nn=length(intersect(target_A,target_B))
    sco=nn/min(length(target_A),length(target_B))
    High_PPI=rbind(High_PPI,cbind(highco$V1[i],highco$V1[j],nn,sco))
  }
}

Low_PPI=c()
for (i in 1:20) {
  for (j in 1:20) {
    x1=which(PPI$Symbol.for.A==highco$V1[i])
    x2=which(PPI$Symbol.for.B==highco$V1[i])
    target1=unique(PPI$Symbol.for.B[x1])
    target2=unique(PPI$Symbol.for.A[x2])
    target_A=union(target1,target2)
    
    x1=which(PPI$Symbol.for.A==highco$V1[j])
    x2=which(PPI$Symbol.for.B==highco$V1[j])
    target1=unique(PPI$Symbol.for.B[x1])
    target2=unique(PPI$Symbol.for.A[x2])
    target_B=union(target1,target2)
    
    nn=length(intersect(target_A,target_B))
    sco=nn/min(length(target_A),length(target_B))
    Low_PPI=rbind(Low_PPI,cbind(lowco$V1[i],lowco$V1[j],nn,sco))
  }
}
x1=which(High_PPI[,4]!="NaN")
x2=which(Low_PPI[,4]!="NaN")
AA1=as.numeric(as.character(High_PPI[x1,4]))
BB1=as.numeric(as.character(Low_PPI[x2,4]))
library(vioplot)
vioplot(AA1,BB1)
wilcox.test(AA1,BB1)
hist(AA1,10)
hist(BB1,10)#############no significant result
######################################################CDKN1A survival
setwd("C:/projects/RBP/newidea")
CDKN1A=read.csv("CDKN1A_sur.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
library(survival)
fit<-survdiff(Surv(CDKN1A$V1,CDKN1A$V2=="dead")~CDKN1A$V3)
fit1=survfit(Surv(CDKN1A$V1,CDKN1A$V2=="dead")~CDKN1A$V3)
plot(fit1)

################################predict the RBP-ORF interactions
###############get the unique sequence of ORF
setwd("C:/projects/RBP/predictRBP")
ORF=read.csv("FULL-ORFeome-011217.txt",stringsAsFactors=F,sep="\t",skip=0,header = TRUE)
U_gene=unique(ORF$ENTREZ_GENE_ID)
Gene_seq=c()
xx=which(U_gene!="NULL")
U_gene=U_gene[xx]
for(i in 1:length(U_gene)){
  print(i)
  xx=which(ORF$ENTREZ_GENE_ID==U_gene[i])
  Gseq=ORF$cds_seq[xx]
  Lseq=c()
  for(j in 1:length(xx)){
    Lseq=rbind(Lseq,nchar(Gseq[j]))
  }
  aa=which(Lseq==max(Lseq))
  aa=aa[1]
  Myseq=Gseq[aa]
  Gene_seq=rbind(Gene_seq,cbind(U_gene[i],Myseq))
}
write.table(Gene_seq,file = "ORF_Seq.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
############################RBP mutation pan-cancer analysis
######average expression and var comparsion
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
FC_mean=c()
FC_var=c()
P_mean=c()
P_var=c()
for(kc in 1:33){
  print(kc)
  setwd("C:/projects/AS/geneexpression")
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep=""))
  
  library(SummarizedExperiment)
  LIHCmatrix<-assay(data,1,"FPKM")############gene exprssion
  ###################filter the gene expression profile
  N=dim(LIHCmatrix)[2]
  filter_line<-c()
  for (i in 1:dim(LIHCmatrix)[1]) {
    n=which(LIHCmatrix[i,]<0.1)
    if(length(n)>N*0.3){
      filter_line<-rbind(filter_line,i)
    }
  }
  LIHCmatrix<-LIHCmatrix[-filter_line,]
  LIHCclinal<-data@colData@listData$definition#############sample information tumor vs normal
  LIHCclinal=t(t(LIHCclinal))
  setwd("C:/projects/RBP/Pancancer")
  RBP<-read.csv("EnsemRBP.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
  LIHCmatrix_nor=log2(LIHCmatrix+0.05)
  M_exp=c()
  V_exp=c()
  for (i in 1:dim(LIHCmatrix_nor)[1]) {
    #print(i)
    M_exp=rbind(M_exp,mean(LIHCmatrix_nor[i,]))
    V_exp=rbind(V_exp,var(LIHCmatrix_nor[i,]))
  }
  GeneMV=cbind(rownames(LIHCmatrix),M_exp,V_exp)
  RBPMV=merge(RBP,GeneMV,by.x="V1",by.y="V1")
  othergene=setdiff(GeneMV[,1],RBPMV[,1])
  othergene=t(t(othergene))
  OtherMV=merge(othergene,GeneMV,by.x="V1",by.y="V1")
  X1=as.numeric(as.character(RBPMV[,6]))
  X2=as.numeric(as.character(OtherMV$V2))
  X3=as.numeric(as.character(RBPMV[,7]))
  X4=as.numeric(as.character(OtherMV$V3))
  library(vioplot)
  pdf(paste(cancer[kc],".pdf",sep = ""))
  par(mfrow=c(3,3))
  vioplot(X1,X2)
  vioplot(X3,X4)
  dev.off()
  FC1=mean(X1)/mean(X2)
  FC2=mean(X3)/mean(X4)
  P1=wilcox.test(X1,X2,alternative = "greater")$p.value
  P2=wilcox.test(X3,X4,alternative = "less")$p.value#########p<2.2e-16 RBPs have higher expression and lower varation
  FC_mean=rbind(FC_mean,FC1)
  FC_var=rbind(FC_var,FC2)
  P_mean=rbind(P_mean,P1)
  P_var=rbind(P_var,P2)
}

##########################Mutation sample number
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
Sample_ratio=c()
for(kc in 1:33){
  print(kc)
  setwd("C:/projects/RBP/newidea")
  RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
  setwd("C:/projects/TCGAmutation")
  mutation<-read.csv(paste(cancer[kc],".maf",sep = ""),stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
  N=length(unique(mutation$Tumor_Sample_Barcode))
  xx=merge(mutation,RBP,by.x="Hugo_Symbol",by.y="HGNC.symbol")
  n=length(unique(xx$Tumor_Sample_Barcode))
  Sample_ratio=rbind(Sample_ratio,n/N)
}
setwd("C:/projects/RBP/Pancancer")

####################CADD and conservation
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
CADD_FC=c()
CADD_P=c()
Conser_FC=c()
Conser_P=c()

for(kc in 1:33){
  print(kc)
  setwd("C:/projects/RBP/newidea")
  RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
  setwd("C:/projects/TCGAmutation")
  mutation<-read.csv(paste(cancer[kc],".maf",sep = ""),stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
  setwd("C:/projects/RBP/Pancancer/mut")
  mutation_score<-read.csv(paste("TCGA.",cancer[kc],".maf.hg38_multianno.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header =T)
  RBP_mut<-merge(RBP,mutation,by.x = "HGNC.symbol",by.y = "Hugo_Symbol")
  RBPmut_score<-merge(RBP_mut,mutation_score,by.x = c("Chromosome","Start_Position","End_Position","Tumor_Seq_Allele1","Tumor_Seq_Allele2"),
                      by.y = c("Chr","Start","End","Ref","Alt"))
  Real_score=RBPmut_score$CADD_phred
  Real_score=as.numeric(as.character(Real_score))
  X=which(is.na(Real_score)==FALSE)
  Real_score=Real_score[X]
  n=length(Real_score)
  N=dim(mutation_score)[1]
  Rand_score=c()
  for (i in 1:1000) {
   # print(i)
    xx=sample(N,n)
    rr=mutation_score$CADD_phred[xx]
    rr=as.numeric(as.character(rr))
    y=which(is.na(rr)==FALSE)
    rr=rr[y]
    rr=t(t(rr))
    Rand_score=rbind(Rand_score,rr)
  }
  P1=wilcox.test(Real_score,Rand_score,alternative = "greater")$p.value
  FC1=mean(Real_score)/mean(Rand_score)
  CADD_P=rbind(CADD_P,P1)
  CADD_FC=rbind(CADD_FC,FC1)
  library(vioplot)
  setwd("C:/projects/RBP/Pancancer")
  pdf(paste(cancer[kc],"_cadd.pdf",sep = ""))
  par(mfrow=c(3,3))
  vioplot(Real_score,Rand_score)
  ####conservation
  Real_score=RBPmut_score$integrated_fitCons_score
  Real_score=as.numeric(as.character(Real_score))
  X=which(is.na(Real_score)==FALSE)
  Real_score=Real_score[X]
  n=length(Real_score)
  N=dim(mutation_score)[1]
  Rand_score=c()
  for (i in 1:100) {
    #print(i)
    xx=sample(N,n)
    rr=mutation_score$integrated_fitCons_score[xx]
    rr=as.numeric(as.character(rr))
    y=which(is.na(rr)==FALSE)
    rr=rr[y]
    rr=t(t(rr))
    Rand_score=rbind(Rand_score,rr)
  }
  P2=wilcox.test(Real_score,Rand_score,alternative = "greater")$p.value
  FC2=mean(Real_score)/mean(Rand_score)
  Conser_P=rbind(Conser_P,P2)
  Conser_FC=rbind(Conser_FC,FC2)
  vioplot(Real_score,Rand_score)
  dev.off()
}
##############surface analysis
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
R_surface=c()
P1_surface=c()
P2_surface=c()
setwd("C:/projects/ORFchar")
# subdir <- list.files(pattern="result",recursive=T)
# write.table(subdir,file = "file.txt",sep = "\t",quote = F,row.names = F,col.names = F)
Surfscore=read.csv("Surf_score.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)

for(kc in 2:33){
  setwd("C:/projects/ORFchar")
  TCGA=read.csv("TCGA-hg38-orfeome.tsv",stringsAsFactors = FALSE,sep = "\t",header = T)
  xx=which(TCGA$cancer==cancer[kc])
  LIHC_orf=TCGA[xx,]
  LIHC_surf=merge(LIHC_orf,Surfscore,by.x=c("ORF_ID","residue_."),by.y = c("V3","V4"))
  setwd("C:/projects/TCGAmutation")
  mutation<-read.csv(paste(cancer[kc],".maf",sep=""),stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
  mutation_surf<-merge(mutation,LIHC_surf,by.x = c("Chromosome","Start_Position","HGVSc"),by.y = c("chromosome","position","HGVSc"))
  setwd("C:/projects/RBP/newidea")
  RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
  RBPmut_score<-merge(RBP,mutation_surf,by.x = "HGNC.symbol",
                      by.y = "Hugo_Symbol")
  Real_score=RBPmut_score$V5
  Real_score=as.numeric(as.character(Real_score))
  X=which(is.na(Real_score)==FALSE)
  Real_score=Real_score[X]
  n=length(Real_score)
  
  gene=unique(mutation_surf$Hugo_Symbol)
  othergene=setdiff(gene,RBP$HGNC.symbol)
  othergene=t(t(othergene))
  otherscore=merge(mutation_surf,othergene,by.x = "Hugo_Symbol",by.y = "V1")
  N=dim(otherscore)[1]
  Rand_score=c()
  KK=c()
  for (i in 1:100) {
    print(i)
    xx=sample(N,n)
    rr=otherscore$V5[xx]
    rr=as.numeric(as.character(rr))
    y=which(is.na(rr)==FALSE)
    rr=rr[y]
    rr=t(t(rr))
    Rand_score=rbind(Rand_score,rr)
    KK=rbind(KK,sum(rr>0.25))
  }
  NN=length(which(Real_score>0.25))
  N_all=length(Real_score)
  R_surface=rbind(R_surface,NN/N_all)
  P1=wilcox.test(Real_score,Rand_score,alternative = "greater")$p.value
  P1_surface=rbind(P1_surface,P1)
  N2=length(which(KK>NN))
  P2=N2/1000
  P2_surface=rbind(P2_surface,P2)
  setwd("C:/projects/RBP/Pancancer")
  pdf(paste(cancer[kc],"_surf.pdf",sep = ""))
  par(mfrow=c(3,3))
  hist(KK,20, prob=TRUE) 
  lines(density(KK))
  points(NN,0)
  dev.off()
}

######################predicted RBP-ORF regulations
##############Running on the cluster
setwd("C:/projects/RBP/predictRBP/RBP")
motifid=read.csv("ATtRACT_db.txt",stringsAsFactors=F,sep="\t",skip=0,header = TRUE)
motifid=cbind(motifid$Gene_name,motifid$Organism,motifid$Matrix_id)
motifid=motifid[!duplicated(motifid),]
xx=which(motifid[,2]=="Homo_sapiens")
Human=motifid[xx,]###########1198 motif matrix ID; 160 RBPs
RBP_T=read.csv("RBP_result4.txt",stringsAsFactors=F,sep="\t",skip=0,header = TRUE)
Human_RBP=merge(RBP_T,Human,by.x="X..motif_id",by.y="V3")
Human_RBP_Target=cbind(as.character(Human_RBP$V1),Human_RBP$sequence_name,Human_RBP$start,Human_RBP$stop,Human_RBP$p.value)
colnames(Human_RBP_Target)=c("RBP_name","Target","start","end","P")
write.table(Human_RBP_Target,"RBP_gene_binding.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)


RBP_gene=c()
Biofunction <- file("RBP_result4.txt", "r")
line=readLines(Biofunction,n=1)
i=0
while( length(line) != 0 ) {
  i <- i+1
  print(i)
  msig.go <- strsplit(line,"\t")
  tt=msig.go[[1]]
  ID=msig.go[[1]][1]
  xa=which(Human[,3]==ID)
  if(length(xa)>0){
    RT=cbind(Human[xa,1],tt[3],tt[4],tt[5],tt[8])
    RBP_gene=rbind(RBP_gene,RT)
  }
  line=readLines(Biofunction,n=1);
}
close(Biofunction)
colnames(RBP_gene)=c("RBP_gene","Entrez","start","end","p_value")
write.table(RBP_gene,"RBP_gene_binding.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)


setwd("C:/projects/RBP/newidea")
HepG2=read.csv("RBP_gene_inter_HepG2.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBP.cell=unique(HepG2$V1) ####22 RBP overlap with prediction 22-65-160
##################################
##################expression filters
##################correlation coefficient
install.packages("devtools")
devtools::install_github("stephenturner/annotables")
library(dplyr)
library(annotables)
xx=grch38
pg=which(xx$biotype=="protein_coding")
Protein.gene=cbind(xx$ensgene[pg],xx$entrez[pg],xx$symbol[pg])
colnames(Protein.gene)=c("Ensg","Entrez","Symbol")
write.table(Protein.gene,"Protein_gene.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)

############### for each cancer type, compute the correlation coefficient
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
setwd("C:/projects/RBP/predictRBP/RBP")
RBP.target=read.csv("RBP_gene_binding.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
RBP.target=RBP.target[!duplicated(RBP.target[,1:2]),1:2]
Pgene=read.csv("Protein_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
for(kc in 1:4){
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep = ""))
  library(SummarizedExperiment)
  cancer.matrix<-assay(data,1,"FPKM")
  C_N=data@colData@listData$definition#############sample information tumor vs normal
  Cancer_s=which(C_N!="Solid Tissue Normal")
  cancer.matrix=cancer.matrix[,Cancer_s]
  allgene=rownames(cancer.matrix)
  R_corr_cancer=c()
  for(i in 1:dim(RBP.target)[1]){
    print(i)
    RBP.ensg=which(Pgene$Symbol==RBP.target$RBP_gene[i])
    T.ensg=which(Pgene$Entrez==RBP.target$Entrez[i])
    if(length(RBP.ensg)>0&&length(T.ensg)>0){
      xxa=which(allgene==Pgene$Ensg[RBP.ensg])
      xxb=which(allgene==Pgene$Ensg[T.ensg])
      if(length(xxa)>0&&length(xxb)>0){
        expa=cancer.matrix[xxa[1],]
        expb=cancer.matrix[xxb[1],]
        expa=log2(as.numeric(as.character(expa)))
        expb=log2(as.numeric(as.character(expb)))
        R=cor(expa, expb, method = "spearman")
        pp=cor.test(expa, expb, method = "spearman")$p.value
        paircor=cbind(Pgene$Ensg[RBP.ensg],Pgene$Ensg[T.ensg],RBP.target$RBP_gene[i],
                      RBP.target$Entrez[i],Pgene$Symbol[T.ensg],R,pp)
        R_corr_cancer=rbind(R_corr_cancer,paircor)
      } 
    }
  }
  FDR=p.adjust(R_corr_cancer[,7],method = "BH")
  R_corr_cancer=cbind(R_corr_cancer,FDR)
  colnames(R_corr_cancer)=c("RBP_ensg","Tar_ensg","RBP_sym","Tar_entrez","Tar_sym","R","P","FDR")
  write.table(R_corr_cancer,paste(cancer[kc],"_RBP_corr.txt"),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
  xx=which(R_corr_cancer[,8]<0.01)
  Sig_R=R_corr_cancer[xx,]
  write.table(Sig_R,paste(cancer[kc],"_RBP_corr_sig.txt"),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
}

###################################################
###################################################
####RBP network analysis in 33 types of cancer
#########network components ratio
###data structure: LGG, RBP, 0.8; LGG, Gene, 0.7
setwd("C:/projects/RBP/Pancancer/RBP_network")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
Network.pro=c()
N.net=751539
N.RBP=125
N.gene=15429
Net_r=c()
for(kc in 1:33){
  print(kc)
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,3:4]),3:4]
  RBP.gene=unique(RBP.target$RBP_sym)
  T.gene=unique(RBP.target$Tar_entrez)
  n1=dim(RBP.target)[1]
  n2=length(RBP.gene)
  n3=length(T.gene)
  C1=cbind(kc,"cnet",n1/N.net)
  C2=cbind(kc,"RBP",n2/N.RBP)
  C3=cbind(kc,"Gene",n3/N.gene)
  Network.pro=rbind(Network.pro,C1)
  Network.pro=rbind(Network.pro,C2)
  Network.pro=rbind(Network.pro,C3)
  Net_r=rbind(Net_r,cbind(n1,n2,n3,n1/N.net,n2/N.RBP,n3/N.gene))
}
barplot(Net_r[,4])
colnames(Network.pro)=c("ctype","comtype","Pro")
Network.pro1=data.frame(Pro=as.numeric(Network.pro[,3]),ctype=as.numeric(as.character(Network.pro[,1])),comtype=Network.pro[,2])
library(latticeExtra)
cloud(Pro~comtype+ctype, Network.pro1, panel.3d.cloud=panel.3dbars, col.facet='lightblue', 
      xbase=0.3, ybase=0.3, scales=list(arrows=FALSE, col=1), 
      par.settings = list(axis.line = list(col = "transparent")))
xx=which(Network.pro1$comtype=="cnet")
library(plotly)
plot_ly(x=Network.pro1$ctype[xx],y=Network.pro1$Pro[xx],type = 'scatter', mode = 'lines',fill = 'tozeroy')
#####################################################RBP_network recurrent and cancer similarity
RBP.net.all=read.csv("RBP_gene_binding.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
RBP.net.all=RBP.net.all[!duplicated(RBP.net.all[,1:2]),1:2]
Reg01=c()
for(kc in 1:33){
  print(kc)
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,3:4]),3:4]
  K=rep(0,dim(RBP.net.all)[1])
  WW=match(interaction(RBP.net.all$RBP_gene, RBP.net.all$Entrez), 
           interaction(RBP.target$RBP_sym, RBP.target$Tar_entrez))
  aa=which(is.na(WW)==F)
  K[aa]=1
  Reg01=cbind(Reg01,K)
}  

Reg.number=apply(Reg01,1,sum)
Cancer.number=c()
for(i in 1:33){
  xx=which(Reg.number==i)
  Cancer.number=rbind(Cancer.number,length(xx))
}
sum(Cancer.number)
Cancer.sim=c()
for(i in 1:33){
  Sim1=c()
  for(j in 1:33){
    x1=which(Reg01[,i]==1)
    x2=which(Reg01[,j]==1)
    n1=length(x1)
    n2=length(x2)
    x3=intersect(x1,x2)
    x4=union(x1,x2)
    Sim1=rbind(Sim1,length(x3)/min(n1,n2))
  }
  Cancer.sim=cbind(Cancer.sim,Sim1)
}
KS=c(
  151635,
  377494,
  59984,
  10063,
  2440
)
pie(KS)
xx=which(Reg.number==33)
Recure.reg=RBP.net.all[xx,]
library(pheatmap)
pheatmap(Cancer.sim,show_colnames = T,labels_col = cancer,clustering_distance_rows = "correlation",clustering_distance_cols  = "correlation")
colnames(Cancer.sim)=cancer

hc = hclust(dist(t(Cancer.sim)),method = "average")
library(ape)
plot(as.phylo(hc), type = "fan")





##############################RBP-network analysis
setwd("C:/projects/RBP/Pancancer/RBP_network")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
HI3=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
All.coreg=c()
par(mfrow=c(3,3))
p.int=c()############whether the interacting RBP co-regulation
P.all=c()
for(kc in 1:33){
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,3:4]),3:4]
  RBP.gene=unique(RBP.target$RBP_sym)
  T.gene=unique(RBP.target$Tar_entrez)
  # RBP.degree=c()
  # for(i in 1:length(RBP.gene)){
  #   print(i)
  #   xx=which(RBP.target$RBP_sym==RBP.gene[i])
  #   RBP.degree=rbind(RBP.degree,length(xx))
  # }
  # T.degree=c()
  # for(i in 1:length(T.gene)){
  #   print(i)
  #   xx=which(RBP.target$Tar_entrez==T.gene[i])
  #   T.degree=rbind(T.degree,length(xx))
  # }
  # RBP.degree=cbind(RBP.gene,RBP.degree)
  # T.degree=cbind(T.gene,T.degree)
  ###############################
  ######co-regulation analysis
  RR_g=c()
  n=length(RBP.gene)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      print(c(i,j))
      xx1=which(RBP.target$RBP_sym==RBP.gene[i])
      xx2=which(RBP.target$RBP_sym==RBP.gene[j])
      TT1=unique(RBP.target$Tar_entrez[xx1])
      TT2=unique(RBP.target$Tar_entrez[xx2])
      OT=intersect(TT1,TT2)
      UT=union(TT1,TT2)
      jad=length(OT)/length(UT)
      R=cbind(RBP.gene[i],RBP.gene[j],jad)
      RR_g=rbind(RR_g,R)
    }
  }
  # plot(density(as.numeric(RR_g[,3])),col="red",ylim=c(0,40))
  # N=dim(RBP.target)[1]
  # xx=sample(c(1:N),N)
  # RBP.target.r=cbind(RBP.target$RBP_sym,RBP.target$Tar_entrez[xx])
  # RR_g_r=c()
  # n=length(RBP.gene)
  # for(i in 1:(n-1)){
  #   for(j in (i+1):n){
  #     print(c(i,j))
  #     xx1=which(RBP.target.r[,1]==RBP.gene[i])
  #     xx2=which(RBP.target.r[,1]==RBP.gene[j])
  #     TT1=unique(RBP.target.r[xx1,2])
  #     TT2=unique(RBP.target.r[xx2,2])
  #     OT=intersect(TT1,TT2)
  #     UT=union(TT1,TT2)
  #     jad=length(OT)/length(UT)
  #     R=cbind(RBP.gene[i],RBP.gene[j],jad)
  #     RR_g_r=rbind(RR_g_r,R)
  #   }
  # }
  # lines(density(as.numeric(RR_g_r[,3])),col="blue")
  # p.value=ks.test(as.numeric(RR_g[,3]),as.numeric(RR_g_r[,3]),alternative = "greater")$p.value
  # P.all=rbind(P.all,p.value)
  I.lab=c()#########whether the RBP-RBP interaction
  for(i in 1:dim(RR_g)[1]){
    x1=which(HI3$V1==RR_g[i,1]&HI3$V2==RR_g[i,2])
    x2=which(HI3$V2==RR_g[i,1]&HI3$V1==RR_g[i,2])
    if(length(x1)>0|length(x2)>0){
      I.lab=rbind(I.lab,1)
    } else {
      I.lab=rbind(I.lab,0)
    }
  }
  boxplot(as.numeric(RR_g[,3])~I.lab,outline = F)
  stripchart(RR_g[,3] ~ I.lab, vertical = TRUE, data = RR_g,
             method = "jitter", add = TRUE, pch = 3, col = c('gray','lightblue'))
  X1=which(I.lab==1)
  X2=which(I.lab==0)
  pp=wilcox.test(as.numeric(as.character(RR_g[X1,3])),as.numeric(RR_g[X2,3]),alternative = "greater")$p.value
  p.int=rbind(p.int,pp)
  aa=which(I.lab==1)
  All.coreg=rbind(All.coreg,RR_g[aa,])
}
xx=which(All.coreg[,3]>0.2)
write.table(All.coreg[xx,],"Pan_RBP_co.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
######################################
###RBP-RBP motif analysis
setwd("C:/projects/RBP/Pancancer/RBP_network")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")

ACTA=read.csv("C:/projects/RBP/predictRBP/RBP/ATtRACT_db.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
x=which(ACTA$Organism=="Homo_sapiens")
Human.RBPs=unique(ACTA$Gene_name[x])
Human.RBPs=t(t(Human.RBPs))
All.RBP.RBP=c()
for(kc in 1:33){
  print(kc)
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,3:5]),3:5]
  xrbp=merge(RBP.target,Human.RBPs,by.x = "Tar_sym","V1")
  All.RBP.RBP=rbind(All.RBP.RBP,xrbp)
}
All.RBP.RBP=All.RBP.RBP[!duplicated(All.RBP.RBP[,c(2,1)]),c(2,1)]
write.table(All.RBP.RBP,"rbp_rbp_reg.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
gene1=unique(All.RBP.RBP$RBP_sym)
gene2=unique(All.RBP.RBP$Tar_sym)
gene=union(gene1,gene2)
write.table(gene,"rbp_ID.txt",sep = "\t",row.names = TRUE, col.names = F,quote = FALSE)
#############delete the self regulation
RBP.RBP=read.csv("RBP_RBP_id.txt",stringsAsFactors=F,sep="\t",skip=0,header = FALSE)
xx=which(RBP.RBP$V1==RBP.RBP$V2)
RBP.RBP=RBP.RBP[-xx,]
write.table(RBP.RBP,"rbp_ID_NOSELF.txt",sep = "\t",row.names = F, col.names = F,quote = FALSE)

Z.score=c(-4.21,1.45,1.29,-4.25,3.65,0.76,-1.47,-1.04,-3.86,-3.67,4.05,0,0.97)
Z.score.n=(Z.score-mean(Z.score))/var(Z.score)
par(mfrow=c(3,3))
plot(Z.score.n,type = "l")
points(Z.score.n,col="red")
####################################whether the cancer genes were regulated by more RBPs
setwd("C:/projects/RBP/Pancancer/RBP_network")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
core.fit=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
CGC=read.csv("Census_allSun Nov 13 17-37-13 2016.tsv",stringsAsFactors=F,sep="\t",skip=0,header = T)
HI3=read.csv("SIN.csv",stringsAsFactors=F,sep=",",skip=0,header = T)
par(mfrow=c(3,3))
P.fit=c()
p.CGC=c()
R.PPI=c()
p.PPI=c()
AllT.d=c()
for(kc in 1:33){
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,c(3,5)]),c(3,5)]
  RBP.gene=unique(RBP.target$RBP_sym)
  T.gene=unique(RBP.target$Tar_sym)
  RBP.degree=c()
  for(i in 1:length(RBP.gene)){
    print(i)
    xx=which(RBP.target$RBP_sym==RBP.gene[i])
    RBP.degree=rbind(RBP.degree,length(xx))
  }
  T.degree=c()
  for(i in 1:length(T.gene)){
    print(i)
    xx=which(RBP.target$Tar_sym==T.gene[i])
    T.degree=rbind(T.degree,length(xx))
  }
  RBP.degree=cbind(RBP.gene,RBP.degree)
  T.degree=cbind(T.gene,T.degree)
  hist(as.numeric(T.degree[,2]),50)
  T.PPI=merge(HI3,T.degree,by.x = "name",by.y="T.gene")
  AllT.d=rbind(AllT.d,T.PPI)
  R=cor(as.numeric(T.PPI$Degree), as.numeric(T.PPI$V2), method = "spearman")
  pp=cor.test(as.numeric(T.PPI$Degree), as.numeric(T.PPI$V2), method = "spearman")$p.value
  R.PPI=rbind(R.PPI,R)
  p.PPI=rbind(p.PPI,pp)
  CGC_d=merge(CGC,T.degree,by.x = "Gene.Symbol",by.y = "T.gene")
  other.cgc=setdiff(T.degree[,1],CGC$Gene.Symbol)
  other.cgc=t(t(other.cgc))
  other.cgc.d=merge(other.cgc,T.degree,by.x="V1",by.y="T.gene")
  ##boxplot(as.numeric(CGC_d$V2),as.numeric(other.cgc.d$V2))
  p1=wilcox.test(as.numeric(CGC_d$V2),as.numeric(other.cgc.d$V2),alternative = "greater")$p.value
  p.CGC=rbind(p.CGC,p1)
  fit.d=merge(core.fit,T.degree,by.x = "V1",by.y = "T.gene")
  other.fit=setdiff(T.degree[,1],core.fit)
  other.fit=t(t(other.fit))
  other.fit.d=merge(other.fit,T.degree,by.x="V1",by.y="T.gene")
 ## boxplot(as.numeric(fit.d$V2),as.numeric(other.fit.d$V2))
  p2=wilcox.test(as.numeric(fit.d$V2),as.numeric(other.fit.d$V2),alternative = "greater")$p.value
  P.fit=rbind(P.fit,p2)
}
kk=sort.int(-log10(p.PPI),index.return = T)
xx=c(1:33)
a=which(kk$x>1.3)
b=which(kk$x<1.3)
plot(xx,rep(1.3,33),type = "l",ylim = c(0,35))
points(a,kk$x[a],col="red")
points(b,kk$x[b],col="gray")

kk=sort.int(-log10(P.all),index.return = T)
xx=c(1:33)
a=which(kk$x>1.3)
b=which(kk$x<1.3)
points(a,kk$x[a],col="red",pch=22)
points(b,kk$x[b],col="gray",pch=22)
#############################negative results
##########################################################################
##################################################obtained the mutations in RBP binding sites
setwd("C:/projects/RBP/Pancancer/Targetmut")
hg38.score=read.csv("hg38freq_all.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
ORF.id=read.csv("ORFID.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
ORF.score=merge(ORF.id,hg38.score,by.x="ORF_ID",by.y="ORF_ID")
write.table(ORF.score,"ORF_score.txt",sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)

##############################################
setwd("C:/projects/RBP/Pancancer/Targetmut")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
RBP.target=read.csv("RBP_gene_binding.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
ORF.score=read.csv("ORF_score.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
Cacner.N=c()
P.all=c()
par(mfrow=c(3,3))
for(kc in 33:33){
  setwd("C:/projects/RBP/Pancancer/RBP_network")
  cancer.sig=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  xx=which(cancer.sig$FDR<1.0e-3)
  cancer.sig=cancer.sig[xx,]
  cancer.site=merge(cancer.sig,RBP.target,by.x=c("RBP_sym","Tar_entrez"),by.y=c("RBP_gene","Entrez"))
  Target_mut=c()
  ##########all binding sites in cancer kc
  Target.gene=cancer.site[!duplicated(cancer.site[,c(2,9,10)]),c(2,9,10)]
  xx=which(ORF.score$cancer==cancer[kc])
  Cancer.mut=ORF.score[xx,] #############all mutations in cancer kc
  for(i in 1:dim(Cancer.mut)[1]){
    print(i)
    mut.loc=as.numeric(Cancer.mut$base_.[i])
    tindex=which(Target.gene$Tar_entrez==Cancer.mut$ENTREZ_GENE_ID[i]&Target.gene$start<=mut.loc&Target.gene$end>=mut.loc)
    if(length(tindex)>0){
      KK=cbind(Target.gene[tindex,],Cancer.mut[i,])
      Target_mut=rbind(Target_mut,KK)
    }
  }
  Target.infor=merge(Target_mut,cancer.site,by.x = c("Tar_entrez","start","end"),
                     by.y = c("Tar_entrez","start","end"))
  setwd("C:/projects/RBP/Pancancer/Targetmut")
  write.table(Cancer.mut,paste(cancer[kc],"_cancer_allmut.txt",sep = ""),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
  write.table(Target.infor,paste(cancer[kc],"_Mut_infor.txt",sep = ""),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
  # cadd_score<-as.numeric(Target_mut$integrated_fitCons_score)
  # ks<-which(is.na(cadd_score)==FALSE)
  # cadd_score<-cadd_score[ks]
  # N=dim(Cancer.mut)[1]
  # n=dim(Target_mut)[1]
  # Cacner.N=rbind(Cacner.N,cbind(cancer[kc],N,n))
  # RRscore=c()
  # for(i in 1:100){
  #   xx=sample(c(1:N),n)
  #   RRs=Cancer.mut[xx,]
  #   r.cadd=RRs$integrated_fitCons_score
  #   ks<-which(is.na(r.cadd)==FALSE)
  #   r.cadd<-r.cadd[ks]
  #   RRscore=rbind(RRscore,r.cadd)
  # }
  # RRscore=matrix(RRscore,ncol = 1)
  # xx=as.numeric(cadd_score)
  # xx=xx[which(is.na(xx)==F)]
  # yy=as.numeric(RRscore)
  # yy=yy[which(is.na(yy)==F)]
  # pp=wilcox.test(xx,yy,alternative = "greater")$p.value
  # P.all=rbind(P.all,pp)
  # library(vioplot)
  # vioplot(xx,yy)
}
#############################calculated the number of mutation and genes in each cancer
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC")
Num.mut=c()
setwd("C:/projects/RBP/Pancancer/Targetmut")
for(kc in 1:32){
  print(kc)
  cancer.sig=read.csv(paste(cancer[kc],"_Mut_infor.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  all.cancer=read.csv(paste(cancer[kc],"_cancer_allmut.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  cancer.sig=cancer.sig[!duplicated(cancer.sig),]
  all.cancer=all.cancer[!duplicated(all.cancer),]
  mut.num=dim(cancer.sig)[1]
  mut.all=dim(all.cancer)[1]
  mut.g=length(unique(cancer.sig$Tar_entrez))
  all.g=length(unique(all.cancer$ENTREZ_GENE_ID))
  Num.mut=rbind(Num.mut,cbind(cancer[kc],mut.num,mut.all,mut.g,all.g))
}
X1=as.numeric(as.character(Num.mut[,2]))/as.numeric(as.character(Num.mut[,3]))
X2=as.numeric(as.character(Num.mut[,4]))/as.numeric(as.character(Num.mut[,5]))
X3=as.numeric(as.character(Num.mut[,2]))
X4=as.numeric(as.character(Num.mut[,4]))
x=c(1:32)
par(mfrow=c(3,1))
barplot(t(cbind(log2(X3),log2(X4))),beside=T, 
        col=c("aquamarine3","coral"),axes=FALSE)
axis(2, ylim=c(0,14),col="black",las=1)
par(new=TRUE)
plot(x,X1,pch=15,col="aquamarine3",type ="b",axes=FALSE)
par(new=TRUE)
plot(x,X2,pch=15,col="coral",type ="b",axes=FALSE)
axis(4, ylim=c(0,0.55), col="red",col.axis="red",las=1)
######################CHOL no mutations in target sites
#################The mutation location on the binding sites
setwd("C:/projects/RBP/Pancancer/Targetmut")
cancer<-c("KIRC","KIRP","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","UVM","DLBC")
Bid.loc=c()
for(kc in 1:28){
  print(kc)
  cancer.mut=read.csv(paste(cancer[kc],"_Mut_infor.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  cancer.mut=cancer.mut[!duplicated(cancer.mut),]
  T.loc=c()
  for(i in 1:dim(cancer.mut)[1]){
    print(i)
    mut.loc=cancer.mut$base_.[i]
    st=cancer.mut$start[i]
    ed=cancer.mut$end[i]
    leg=ed-st+1
    Ml=mut.loc-st+1
    T.loc=rbind(T.loc,Ml/leg)
  }
  cancer.mut=cbind(cancer.mut,T.loc)
  x1=sum(T.loc<0.33)
  x2=length(which(T.loc>=0.33&T.loc<=0.66))
  x3=length(which(T.loc>0.66))
  n=length(T.loc)
  Bid.loc=rbind(Bid.loc,cbind(x1/n,x2/n,x3/n))
}
barplot(t(Bid.loc))
boxplot(Bid.loc,boxwex=c(0.5,0.5,0.5),outline = FALSE)
wilcox.test(Bid.loc[,1],Bid.loc[,2],alternative = "less")
wilcox.test(Bid.loc[,2],Bid.loc[,3],alternative = "less")

################################store the files of for the number of RBPs regulated each gene
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
setwd("C:/projects/RBP/Pancancer/RBP_network")
for(kc in 1:33){
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,c(3,5)]),c(3,5)]
  RBP.gene=unique(RBP.target$RBP_sym)
  T.gene=unique(RBP.target$Tar_sym)
  RBP.degree=c()
  for(i in 1:length(RBP.gene)){
    print(i)
    xx=which(RBP.target$RBP_sym==RBP.gene[i])
    RBP.degree=rbind(RBP.degree,length(xx))
  }
  T.degree=c()
  for(i in 1:length(T.gene)){
    print(i)
    xx=which(RBP.target$Tar_sym==T.gene[i])
    T.degree=rbind(T.degree,length(xx))
  }
  RBP.degree=cbind(RBP.gene,RBP.degree)
  T.degree=cbind(T.gene,T.degree)
  write.table(RBP.degree,paste(cancer[kc],"_RBP_degree.txt",sep = ""),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
  write.table(T.degree,paste(cancer[kc],"_T_degree.txt",sep = ""),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
}
########################store the expression and var
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")

for(kc in 1:33){
  print(kc)
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep=""))
  library(SummarizedExperiment)
  LIHCmatrix<-assay(data,1,"FPKM")############gene exprssion
  ###################filter the gene expression profile
  N=dim(LIHCmatrix)[2]
  filter_line<-c()
  for (i in 1:dim(LIHCmatrix)[1]) {
    n=which(LIHCmatrix[i,]<0.1)
    if(length(n)>N*0.3){
      filter_line<-rbind(filter_line,i)
    }
  }
  LIHCmatrix<-LIHCmatrix[-filter_line,]
  LIHCclinal<-data@colData@listData$definition#############sample information tumor vs normal
  LIHCclinal=t(t(LIHCclinal))
  #setwd("C:/projects/RBP/Pancancer")
  #RBP<-read.csv("EnsemRBP.txt",stringsAsFactors=F,sep="\t",skip=0,header =F)
  LIHCmatrix_nor=log2(LIHCmatrix+0.05)
  M_exp=c()
  V_exp=c()
  for (i in 1:dim(LIHCmatrix_nor)[1]) {
    print(i)
    M_exp=rbind(M_exp,mean(LIHCmatrix_nor[i,]))
    V_exp=rbind(V_exp,var(LIHCmatrix_nor[i,]))
  }
  GeneMV=cbind(rownames(LIHCmatrix),M_exp,V_exp)
  write.table(GeneMV,paste(cancer[kc],"_Gene_MV.txt",sep = ""),sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
}


##########################mut target genes enrich in CGC or not
setwd("C:/projects/RBP/Pancancer/Targetmut")
cancer<-c("KIRC","KIRP","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","UVM","DLBC")
CGC=read.csv("Census_allSun Nov 13 17-37-13 2016.tsv",stringsAsFactors=F,sep="\t",skip=0,header = T)
fitness=read.csv("corefitness.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Protein.gene=read.csv("Protein_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
P.all=c()
for(kc in 1:28){
  print(kc)
  allmut=read.csv(paste(cancer[kc],"_cancer_allmut.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  cancer.mut=read.csv(paste(cancer[kc],"_Mut_infor.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  allmut=allmut[!duplicated(allmut),]
  cancer.mut=cancer.mut[!duplicated(cancer.mut),]
  mgene=unique(cancer.mut$Tar_sym)
  O.CGC=intersect(mgene,CGC$Gene.Symbol)
  O.fit=intersect(mgene,fitness$V1)
  all.mutg=unique(allmut$ENTREZ_GENE_ID)
  all.mutg=t(t(all.mutg))
  all.s=merge(all.mutg,Protein.gene,by.x="V1",by.y="Entrez")
  all.s=all.s$Symbol
  R.CGC=c()
  R.fit=c()
  for(i in 1:1000){
    x1=sample(c(1:length(all.s)),length(mgene))
    r.cgc=intersect(all.s[x1],CGC$Gene.Symbol)
    r.fit=intersect(all.s[x1],fitness$V1)
    R.CGC=rbind(R.CGC,length(r.cgc))
    R.fit=rbind(R.fit,length(r.fit))
  }
  p1=sum(R.CGC>length(O.CGC))/1000
  p2=sum(R.fit>length(O.fit))/1000
  XX1=allmut$integrated_fitCons_score
  XX2=cancer.mut$integrated_fitCons_score
  XX1=as.numeric(as.character(XX1))
  XX2=as.numeric(as.character(XX2))
  ks1<-which(is.na(XX1)==FALSE)
  ks2<-which(is.na(XX2)==FALSE)
  XX1=XX1[ks1]
  XX2=XX2[ks2]
  p6=wilcox.test(XX1,XX2,alternative = "less")$p.value
  #######degree
  T.degree=read.csv(paste("C:/projects/RBP/Pancancer/RBP_network/",cancer[kc],"_T_degree.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  mgene=t(t(mgene))
  m.degree=merge(T.degree,mgene,by.x = "T.gene",by.y = "V1")
  other.gene=setdiff(T.degree$T.gene,mgene[,1])
  other.gene=t(t(other.gene))
  other.degree=merge(T.degree,other.gene,by.x = "T.gene",by.y = "V1")
  p3=wilcox.test(as.numeric(as.character(m.degree$X)),as.numeric(as.character(other.degree$X)),alternative = "greater")$p.value
  #####expression
  GeneMV=read.csv(paste(cancer[kc],"_Gene_MV.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = F)
  Sy.gene=merge(GeneMV,Protein.gene,by.x="V1",by.y="Ensg")
  m.mean=merge(Sy.gene,mgene,by.x="Symbol",by.y="V1")
  o.mean=merge(Sy.gene,other.gene,by.x="Symbol",by.y="V1")
  p4=wilcox.test(as.numeric(as.character(m.mean$V2)),as.numeric(as.character(o.mean$V2)),alternative = "greater")$p.value
  p5=wilcox.test(as.numeric(as.character(m.mean$V3)),as.numeric(as.character(o.mean$V3)),alternative = "less")$p.value
  P.all=rbind(P.all,cbind(cancer[kc],length(mgene),p1,p2,p3,p4,p5,p6))
}

XX=P.all[,3:8]
Xa=(as.numeric(XX)<0.05)
Xa=matrix(Xa,nrow = 28)
ya=apply(Xa,1,sum)
yb=apply(Xa,2,sum)
barplot(ya)
barplot(yb)
#################################Hallmark enrichment
setwd("C:/projects/RBP/Pancancer/Targetmut")
cancer<-c("KIRC","KIRP","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","UVM","DLBC")
P.all=c()
for(kc in 1:28){
  print(kc)
  cancer.mut=read.csv(paste(cancer[kc],"_Mut_infor.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  cancer.mut=cancer.mut[!duplicated(cancer.mut),]
  mut_gene=unique(cancer.mut$Tar_sym)
  Biofunction <- file("go2hallmark1.GMT", "r")
  i<-0
  N <-19022
  line=readLines(Biofunction,n=1)
  Result1 <- c()
  while( length(line) != 0 ) {
    i <- i+1
    msig.go <- strsplit(line,"\t")
    msig.go <- msig.go[[1]]
    m <- length(msig.go)
    goterm <- msig.go[1]
    msig.go <- msig.go[3:m]
    n1 <- length(mut_gene)
    intergene1 <- intersect(msig.go,mut_gene)
    k1 <- length(intergene1)
    print(k1)
    if (k1>2){
      p.value1 <- phyper(k1-1,m,N-m,n1,lower.tail=FALSE, log.p = FALSE)
      go <- c(goterm,k1,m,n1,p.value1)
      Result1 <- rbind(Result1,go)
    } else {
      go <- c(goterm,k1,m,n1,1)
      Result1 <- rbind(Result1,go)
    }
    line=readLines(Biofunction,n=1);
  }
  close(Biofunction)
  P.all=rbind(P.all,cbind(kc,c(1:34),Result1[,5]))
}
p.ALL1=P.all
x1=which(P.all[,3]<=0.01)
P.all[x1,3]=0.01
x2=which(P.all[,3]<=0.05&P.all[,3]>0.01)
P.all[x2,3]=0.05
x3=which(P.all[,3]<=0.1&P.all[,3]>0.05)
P.all[x3,3]=0.1
x4=which(P.all[,3]>0.1)
P.all[x4,3]=NA
dfx = data.frame(ev1=as.numeric(P.all[,1]), ev2=as.numeric(P.all[,2]), ev3=-log2(as.numeric(P.all[,3])))
k=with(dfx, symbols(x=ev1, y=ev2, circles=ev3, inches=1/6,
                    ann=F, bg="steelblue2", fg=NULL,ylim =c(0,35)))
ylabel <- seq(0,35, by = 1)
axis(2, at = ylabel)
xlabel <- seq(0,30, by = 1)
axis(1, at = xlabel)
grid(nx = NULL, ny =NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)


#############################LIHC CO-regulation
setwd("C:/projects/RBP/newidea")
RBPnet<-read.csv("Net.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
RBP=unique(RBPnet$V2)
RBP_jac=c()
n=length(RBP)
for(i in 1:(n-1)){
  for(j in (i+1):n){
    x1=which(RBPnet$V2==RBP[i])
    x2=which(RBPnet$V2==RBP[j])
    G1=unique(RBPnet$V4[x1])
    G2=unique(RBPnet$V4[x2])
    n1=length(intersect(G1,G2))
    n2=length(union(G1,G2))
    RBP_jac=rbind(RBP_jac,cbind(RBP[i],RBP[j],n1/n2,n1))
  }
}

HI3=read.csv("C:/projects/RBP/Pancancer/RBP_network/SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)
Lab=c()
for(i in 1:dim(RBP_jac)[1]){
  x1=which(HI3$V1==RBP_jac[i,1]&HI3$V2==RBP_jac[i,2])
  x2=which(HI3$V2==RBP_jac[i,1]&HI3$V1==RBP_jac[i,2])
  x=union(x1,x2)
  if(length(x)>0){
    Lab=rbind(Lab,1)
  } else {
    Lab=rbind(Lab,0)
  }
}
RBP_jac=cbind(RBP_jac,Lab)
x1=which(RBP_jac[,5]==1)
x2=which(RBP_jac[,5]==0)
Y1=as.numeric(as.character(RBP_jac[x1,4]))
Y2=as.numeric(as.character(RBP_jac[x2,4]))
par(mfrow=c(3,3))
library(vioplot)
vioplot(Y1,Y2)
wilcox.test(Y1,Y2,alternative = "greater")
##############RBPnetwork OVERLAP
HCC=RBPnet[!duplicated(RBPnet[,c(2,4)]),c(2,4)]
setwd("C:/projects/RBP/Pancancer/RBP_network")
RBP.target=read.csv("LIHC _RBP_corr_sig.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
XX=which(as.numeric(RBP.target$FDR)<1e-3)
RBP.target=RBP.target[XX,]
RBP.target=RBP.target[!duplicated(RBP.target[,c(3,5)]),c(3,5)]
O.rbp=intersect(RBP.target$RBP_sym,HCC$V2)  
O.gene=intersect(RBP.target$Tar_sym,HCC$V4)
O.rbp=t(t(O.rbp))
O.gene=t(t(O.gene))

HCC.1=merge(HCC,O.rbp,by.x="V2",by.y="V1")
HCC.2=merge(HCC.1,O.gene,by.x="V4",by.y="V1")
LIHC.1=merge(RBP.target,O.rbp,by.x="RBP_sym",by.y="V1")
LIHC.2=merge(LIHC.1,O.gene,by.x="Tar_sym",by.y="V1")
o.Reg=merge(LIHC.2,HCC.2,by.x=c("RBP_sym","Tar_sym"),by.y=c("V2","V4"))
p=phyper(1142-1,6219,105948-6219,19209,lower.tail=FALSE, log.p = FALSE)
  
######################Mutation network rewiring
setwd("C:/projects/RBP/Pancancer/Targetmut")
cancer<-c("KIRC","KIRP","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","UVM","DLBC")
Three.sam=function(fusion_sample){
  fusionthree=c()
  for (i in 1:length(fusion_sample)) {
    aa=strsplit(fusion_sample[i],"-")
    bb=paste(aa[[1]][1],aa[[1]][2],aa[[1]][3],sep="-")
    fusionthree=rbind(fusionthree,bb)
  }
  return(fusionthree)
}

for(kc in 1:28){
  print(kc)
  Mut_network=c()
  cancer.mut=read.csv(paste(cancer[kc],"_Mut_infor.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  cancer.mut=cancer.mut[!duplicated(cancer.mut),]
  Umut=cbind(cancer.mut$Tar_sym,cancer.mut$Tar_ensg,cancer.mut$chromosome,cancer.mut$start,
             cancer.mut$end,cancer.mut$HGVSc,cancer.mut$HGVSp_short,cancer.mut$RBP_sym,
             cancer.mut$RBP_ensg,cancer.mut$R)
  Umut=Umut[!duplicated(Umut),]
  Allmut=read.csv(paste("C:/projects/TCGAmutation/",cancer[kc],".maf",sep=""),stringsAsFactors=F,sep="\t",skip=1,header = T)
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep=""))
  library(SummarizedExperiment)
  LIHCmatrix<-assay(data,1,"FPKM")############gene exprssion
  LIHCclinal<-data@colData@listData$definition#############sample information tumor vs normal
  LIHCclinal=t(t(LIHCclinal))
  Allgene=rownames(LIHCmatrix)
  cancer.index=which(LIHCclinal!="Solid Tissue Normal")
  All.sample=colnames(LIHCmatrix)
  LIHCmatrix=LIHCmatrix[,cancer.index]
  All.sample=All.sample[cancer.index]
  All.sample=Three.sam(All.sample)
  for(i in 1:dim(Umut)[1]){
    #print(i)
    xx=which(Allmut$Hugo_Symbol==Umut[i,1]&Allmut$Chromosome==Umut[i,3]&
             Allmut$HGVSc==Umut[i,6]&Allmut$HGVSp_Short==Umut[i,7])
    if(length(xx)>0){
     mutsam=unique(Allmut$Tumor_Sample_Barcode[xx])
     mutsam=Three.sam(mutsam)
     mutindex=c()
     for(j in 1:length(mutsam)){
       xa=which(All.sample==mutsam[j])
       if(length(xa)>0){
         mutindex=rbind(mutindex,xa)
       }
     }
     print(length(mutindex))
     norindex=setdiff(c(1:length(All.sample)),mutindex)
     ge.index=which(Allgene==Umut[i,2])
    if(length(ge.index)>0&length(mutindex)>0){
      Xmut=as.numeric(as.character(LIHCmatrix[ge.index,mutindex]))
      Xnor=as.numeric(as.character(LIHCmatrix[ge.index,norindex]))
      if(mean(Xnor)>0&mean(Xmut)>0){
        FC=mean(Xnor)/mean(Xmut)
        KW=cbind(t(Umut[i,]),FC)
        Mut_network=rbind(Mut_network,KW)
      }
    }
    }
  }
  colnames(Mut_network)=c("Tar_sym","Tar_ensg","chr","start","end","mut_g","mut_p","RBP_sym","RBP_ensg",
                          "R","FC")
  write.table(Mut_network,paste("C:/projects/RBP/Pancancer/mutnet/",cancer[kc],"mut_network.txt"),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
}
  
###################extract network
setwd("C:/projects/RBP/Pancancer/mutnet")
cancer<-c("KIRC","KIRP","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","UVM","DLBC")
Network=c()
Number.C=c()
for(kc in 1:28){
  print(kc)
  mutnet=read.csv(paste(cancer[kc],"_mut_network.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  net.u=mutnet[!duplicated(mutnet[,c(8,1)]),c(8,1)]
  Number.C=rbind(Number.C,dim(net.u)[1])
  x1=which(mutnet$R>0&mutnet$FC<0.5)
  if(length(x1)>0){
    Dac=mutnet[x1,c(8,1)]
    Dac=Dac[!duplicated(Dac),]
    Dac=cbind(Dac,"Dac",cancer[kc])
    colnames(Dac)=c("RBP","Tar","Type","Cancer")
    Network=rbind(Network,Dac)
  }
  if(length(x2)>0){
    x2=which(mutnet$R<0&mutnet$FC>2)
    Din=mutnet[x2,c(8,1)]
    Din=Din[!duplicated(Din),]
    Din=cbind(Din,"Din",cancer[kc])
    colnames(Din)=c("RBP","Tar","Type","Cancer")
    Network=rbind(Network,Din)
  }
}

cancer.DI=c()
for(kc in 1:28){
  x1=which(Network$Cancer==cancer[kc]&Network$Type=="Dac")
  x2=which(Network$Cancer==cancer[kc]&Network$Type=="Din")
  cancer.DI=rbind(cancer.DI,cbind(length(x1),length(x2)))
}  
Cancer.rew=cbind(Number.C,cancer.DI)

X1=as.numeric(as.character(Cancer.rew[,2]))/as.numeric(as.character(Cancer.rew[,1]))
X2=as.numeric(as.character(Cancer.rew[,3]))/as.numeric(as.character(Cancer.rew[,1]))
X3=as.numeric(as.character(Cancer.rew[,2]))
X4=as.numeric(as.character(Cancer.rew[,3]))
x=c(1:28)
par(mfrow=c(3,1))
barplot(t(cbind(X3,X4)),beside=T, 
        col=c("aquamarine3","coral"),axes=FALSE)
axis(2, col="black",las=1)
par(new=TRUE)
plot(x,X1,pch=15,col="aquamarine3",type ="b",axes=FALSE)
par(new=TRUE)
plot(x,X2,pch=15,col="coral",type ="b",axes=FALSE)
axis(4, ylim=c(0,0.25), col="red",col.axis="red",las=1)  
  
#################rewiring network hub target degree
Net.U=Network[!duplicated(Network[,1:2]),1:2]
RT=unique(Net.U$Tar)
DD=c()
for(i in 1:length(RT)){
  xx=which(Net.U$Tar==RT[i])
  rj=unique(Net.U$RBP[xx])
  DD=rbind(DD,length(rj))
}
RT_d=cbind(RT,DD)
write.table(RT_d,"RT.rnk",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)

Net.U=Network[!duplicated(Network[,1:3]),1:3]
write.table(Net.U,"Rewiringnet.txt",row.names = FALSE, col.names = F,quote = FALSE,sep = "\t")
hist(as.numeric(RT_d[,2]),20)
xx=which(as.numeric(RT_d[,2])>5)
hub=RT_d[xx,]
Hubnet=merge(hub,Net.U,by.x="RT",by.y="Tar")
write.table(Hubnet,"Rewiring_hub_net.txt",row.names = FALSE, col.names = T,quote = FALSE,sep = "\t")

n.c=c()
for(i in 1:28){
  aa=which(Network$Cancer==cancer[i])
  LIHC_net=Network[aa,]
  n.c=rbind(n.c,dim(LIHC_net)[1])
  write.table(LIHC_net,paste(cancer[i],"_R_net.txt",sep=""),row.names = FALSE, col.names = T,quote = FALSE,sep = "\t")
}
hallmark=read.csv("hallmark.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)

adnet=merge(hallmark,Net.U,by.x = "V2",by.y = "Tar")

write.table(adnet,"Rewiring_ad_net.txt",row.names = FALSE, col.names = T,quote = FALSE,sep = "\t")

RR=unique(Net.U$RBP)
RR=cbind(RR,"RBP")
TT=unique(Net.U$Tar)
TT=cbind(TT,"target")
node=rbind(TT,RR)

write.table(node,"Rewiring_node.txt",row.names = FALSE, col.names = F,quote = FALSE,sep = "\t")
Jad=c()
Rec=c()
for(i in 1:27){
  for(j in (i+1):28){
    x1=which(Network$Cancer==cancer[i])
    x2=which(Network$Cancer==cancer[j])
    C1=Network[x1,]
    C2=Network[x2,]
    OO=merge(C1,C2,by.x=c("RBP","Tar"),by.y=c("RBP","Tar"))
    KK=cbind(cancer[i],cancer[j],dim(C1)[1],dim(C2)[1],dim(OO)[1])
    Rec=rbind(Rec,OO)
    Jad=rbind(Jad,KK)
  }
}
Rec=Rec[!duplicated(Rec[,1:2]),1:2]
write.table(Rec,"Rewiring_Rec.txt",row.names = FALSE, col.names = T,quote = FALSE,sep = "\t")

jindex=c()
for(i in 1:dim(Jad)[1]){
  xx=as.numeric(Jad[i,5])/min(as.numeric(Jad[i,3]),as.numeric(Jad[i,4]))
  jindex=rbind(jindex,xx)
}
library(ggplot2)
jindex=data.frame(V1=jindex[,1])
ggplot(jindex,aes(V1))+geom_density()


#####################additional file----tables
####################RBP mutation table
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
Sample_ratio=c()
N.cancer=c()
for(kc in 1:33){
  print(kc)
  setwd("C:/projects/RBP/newidea")
  RBP<-read.csv("RBP_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header =T)
  setwd("C:/projects/TCGAmutation")
  mutation<-read.csv(paste(cancer[kc],".maf",sep = ""),stringsAsFactors=F,sep="\t",skip=1,header = TRUE)
  N=length(unique(mutation$Tumor_Sample_Barcode))
  R=c()
  for(i in 1:dim(RBP)[1]){
    xx=which(mutation$Hugo_Symbol==RBP$HGNC.symbol[i])
    yy=unique(mutation$Tumor_Sample_Barcode[xx])
    R=rbind(R,length(yy))
  }
  N.cancer=rbind(N.cancer,N)
  Sample_ratio=cbind(Sample_ratio,R)
}
RBP_num=cbind(RBP,Sample_ratio)
write.table(RBP_num,"RBP_pan_cancer.txt",row.names = FALSE, col.names = T,quote = FALSE,sep = "\t")

RR=c()
for(i in 1:33){
  xx=rank(-as.numeric(RBP_num[,(i+4)]))
  xx=xx/max(xx)
  RR=cbind(RR,xx)
}
KW=apply(RR,1,mean)
KW=cbind(RBP_num$HGNC.symbol,KW)
write.table(KW,"RBP_MUT_pan.rnk",row.names = FALSE, col.names = F,quote = FALSE,sep = "\t")

######################The number of samples with expression

cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
Sam.n=c()
for(kc in 1:33){
  print(kc)
  load(paste("C:/projects/AS/geneexpression/",cancer[kc],"Expression.rda",sep=""))
  library(SummarizedExperiment)
  LIHCmatrix<-assay(data,1,"FPKM")############gene exprssion
  ###################filter the gene expression profile
  LIHCclinal<-data@colData@listData$definition#############sample information tumor vs normal
  N=length(LIHCclinal)
  XX=which(LIHCclinal=="Solid Tissue Normal")
  kk=cbind(length(XX),N-length(XX))
  Sam.n=rbind(Sam.n,kk)
}

################################RBP-gene network table
setwd("C:/projects/RBP/Pancancer/RBP_network")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
for(kc in 1:33){
  print(kc)
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,c(3,5:8)]),c(3,5:8)]
  write.table(RBP.target,paste("C:/projects/RBP/Pancancer/Net_table/",cancer[kc],"_network.txt",sep = ""),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
}




setwd("C:/projects/RBP/Pancancer/RBP_network")
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCS","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","ACC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","MESO","UVM","DLBC","CHOL")
HI3=read.csv("SCIN.txt",stringsAsFactors=F,sep="\t",skip=0,header = F)

par(mfrow=c(3,3))

P.all=c()
for(kc in 1:33){
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,3:4]),3:4]
  RBP.gene=unique(RBP.target$RBP_sym)
  RR_g=c()
  n=length(RBP.gene)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      print(c(i,j))
      xx1=which(RBP.target$RBP_sym==RBP.gene[i])
      xx2=which(RBP.target$RBP_sym==RBP.gene[j])
      TT1=unique(RBP.target$Tar_entrez[xx1])
      TT2=unique(RBP.target$Tar_entrez[xx2])
      OT=intersect(TT1,TT2)
      UT=union(TT1,TT2)
      jad=length(OT)/length(UT)
      R=cbind(RBP.gene[i],RBP.gene[j],jad)
      RR_g=rbind(RR_g,R)
    }
  }
  plot(density(as.numeric(RR_g[,3])),col="red",ylim=c(0,40))
  N=dim(RBP.target)[1]
  xx=sample(c(1:N),N)
  RBP.target.r=cbind(RBP.target$RBP_sym,RBP.target$Tar_entrez[xx])
  RR_g_r=c()
  n=length(RBP.gene)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      print(c(i,j))
      xx1=which(RBP.target.r[,1]==RBP.gene[i])
      xx2=which(RBP.target.r[,1]==RBP.gene[j])
      TT1=unique(RBP.target.r[xx1,2])
      TT2=unique(RBP.target.r[xx2,2])
      OT=intersect(TT1,TT2)
      UT=union(TT1,TT2)
      jad=length(OT)/length(UT)
      R=cbind(RBP.gene[i],RBP.gene[j],jad)
      RR_g_r=rbind(RR_g_r,R)
    }
  }
  lines(density(as.numeric(RR_g_r[,3])),col="blue")
  p.value=ks.test(as.numeric(RR_g[,3]),as.numeric(RR_g_r[,3]),alternative = "greater")$p.value
  P.all=rbind(P.all,p.value)
}

################cancer similarity
setwd("C:/projects/RBP/Pancancer/mutnet")
cancer<-c("KIRC","KIRP","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC","PCPG","SARC","LAML","PAAD","ESCA","TGCT","THYM","UVM","DLBC")
Jad=c()
for(i in 1:27){
  for(j in (i+1):28){
    print(c(i,j))
    setwd("C:/projects/RBP/Pancancer/mutnet")
    mutnet1=read.csv(paste(cancer[i],"_mut_network.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
    net.u1=mutnet1[!duplicated(mutnet1[,c(8,1)]),c(8,1)]
    mutnet2=read.csv(paste(cancer[j],"_mut_network.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
    net.u2=mutnet2[!duplicated(mutnet2[,c(8,1)]),c(8,1)]
    O.net=merge(net.u1,net.u2,by.x = c("RBP_sym","Tar_sym"),by.y = c("RBP_sym","Tar_sym"))
    n1=dim(net.u1)[1]
    n2=dim(net.u2)[1]
    n.o=dim(O.net)[1]
    setwd("C:/projects/RBP/Pancancer/Net_table")
    net1=read.csv(paste(cancer[i],"_network.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
    net1=net1[!duplicated(net1[,1:2]),1:2]
    net2=read.csv(paste(cancer[j],"_network.txt",sep = ""),stringsAsFactors=F,sep="\t",skip=0,header = T)
    net2=net2[!duplicated(net2[,1:2]),1:2]
    RRO=merge(net1,net2,by.x=c("RBP_sym","Tar_sym"),by.y=c("RBP_sym","Tar_sym"))
    m.o=dim(RRO)[1]
    N1=dim(net1)[1]
    N2=dim(net2)[1]
    x1=sample(c(1:N1),n1)
    x2=sample(c(1:N2),n2)
    R1=net1[x1,]
    R2=net2[x2,]
    RR=merge(R1,R2,by.x=c("RBP_sym","Tar_sym"),by.y=c("RBP_sym","Tar_sym"))
    n.r=dim(RR)[1]
    Jad=rbind(Jad,cbind(cancer[i],cancer[j],n1,n2,n.o,n.r,N1,N2,m.o))
  }
} 
Jadd_n=c()
for(i in 1:dim(Jad)[1]){
  J1=as.numeric(Jad[i,5])/min(as.numeric(Jad[i,3]),as.numeric(Jad[i,4]))
  J2=as.numeric(Jad[i,9])/min(as.numeric(Jad[i,7]),as.numeric(Jad[i,8]))
  Jadd_n=rbind(Jadd_n,cbind(J1,J2))
}
X1=as.numeric(Jadd_n[,1]) #####real
Y1=as.numeric(Jadd_n[,2]) ########random
hist(Y1,100,probability = T,col = "gray")
lines(density(Y1),col="gray")
hist(X1,50,probability = T,col = "skyblue",add=T)
lines(density(X1),col="skyblue")

boxplot(Y1,X1)
ks.test(X1,Y1,alternative = "greater")
X1=cbind(X1,"real")
Y1=cbind(Y1,"rand")
Z=rbind(X1,Y1)
jindex=data.frame(V1=Z[,2],V2=as.numeric(Z[,1]))
library(ggplot2)
ggplot(jindex,aes(V2, fill = V1, colour = V1))+geom_density()


#####################################lincRNA coexpressed genes
load(paste("C:/projects/AS/geneexpression/LIHCExpression.rda",sep=""))
library(SummarizedExperiment)
LIHCmatrix<-assay(data,1,"FPKM")############gene exprssion
P.gene=read.csv("Protein_gene.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
allgene=rownames(LIHCmatrix)
xx=which(allgene=="ENSG00000212694")
Lincexp=LIHCmatrix[xx,]
RR=c()
P=c()
for(i in 1:dim(P.gene)[1]){
  x1=which(allgene==P.gene$Ensg[i])
  print(i)
  if(length(x1)>0){
    g.ep=LIHCmatrix[x1,]
    R=cor(as.numeric(Lincexp), as.numeric(g.ep), method = "spearman")
    pp=cor.test(as.numeric(Lincexp), as.numeric(g.ep), method = "spearman")$p.value
    RR=rbind(RR,cbind(P.gene$Symbol[i],R,pp))
  }
}
x1=which(is.na(RR[,2])==T)
RR=RR[-x1,1:3]
write.table(RR,"linc_coexp.txt",sep = "\t",row.names = FALSE, col.names = F,quote = FALSE)
########################compute the degree of genes in TF-gene network
cancer<-c("KIRC","KIRP","KICH","LGG","GBM","BRCA","LUSC","LUAD","READ","COAD","UCEC","OV","HNSC","THCA","PRAD","STAD","SKCM","BLCA",
          "LIHC","CESC")
All.TF.d=c()
R.all=c()
P.all=c()
for(kc in 1:20){
  setwd(paste("C:/projects/RBP/Pancancer/TF_gene/",cancer[kc],sep = ""))
  TF.reg=read.csv("result.txt",stringsAsFactors=F,sep="\t",skip=0,header = T)
  xx=which(TF.reg$FDR<0.01)
  TF.reg=TF.reg[xx,]
  T.gene=unique(TF.reg$gene_name)
  T.degree=c()
  for(i in 1:length(T.gene)){
    print(i)
    x1=which(TF.reg$gene_name==T.gene[i])
    tf=unique(TF.reg$TF_name[x1])
    T.degree=rbind(T.degree,length(tf))
  }
  T.degree=cbind(T.gene,T.degree)
  colnames(T.degree)=c("gene","degree")
  #write.table(T.degree,paste(cancer[kc],"T_degree.txt",sep=""),sep = "\t",row.names = FALSE, col.names = T,quote = FALSE)
  setwd("C:/projects/RBP/Pancancer/RBP_network")
  RBP.target=read.csv(paste(cancer[kc]," _RBP_corr_sig.txt",sep=""),stringsAsFactors=F,sep="\t",skip=0,header = T)
  XX=which(as.numeric(RBP.target$FDR)<1e-3)
  RBP.target=RBP.target[XX,]
  RBP.target=RBP.target[!duplicated(RBP.target[,c(3,5)]),c(3,5)]
  RBP.gene=unique(RBP.target$RBP_sym)
  Tb.gene=unique(RBP.target$Tar_sym)
  Tb.degree=c()
  for(i in 1:length(Tb.gene)){
    print(i)
    xx=which(RBP.target$Tar_sym==Tb.gene[i])
    Tb.degree=rbind(Tb.degree,length(xx))
  }
  Tb.degree=cbind(Tb.gene,Tb.degree)
  KK.degree=merge(T.degree,Tb.degree,by.x="gene",by.y="Tb.gene")
  R=cor(as.numeric(KK.degree$degree), as.numeric(KK.degree$V2), method = "spearman")
  pp=cor.test(as.numeric(KK.degree$degree), as.numeric(KK.degree$V2), method = "spearman")$p.value
  R.all=rbind(R.all,R)
  P.all=rbind(P.all,pp)
  All.TF.d=rbind(All.TF.d,KK.degree)########store all the degree of genes TF/RBP
}
plot(as.numeric(KK.degree$degree), as.numeric(KK.degree$V2),pch=19)


