options(stringsAsFactors=F)
options(max.print=1000)
library(spls)
library(limma)
tab=read.table("AffinityMatrix_diffBase_DBA_SCORE_READS_1_pVal_FALSE_0.05_FALSE_DESeq2.txt",sep="\t",header=T)
fpMat=tab[,c(1:4)]
d=fpMat[,3]-fpMat[,2]
m=floor(fpMat[,2]+(d)/2)
fpMat[,2]=m-300
fpMat[,3]=m+300
fpMat[,1]=paste0("chr",fpMat[,1])
cMat=as.matrix(tab[,5:ncol(tab)])
cMat=t(t(cMat)/colSums(cMat))*1000000
cMat=cMat/d*1000#RPKM normalized
#remove outliers
r=apply(cMat,2,quantile)
dRange=(r[4,]-r[2,])*2.5
outThres=r[4,]+dRange
idOut=rowSums(cMat>outThres)>0
cMat=cMat[!idOut,]
write.table(fpMat[!idOut,],"hg19.1.ATAC_v6_600bp.region",sep="\t",quote=F,row.names=F,col.names=F)
eMat3=normalizeBetweenArrays(cMat,method="quantile")
eMat3=cbind(rowMeans(eMat3[,c(2,3)]),rowMeans(eMat3[,c(1,4)]))
colnames(eMat3)=c("hESC","hDE")
resFrame=data.frame(fpMat[!idOut,c(1:4)],eMat3)
write.table(resFrame,"hg19_1_ATAC_v6_600bp_TERA_ATAC_norm_full.txt",sep="\t",row.names=F,quote=F)
#finished preprocessing
#####################################
#load motif matrix
tab=read.table("hg19_1_ATAC_v6_600bp_fimo_out_LogPvalRegion_matrix_top_TFmatrix.txt",sep="\t",header=T)
Xmat=as.matrix(tab[,c(5:ncol(tab))])
#get normalized Y matrix
Ymat=as.matrix(resFrame[,c(5:ncol(resFrame))])
Ymat=log2(Ymat+1)
#cap at max 5
Ymat[Ymat>5]=5
#filter out only non-zero entry motifs
Xmat=Xmat[,colSums(Xmat!=0)>50]
set.seed(1)
cv <- cv.spls( Xmat, Ymat, eta = seq(0.1,0.9,0.1),kappa=0.5, K = c(2:20),scale.y = T,scale.x=T )
dev.off()
#############
#process results
print(paste("K opt",cv$K.opt))#K =16
print(paste("eta opt",cv$eta.opt))# eta 0.1
#get final prediction
res=spls( Xmat, Ymat, K=cv$K.opt,eta=cv$eta.opt, kappa=0.5, select="pls2", fit="simpls",   scale.x=T, scale.y=T)
sR=res$betahat
write.table(sR,paste("TFA",cv$K.opt,eta=cv$eta.opt,"raw_final.txt",sep="_"),sep="\t",quote=F)
delta=sR[,2]-sR[,1]
names(delta)=rownames(sR)
gNam=rownames(sR)[order(delta,decreasing=T)]
#associate with expression data
expDat=read.table("EndoDiff_allGenes_FPKM_repCollapsed.txt",sep="\t",header=T)
#filter for minimal TF expression
expDat=expDat[rowSums(expDat[,c(3,4)]>=5)>0,]
gNam=gNam[is.element(gNam,expDat[,2])]#get TF with motifs that are expressed in at least one condition
#export selected
eDat=as.matrix(expDat[,c(3,4)])
rownames(eDat)=expDat[,2]
resFrame=data.frame(geneName=gNam,sR[gNam,],delta[gNam],eDat[gNam,],deltaExp=eDat[gNam,2]-eDat[gNam,1])
write.table(resFrame,"AnnotatedTFAP_resultFrame.txt",sep="\t",quote=F,row.names=F)
selCand=read.table("SelectedCandidates.txt",sep="\t",header=T)
resFrame=data.frame(geneName=gNam,sR[gNam,],delta[gNam],selected=as.numeric(is.element(gNam,selCand[,1])),eDat[gNam,],deltaExp=eDat[gNam,2]-eDat[gNam,1])
write.table(resFrame[resFrame[,"selected"]==1,],"AnnotatedTFAP_resultFrame_selectedCandidates.txt",sep="\t",quote=F,row.names=F)
colGradient = list()
colGradient[["tera"]] = colorRampPalette(c("blue","white","red"),bias=1,space="Lab")
library(fields)
library(gplots)
sel=resFrame[,"selected"]==1
pdf(file="atacTFAP.pdf", width=12, height=10)
par(mar=c(4,4,4,4),mfcol=c(2,1))
val=as.matrix(resFrame[sel,5])
m=max(abs(val))
m=0.05
val[val>m]=m
val[val<(-m)]=m
image(seq(1,nrow(val)),seq(1,ncol(val)),as.matrix(val),col=colGradient[["tera"]](20),zlim=c(-m,m),xaxt="n",yaxt="n",bty="or",ylab="",xlab="",main="All combined")
axis(1,at=seq(1,nrow(val)),labels=resFrame[sel,1],las=2,cex.axis=0.6)
image.plot(seq(1,nrow(val)),seq(1,ncol(val)),as.matrix(val),col=colGradient[["tera"]](20),zlim=c(-m,m),xaxt="n",yaxt="n",bty="or",ylab="",xlab="",legend.only=T)
val=as.matrix(log2(resFrame[sel,8]+1)-log2(resFrame[sel,7]+1))
par(mar=c(12,4,4,4),cex=1)
m=max(abs(val))
m=5
val[val>m]=m
val[val<(-m)]=m
image(seq(1,nrow(val)),seq(1,ncol(val)),as.matrix(val),col=greenred(75),zlim=c(-m,m),xaxt="n",yaxt="n",bty="or",ylab="",xlab="",main="All combined")
axis(1,at=seq(1,nrow(val)),labels=resFrame[sel,1],las=2,cex.axis=0.6)
image.plot(seq(1,nrow(val)),seq(1,ncol(val)),as.matrix(val),col=greenred(75),zlim=c(-m,m),xaxt="n",yaxt="n",bty="or",ylab="",xlab="",legend.only=T)
dev.off()
#######################

