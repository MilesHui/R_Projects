library(GEOquery)
library(oligo)
library(limma)
library(VIM)  # For handling missing values
library(mice) # For handling missing values
library(YuGene)
library(metaMA)

############################################## RNA ##################3####################################
################################# GSE74629: Illumina Platform#############################################
 studyID = "GSE74629"
 destdir = "C:/Users/think/Desktop/Gene/raw data/RNA new/GSE74629"
 setwd(destdir)
 eSet_1 <- getGEO(studyID, destdir = ".")
 eSet_1 <- eSet_1[[1]]  # This is eSet after author preprocess
dim(eSet_1)
# Get phenoData
 pdata = pData(eSet_1)
 write.csv(pdata,paste(studyID,"_matadata.csv")) 

# Construct grouping information 
 group_list_1 =subset(pdata,select=geo_accession) # Sample name
 group_list_1$condition = c(rep(c("T"),36),rep(c("N"),each=14))  # T:disease, N:normal grouping information
 group_list_1¡¡

# start with the expression matrix
# Manual download of non standardized raw data GSE74629_non-normalized.txt
 Elist_1 <- read.ilmn(files="GSE74629_non-normalized.txt",expr="SAMPLE",probeid="ID_REF",other.columns="Detection.Pval")

# See missing data
 md.pattern(Elist_1$E)       # No missing data

# background regulation, nec function, use normal-exponential convolution model
 Elist_1_bg <- nec(Elist_1,detection.p="Detection.Pval") #Elist structure
 matrix_1_raw = Elist_1_bg$E

# log2transformation
 matrix_1_log = log2(matrix_1_raw)

# Yugene normalization
 matrix_1_yugene = YuGene(matrix_1_log)
 dim(matrix_1_yugene)

# matrix_1_yugene is YuGene Class, and YuGene Class can't change to data.frame directly
 write.csv(matrix_1_yugene,"matrix_1_yugene.csv")
 matrix_1 = read.csv("matrix_1_yugene.csv")
 colnames(matrix_1) <- c("ID", group_list_1$geo_accession)  
 head(matrix_1)

# Get platform annotation information
 GPL10558 <- getGEO(filename ='GPL10558.soft')
 colnames(Table(GPL10558))
# Store ID,gene Symbol,  
 write.csv(Table(GPL10558)[,c(1,13)],"GPL10558.csv",row.names =F) 
 
# ID mapping (probeID -> EntrezID)
 matrix_1 <- as.data.frame(matrix_1) 
 genename_1 =read.csv("GPL10558.csv")¡¡
 express_1 = merge(x=matrix_1, y=genename_1,by="ID",all.x =T)
 express_1$ID =NULL   
        
# Remove the duplicate gene and retain the maximum expression of each gene(collapse)
 rowMeans = apply(express_1,1,function(x)mean(as.numeric(x),na.rm=T))
 express_1 = express_1[order(rowMeans,decreasing =T),]
 express_1 = express_1[!duplicated(express_1[,51]),]  #express_1 $ column 51 is gene symbol
 rownames(express_1)=express_1[,51]
 express_1=express_1[,-51]
 head(express_1)
 dim(express_1)
 
 genenames_1 = rownames(express_1)

################################## GSE15932: Affymetrix platform #########################################
 studyID = "GSE15932"
 destdir = "C:/Users/think/Desktop/Gene/raw data/RNA new/GSE15932"
 setwd(destdir)
 eSet_2 <- getGEO(studyID, destdir = ".")
 eSet_2 <- eSet_2[[1]]
 class(eSet_2)
 dim(eSet_2) 
 head(eSet_2)
 pdata = pData(eSet_2)
 write.csv(pdata,paste(studyID,"_matadata.csv")) # phenoData

# Construct grouping information 
 group_list_2 = subset(pdata,select=geo_accession) # Sample name
 group_list_2$condition = c(rep(c("T"),24),rep(c("N"),8))  # T:disease, N:normal grouping information
 group_list_2 = group_list_2[-c(9:16),] # Keep samples related to PC
 group_list_2¡¡


# Celfiles folder needs to be created in advance
 for (celFile in as.character(pData(eSet_2)$supplementary_file)) {
    storedFile <- sprintf("./CELfiles/%s", basename(celFile)) #"."very important!
    if (!file.exists(storedFile)) {
        download.file(url = celFile, destfile = storedFile)
    }
}

# Import CEL files
 storedFiles <- basename(as.character(pData(eSet_2)$supplementary_file))
 storedFiles <- sprintf("CELfiles/%s", storedFiles)
 affy_2 <- read.celfiles(filenames = storedFiles,phenoData = phenoData(eSet_2), sampleNames = rownames(pData(eSet_2)))
 class(affy_2)
 dim(affy_2)

# RMA four in one method (including log2 conversion) 
 affy_2_rma = rma(affy_2) # Only in this way can there be a probe name. RMA is a four in one method, which has been transformed into log2
 class(affy_2_rma)
 dim(affy_2_rma)

 md.pattern(exprs(affy_2_rma)) #No missing data

 matrix_2_raw <- exprs(affy_2_rma)

# Yugene normalization
 matrix_2_yugene = YuGene(matrix_2_raw)

# delete other cancers' samples, columns
 matrix_2_yugene = matrix_2_yugene[,-c(9:16)]
 head(matrix_2_yugene)
 dim(matrix_2_yugene)
 class(matrix_2_yugene)

 write.csv(matrix_2_yugene,"matrix_2_yugene.csv")
 matrix_2 = read.csv("matrix_2_yugene.csv")
 colnames(matrix_2) <- c("ID", group_list_2$geo_accession)  
 head(matrix_2)
 class(matrix_2)
 dim(matrix_2)

# Get platform annotation information
 GPL570 <- getGEO(filename ='GPL570.soft')
 colnames(Table(GPL570))
# Store ID,gene Symbol, column 1 and 11
 write.csv(Table(GPL570)[,c(1,11)],"GPL570.csv",row.names =F) 
 
# ID mapping (probeID -> gene Symbol)
 matrix_2 <- as.data.frame(matrix_2) 
 genename_2 =read.csv("GPL570.csv")¡¡
 express_2 = merge(x=matrix_2, y=genename_2, by="ID", all.x =T)
 head(express_2)
 express_2$ID = NULL   # delete probeID
 express_2 = na.omit(express_2) # There is a situation where the symbol does not exist
        
# Remove the duplicate gene and retain the maximum expression of each gene (collapses)
 rowMeans = apply(express_2,1,function(x)mean(as.numeric(x),na.rm=T))
 express_2 = express_2[order(rowMeans,decreasing =T),]
 express_2 = express_2[!duplicated(express_2[,25]),]  #express_2$ column 25 is gene symbol
 rownames(express_2)=express_2[,25]
 express_2 = express_2[,-25]
 head(express_2)
 dim(express_2)

 genenames_2 = rownames(express_2)

################################### GSE49515: Affymetrix #############################################
 studyID = "GSE49515"
 destdir = "C:/Users/think/Desktop/Gene/raw data/RNA new/GSE49515"
 setwd(destdir)
 eSet_3 <- getGEO(studyID, destdir = ".")
 eSet_3 <- eSet_3[[1]]
 class(eSet_3)
 dim(eSet_3)
 head(eSet_3)
 pdata = pData(eSet_3)
 write.csv(pdata,paste(studyID,"_matadata.csv")) # phenoData

# Construct grouping information  
 group_list_3 = subset(pdata,select=geo_accession) # Sample name
 group_list_3$condition = c(rep(c("N"),23),rep(c("T"),3))  # T:disease, N:normal grouping information
 group_list_3 = group_list_3[-c(1:13),] # Keep samples related to PC
 group_list_3¡¡

# Fetch the raw data (CEL files)
 for (celFile in as.character(pData(eSet_3)$supplementary_file)) {
    storedFile <- sprintf("./CELfiles/%s", basename(celFile)) 
    if (!file.exists(storedFile)) {
        download.file(url = celFile, destfile = storedFile)
    }
}

# Import CEL files
 storedFiles <- basename(as.character(pData(eSet_3)$supplementary_file))
 storedFiles <- sprintf("CELfiles/%s", storedFiles)
 affy_3 <- read.celfiles(filenames = storedFiles,phenoData = phenoData(eSet_3), sampleNames = rownames(pData(eSet_3)))
 class(affy_3)
 dim(affy_3)

# RMA four in one method (including log2 conversion)   
 affy_3_rma = rma(affy_3) 
 class(affy_3_rma)
 dim(affy_3_rma)

 md.pattern(exprs(affy_3_rma)) #No missing data

 matrix_3_raw <- exprs(affy_3_rma)

# Yugene normalization
 matrix_3_yugene = YuGene(matrix_3_raw)

# delete other cancers' samples, columns
 matrix_3_yugene = matrix_3_yugene[,-c(1:13)]
 head(matrix_3_yugene)
 dim(matrix_3_yugene)
 class(matrix_3_yugene)

 write.csv(matrix_3_yugene,"matrix_3_yugene.csv")
 matrix_3 = read.csv("matrix_3_yugene.csv")
 colnames(matrix_3) <- c("ID", group_list_3$geo_accession)  
 head(matrix_3)
 class(matrix_3)
 dim(matrix_3)

# Get platform annotation information
 GPL570 <- getGEO(filename ='GPL570.soft')
 colnames(Table(GPL570))
# Store ID,gene Symbol 
 write.csv(Table(GPL570)[,c(1,11)],"GPL570.csv",row.names =F) 
 
# ID mapping (probeID -> gene Symbol)
 matrix_3 <- as.data.frame(matrix_3) 
 genename_3 =read.csv("GPL570.csv")¡¡
 express_3 = merge(x=matrix_3, y=genename_3, by="ID", all.x =T)
 head(express_3)
 express_3$ID = NULL   # Remove probeID
 express_3 = na.omit(express_3)
        
# Remove the duplicate gene and retain the maximum expression of each gene (collapses)
 rowMeans = apply(express_3,1,function(x)mean(as.numeric(x),na.rm=T))
 express_3 = express_3[order(rowMeans,decreasing =T),]
 express_3 = express_3[!duplicated(express_3[,14]),]  #express_34 column 14 is gene symbol
 rownames(express_3)=express_3[,14]
 express_3 = express_3[,-14]
 head(express_3)
 dim(express_3)

 genenames_3 = rownames(express_3)

############################# Intersect common genes ##################################################
# intersect common genes
 temp = intersect(genenames_1,genenames_2)
 genenames = intersect(temp,genenames_3)
 length(genenames)

 m1 = express_1[genenames,]
 m2 = express_2[genenames,]
 m3 = express_3[genenames,]
 
 head(m1)
 head(m2)
 head(m3) # In this case, the number of rows of the three expression matrices is the same, and the probes represented by each row are the same

 c1 = group_list_1[,2]
 c1[c1=="T"] = 1
 c1[c1=="N"] = 0
 c2 = group_list_2[,2]
 c2[c2=="T"] = 1
 c2[c2=="N"] = 0
 c3 = group_list_3[,2]
 c3[c3=="T"] = 1
 c3[c3=="N"] = 0



###################################### Meta-analysis ##################################################
 RNA_esets = list(m1,m2,m3)
 RNA_classes = list(c1,c2,c3)
 RNA_genenames = as.factor(genenames)

# p-value combination
 RNA_res_1 = pvalcombination(RNA_esets,RNA_classes,BH = 0.01)
 str(RNA_res_1)
 length(RNA_res_1$Meta)
 DE_index_1 = RNA_res_1$Meta
 temp_1 = abs(RNA_res_1$TestStatistic[DE_index_1]) # Test statistics
 top_index_1= order(temp_1,decreasing=TRUE)[1:10]
 DEgenes_1 = genenames[top_index_1]
 DEgenes_1
#Add one step to get all the gene names
 top_index_11= order(temp_1,decreasing=TRUE)[1:278]
 DEgenes_11 = genenames[top_index_11]
 DEgenes_11
 data.frame(DEgenes_11)
 write.csv(DEgenes_11,"DEgenes_11.csv")
# Effect size combination
 RNA_res_2 = EScombination(RNA_esets,RNA_classes)
 str(RNA_res_2)
 length(RNA_res_2$Meta)
 DE_index_2 = RNA_res_2$Meta
 temp_2 = abs(RNA_res_2$TestStatistic[DE_index_2]) # Test statistics
 top_index_2 = order(temp_2,decreasing=TRUE)[1:10]
 DEgenes_2 = genenames[top_index_2]
 DEgenes_2
 #Add one step to get all the gene names
 top_index_21= order(temp_2,decreasing=TRUE)[1:328]
 DEgenes_21 = genenames[top_index_21]
 DEgenes_21
 data.frame(DEgenes_21)
 write.csv(DEgenes_21,"DEgenes_21.csv")

# intersect
 t1 = order(temp_1,decreasing=TRUE)[1:1000]
 t2 = order(temp_2,decreasing=TRUE)[1:1000]
 DE_intersect = intersect(t1,t2)
 DE_intersect
 DEgenes = genenames[DE_intersect]
 DEgenes
#Find all the difference genes
 
######################################## Merge #######################################

 write.csv(m1,"m1.csv")
 write.csv(m2,"m2.csv")
 write.csv(m3,"m3.csv")
 mc1=read.csv("m1.csv")
 mc2=read.csv("m2.csv")
 mc3=read.csv("m3.csv")
 mm=merge(x=mc1,y=mc2,by="X",all.x = T)
 m=merge(x=mm,y=mc3,by="X",all.x = T)
 c=c(c1,c2,c3)
 
 ################################### limma Degene #######################################
 

 #DGE for microarray by limma
 library('gplots')
 library('limma')
 foldChangelimma=0 #fold change=1 which means that the difference is twice
 padjlimma=0.01#padj=0.05 which means P value after correction is less than 0.05
 rawexprSetlimma=read.csv("m.csv",header=TRUE,row.names=1,check.names = FALSE)  
 #Read the matrix file. This is the input data path. Change it to your own file name#
 dim(rawexprSetlimma)
 #exprSet=log2(rawexprSet)
 par(mfrow=c(1,2))
 boxplot(data.frame(rawexprSetlimma),col="blue") ## Draw box chart and compare data distribution
 rawexprSetlimma[1:5,1:5]
 grouplimma <- read.csv("m.delete.csv",header=TRUE,row.names=1,check.names = FALSE)
 grouplimma <- grouplimma[,1] #Define the comparison group and modify it according to the number of cancer and normal samples#
 designlimma <- model.matrix(~0+factor(grouplimma))#Set group as a model matrix#
 colnames(designlimma)=levels(factor(grouplimma))
 rownames(designlimma)=colnames(rawexprSetlimma)
#computational procedure
 fit <- lmFit(rawexprSetlimma,designlimma)
 cont.matrix<-makeContrasts(paste0(unique(grouplimma),collapse = "-"),levels = designlimma)
 fit2=contrasts.fit(fit,cont.matrix)
 fit2 <- eBayes(fit2)  ## default no trend !!!
 ##eBayes() with trend=TRUE
 tempOutput = topTable(fit2,coef=1,n=Inf,adjust="BH") 
 nrDEG = na.omit(tempOutput) 
 
 #Output result
 allDiff <- nrDEG
 diff=allDiff
 write.csv(diff, "limmaOut.csv")
 diffSig = diff[(diff$P.Value < padjlimma & (diff$logFC>foldChangelimma | diff$logFC<(-foldChangelimma))),]#Select significant differences in screening#
 #write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)#Output to diffsig file with significant difference expression#
 write.csv(diffSig, "diffSig.csv")
 diffUp = diff[(diff$P.Value < padjlimma & (diff$logFC>foldChangelimma)),]#foldchange>0 is up£¬foldchange<0 is down#
 #write.table(diffUp, file="up.xls",sep="\t",quote=F)#39-42Input up and down files into up and down files respectively#
 write.csv(diffUp, "diffUp.csv")
 diffDown = diff[(diff$P.Value < padjlimma & (diff$logFC<(-foldChangelimma))),]
 #write.table(diffDown, file="down.xls",sep="\t",quote=F)
 write.csv(diffDown, "diffDown.csv")
 
 

 ################################### edgeR Degene #######################################
 
 library("edgeR")
 foldChangeedgeR=0
 padjedgeR=0.05#Lines 6 and 7 represent the filtering standard, and fold change = 1 means the difference is twice£¬
 rt=read.csv("m.csv",header=TRUE,row.names=1,check.names = FALSE)
 #rownames(rt)
 #exp=rt[,1:ncol(rt)] #The first column to the last column is the data expressed#
 #dimnames=list(rownames(exp),colnames(exp))
 data=rt
 #data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
 #It means to convert quoted data into numerical value, useless#
 #data=data[rowMeans(data)>1,] #Remove low expression data#
 groupedgeR <- read.csv("m.delete.csv",header=TRUE,row.names=1,check.names = FALSE)
 groupedgeR <- groupedgeR[,1] #Define the comparison group and modify it according to the number of cancer and normal samples#
 designedgeR <- model.matrix(~0+factor(groupedgeR))#Set group as a model matrix#
 data[data < 0] <- 0
 y <- DGEList(counts=data,group=groupedgeR) #Which are normal and which are cancer samples for edger to identify#
 y <- calcNormFactors(y) #Factor correction#
 y <- estimateCommonDisp(y)
 #This step and the next step estimate the coefficient of variation, i.e. estimate the variance; estimate the degree of internal difference, 
 #to see whether the difference between groups is greater than the internal difference, if so, select the difference gene#
 y <- estimateTagwiseDisp(y)
 et <- exactTest(y,pair = c("h","d"))
 topTags(et)
 ordered_tags <- topTags(et, n=100000) #show first 100000#
 allDiff=ordered_tags$table
 allDiff=allDiff[is.na(allDiff$FDR)==FALSE,]
 diff=allDiff
 newData=y$pseudo.counts
 write.table(diff,file="edgerOut.xls",sep="\t",quote=F)#First, put all the differences into the file called edgerout#
 
 
 diffSig = diff[(diff$FDR < padjedgeR & (diff$logFC>foldChangeedgeR | diff$logFC<(-foldChangeedgeR))),]
 write.table(diffSig, file="diffSig.xls",sep="\t",quote=F)#Output to diffsig file with significant difference expression#
 
 
 
 ################################### t-test Degene #######################################
 
 
 t=read.csv("m.csv",header=TRUE,row.names=1,check.names = FALSE)
 #Pre generate 2 vectors with the same length and the same number of lines of input file, which will be used to store p value and difference multiple (log2fc)
 
 Pvalue<-c(rep(0,nrow(t)))
 log2_FC<-c(rep(0,nrow(t)))
 
 for(i in 1:nrow(t)){
   if(sd(t[i,1:66])==0&&sd(t[i,67:87])==0){
     Pvalue[i] <-"NA"
     log2_FC[i]<-"NA"
   }else{
     y=t.test(as.numeric(t[i,1:66]),as.numeric(t[i,67:87]))
     Pvalue[i]<-y$p.value
     log2_FC[i]<-log2((mean(as.numeric(t[i,1:66]))+0.001)/(mean(as.numeric(t[i,67:87]))+0.001)) 
   }
 }
 # FDR correction of P value
 fdr<-c(rep(0,nrow(t)))
 fdr[i]=p.adjust(Pvalue[i], 'BH')
 
 # Add log2fc, P value and FDR after the original file, 3 columns in total;
 out<-cbind(t,log2_FC,Pvalue,fdr)
 write.table(out,file="ttest.out.xls",quote=FALSE,sep="\t",row.names=FALSE)
 
 ################################### ANOVA Degene #######################################
 
 a=read.csv("m.csv",header=TRUE,row.names=1,check.names = FALSE)
 
 #Pre generate 2 vectors with the same length and the same number of lines of input file, which will be used to store p value and difference multiple (log2fc)
 
 Pvalue<-c(rep(0,nrow(a)))
 
 log2_FC<-c(rep(0,nrow(a)))
 
 #Determine grouping information
 
 type<-factor(c(rep('c',66),rep('m',21)))
 
 #2-4 columns were treatment group 1, 5-7 columns were treatment group 2;
 #Each row will be tested for variance using a cycle
 #The expression level of both groups was 0, not tested;
 #P value and log2fc calculated from each line will be added to the last two columns of the original file;
 #When calculating log2fc, add 0.001 to the mean value of each group to prevent bugs caused by 0 denominator;
 
 for(i in 1:nrow(a)){
   if(sum(a[i,1:66])==0&&sum(a[i,67:87])==0){
     Pvalue[i] <- 'NA'
     log2_FC[i]<- 'NA'
   }else{
     y=aov(as.numeric(a[i,1:87])~type)
     Pvalue[i]<-summary(y)[[1]][,5][1]
     log2_FC[i]<-log2((mean(as.numeric(a[i,1:66]))+0.001)/(mean(as.numeric(a[i,67:87]))+0.001))
   }
 }
 # FDR correction of P value
 fdr=p.adjust(Pvalue, 'BH')
 
 # Add log2fc, P value and FDR after the original file, 3 columns in total;
 out<-cbind(a,log2_FC,Pvalue,fdr)
 #write.table(out,file=¡±anova.out.xls¡±,quote=FALSE,sep=¡±\t¡±,row.names=FALSE)
 
 ################################### DESeq2 Degene #######################################
 
 
 us_count<-read.csv("m.csv",header=TRUE,row.names=1,check.names = FALSE)

 
 condition<-factor(c(rep('c',66),rep('m',21)))
 
 coldata<-data.frame(row.names=colnames(us_count),condition)  
 library(DESeq2) 
 ##Building DDS matrix
 us_count[us_count < 0] <- 0
 us_count<-round(us_count,digits=0) #Rounding input data
 
 dds<-DESeqDataSetFromMatrix(us_count,coldata,design=~condition)
 head(dds) #View the constructed matrix
 
 ##Conduct variance analysis
 #DDS < - deseq (DDS) standardize the original DDS
 #Resultsnames (DDS) - view result names
 #Results < - results (DDS) extract the results with the results function and assign it to the res variable
 #Summary (RES) view results
 plotMA(res,ylim=c(-2,2)) 
 mcols(res,use.names=TRUE)
 
 plot(res$log2FoldChange,res$pvalue) #Mapping volcanoes
 
 #Extraction of differential genes
 res <- res[order(res$padj),]
 resdata <-merge(as.data.frame(res),as.data.frame(counts(dds,normalize=TRUE)),by="row.names",sort=FALSE)
 deseq_res<-data.frame(resdata)
 #up_diff_result<-subset(deseq_res,padj < 0.05 & (log2FoldChange > 1)) #Extraction of up regulated differentially expressed genes
 #down_diff_result<-subset(deseq_res,padj < 0.05 & (log2FoldChange < -1)) #Extraction of down regulated differentially expressed genes
 
 #write.csv(up_diff_result,"C:\\Users\\admin\\Desktop\\up_day0_VS_day1_diff_results.csv") 
# write.csv(down_diff_result,"C:\\Users\\admin\\Desktop\\down_day0_VS_day1_diff_results.csv") 

 
 
 ################################### samr Degene #######################################
 
 c=read.csv("c.csv",header=TRUE,row.names=1,check.names = FALSE)
 c <- c[,1]
 c=as.numeric(c)
 sam=read.csv("m.csv",header=TRUE,row.names=1,check.names = FALSE)
 sam[sam < 0] <- 0
 data<- data.matrix(sam)
 
 suppressPackageStartupMessages(library(samr))
 #data1=list(x=data,y=as.numeric(as.factor(c)), 
           #geneid=as.character(1:nrow(data)),
          # genenames=rownames(data), 
          # logged2=FALSE
 #)
 #samr.obj<-samr(data1, resp.type="Two class unpaired", nperms=100)
 samfit<-SAM(data,c,resp.type="Two class unpaired")
 #pv=samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
 #data.frame(pv)
 #rite.csv(pv,"pv.csv")
 options(max.print=2000) 
 print(samfit)
 sa=print(samfit)
 samfit$siggenes.table$genes.up
 write.csv(samfit$siggenes.table$genes.up,"up.csv")
 write.csv(samfit$siggenes.table$genes.lo,"down.csv")
 
 ################################### WGCNA #######################################
 
 library(affy)
 library(affyPLM)
 library(RColorBrewer)
 library(impute)
 library(limma)
 library(pheatmap)
 library(ggplot2)
 library(WGCNA)
 #Load the expression spectrum and traits file into R and prepare the format required by WGCNA
 datExpr=read.csv("DEgenenew.csv")
 # Manipulate file so it matches the format WGCNA needs,Format conforming to WGCNA under processing
 row.names(datExpr) = datExpr$DEgenes
 datExpr$DEgenes=NULL
 datExpr = as.data.frame(t(datExpr))# now samples are rows and genes are columns
 dim(datExpr)
 gsg = goodSamplesGenes(datExpr, verbose = 3)
 gsg$allOK
 #This must be true, otherwise execute the following code
 #if (!gsg$allOK)
 #{if (sum(!gsg$goodGenes)>0)
  # printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
   #if (sum(!gsg$goodSamples)>0)
    # printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
   #datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
 #}
 #Create an object called "datTraits" that contains your trait data
 datTraits = read.csv("c.csv")
 rownames(datTraits) = datTraits$gsm
 datTraits$gsm = NULL
 table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly,otherwise your names are out of order
 #This must be true, otherwise the following cannot be executed. The reason for the error should be on the naming.
 #So far, all the files have been uploaded, so wait to enter the analysis process.
 #Remember to write the file to the local disk
 save(datExpr, datTraits, file="SamplesAndTraits.RData")
 #load("SamplesAndTraits.RData")
 
 
 # Cluster samples by expression
 A = adjacency(t(datExpr),type="signed") # this calculates the whole network connectivity
 k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
 Z.k = scale(k)
 thresholdZ.k = -2.5 # often -2.5
 outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
 sampleTree = flashClust(as.dist(1-A), method = "average")#Unable to output
 # Convert traits to a color representation where red indicates high values
 traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
 dimnames(traitColors)[[2]] = paste(names(datTraits))
 datColors = data.frame(outlier = outlierColor,traitColors)
 plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                     colors=datColors,main="Sample Dendrogram and Trait Heatmap")
 
 # Choose a soft threshold power- USE A SUPERCOMPUTER IRL -----------
 powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
 sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function
 sizeGrWindow(9,5)
 par(mfrow= c(1,2))
 cex1=0.9
 plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
 text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
 #In some papers, R2 is greater than 0.85
 abline(h=0.90, col="red")
 plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
 text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
 #Check whether the network approaches scale free under the selected threshold
 #Error, unable to run
 ADJ1=abs(cor(datExpr,use="p"))^8
 k=as.vector(apply(ADJ1,2,sum,na.rm=T))
 sizeGrWindow(10,5)
 par(mfrow=c(1,2))
 hist(k)
 scaleFreePlot(k,main="Check scale free topology")
 enableWGCNAThreads()
 #Use power=1**8****Adjacency Matrix Building**
 softPower = 18
 #Adjacency Matrix Building
 adjacency = adjacency(datExpr, power = softPower, type = "signed") #specify network type
 #translate the adjacency into topological overlap matrix and calculate the corresponding #dissimilarity:½«ÁÚ½Ó¾ØÕó×ª»»ÎªÍØÆË¾ØÕó
 TOM = TOMsimilarity(adjacency, TOMType="signed") # specify network type OGM2hours
 dissTOM = 1-TOM
 # Generate Modules --------------------------------------------------------
 # Generate a clustered gene tree
 geneTree = flashClust(as.dist(dissTOM), method="average")
 plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)
 #This sets the minimum number of genes to cluster into a module
 #Dynamic cutting tree cutting module, choose dynamic mixed cutting method
 minModuleSize = 30
 dynamicMods = cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize = minModuleSize)
 dynamicColors= labels2colors(dynamicMods)
 #Merging similar coexpression networks is equivalent to clustering the modules selected from the previous part. 
 #After the pruning height is determined by using the functionmoduleeigegnes(), merge models by using the function mergeclosemodules()
 #set a threhold for merging modules. In this example we are not merging so MEDissThres=0.0
 MEList= moduleEigengenes(datExpr, colors= dynamicColors,softPower = 18)
 MEs= MEList$eigengenes
 MEDiss= 1-cor(MEs)
 METree= flashClust(as.dist(MEDiss), method= "average")
 #plots tree showing how the eigengenes cluster together
 plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
 #Select those with a correlation coefficient greater than 0.8 to merge (Wang Pan's 0.75, below is not 0.2 but 0.25)
 MEDissThres = 0.2
 merge = mergeCloseModules(datExpr, dynamicColors, cutHeight= MEDissThres, verbose =3)
 mergedColors = merge$colors
 sizeGrWindow(12,9)
 mergedMEs = merge$newMEs
 save(dynamicMods, MEList, MEs, MEDiss, METree, file= "Network_allSamples_signed_RLDfiltered.RData")
 #plot dendrogram with module colors below it
 plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
 moduleColors = mergedColors
 colorOrder = c("grey", standardColors(50))
 moduleLabels = match(moduleColors, colorOrder)-1
 MEs = mergedMEs
 save(MEs, moduleLabels, moduleColors, geneTree, file= "Network_allSamples_signed_nomerge_RLDfiltered.RData")
 
 
 #Build a network (build modules and connect with external information)
 
 #Define number of genes and samples
 nGenes = ncol(datExpr)
 nSamples = nrow(datExpr)
 
 #Recalculate MEs with color labels
 MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
 MEs = orderMEs(MEs0)
 moduleTraitCor = cor(MEs, datTraits, use= "p")
 moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
 
 #Print correlation heatmap between modules and traits
 textMatrix= paste(signif(moduleTraitCor, 2), "\n(", 
                   signif(moduleTraitPvalue, 1), ")", sep= "")
 dim(textMatrix)= dim(moduleTraitCor)
 par(mar= c(6, 8.5, 3, 3))
 #display the corelation values with a heatmap plot
 labeledHeatmap(Matrix= moduleTraitCor, 
                xLabels= names(datTraits), 
                yLabels= names(MEs), 
                ySymbols= names(MEs), 
                colorLabels= FALSE, 
                colors= blueWhiteRed(50), 
                textMatrix= textMatrix, 
                setStdMargins= FALSE, 
                cex.text= 0.5, 
                zlim= c(-1,1), 
                main= paste("Module-trait relationships"))
 
 #Thermogram of different modules and expression of key genes
 person=cor(datExpr,use = 'p')
 corr<-TOM
 Colors<-mergedColors
 colnames(corr)<-colnames(datExpr)
 rownames(corr)<-colnames(datExpr)
 names(Colors)<-colnames(datExpr)
 colnames(person)<-colnames(datExpr)
 rownames(person)<-colnames(datExpr)
 umc = unique(mergedColors)
 lumc = length(umc)
 #change the i to be the number which you want to plot
 ##Change the number on the yellow background to plot the diagram of each module
 for (i in c(1:lumc)){
   if(umc[i]== "grey"){
     next
   }
   ME=MEs[, paste("ME",umc[3], sep="")]
   par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
   plotMat(t(scale(datExpr[,Colors==umc[3]])),nrgcols=30,rlabels=F,rcols=umc[3], main=umc[3], cex.main=3)
   par(mar=c(5, 4.2, 0, 0.7))
   barplot(ME, col=umc[3], main="", cex.main=2,ylab="eigengene expression",xlab="array sample")
 }
 
 #Heat map of gene co expression network
 kME=signedKME(datExpr, mergedMEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use = 'p'")
 
 if (dim(datExpr)[2]>=15000) nSelect=15000 else nSelect=dim(datExpr)[2]
 set.seed(1)
 select = sample(nGenes, size = nSelect)
 selectTOM = dissTOM[select, select]
 selectTree = hclust(as.dist(selectTOM), method = "average")
 selectColors = moduleColors[select]
 plotDiss = selectTOM^7
 #the following code will take about   mins
 #The picture below will take a long time, more than 3 hours
 TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot")
 
 
 #Module correlation heat map
 MEs = moduleEigengenes(datExpr, Colors)$eigengenes
 MET = orderMEs(MEs)
 par(mfrow=c(1,1))
 par(mar= c(5, 10.5, 3, 3))
 plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
 
 
 MM= as.data.frame(cor(datExpr, MEs, use ="p"))
 write.csv(MM,"MM.csv")
 
 
 
 
 
 
 
 
 