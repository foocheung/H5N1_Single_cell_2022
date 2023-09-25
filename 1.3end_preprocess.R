library("Seurat") #load Seurat 4
library("dplyr")
library("matrixStats")
library('tidyverse')
library('dsb')

##qrsh -l mem_free=160G,h_vmem=160G
###/sysapps/cluster/software/Anaconda2/5.3.0/bin/R
rawdir="/hpcdata/sg/sg_data/illumina_NCI_runs/CHI_H5citeseq/Batch1/3primeCITE/"
rawdir2="/hpcdata/sg/sg_data/illumina_NCI_runs/CHI_H5citeseq/Batch2/3primeCITE2/"
rawdir3 = "/hpcdata/sg/sg_data/illumina_NCI_runs/CHI_H5citeseq/Batch3and4/Flowcell_HHY33DSXX/"
rawdir4 = "/hpcdata/sg/sg_data/illumina_NCI_runs/CHI_H5citeseq/Batch3and4/Flowcell_HTGC2DSXX/"
##CHECK IS THIS UNSORTED

B1LanesUnsorted = c("01","02","03","04","05","06","07","08")
B2LanesUnsorted = c("01","02","03","04","05","06","07","08")
B3LanesUnsorted = c("01","02","03","04","05","06","07","08", "09", "10", "11", "12")
B4LanesUnsorted = c("01","02","03","04","05","06","07","08", "09", "10", "11", "12")




##Process Batch 1 first

B1_US_data = list()
B1_US_SeuratObj = list()

B2_US_data = list()
B2_US_SeuratObj = list()

B3_US_data = list()
B3_US_SeuratObj = list()

B4_US_data = list()
B4_US_SeuratObj = list()

######################################
##BATCH 1
for(i in 1:length(B1LanesUnsorted)){
  B1_US_data[[i]] = Read10X_h5(paste(rawdir,"multi_config_",B1LanesUnsorted[[i]],
                                     "/outs/per_sample_outs/","multi_config_",
                                     B1LanesUnsorted[[i]],
                                     "/count/sample_feature_bc_matrix.h5", 
                                     sep=""))
}


for(i in 1:length(B1LanesUnsorted)){
  B1_US_SeuratObj[[i]] <- CreateSeuratObject(counts = B1_US_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 0)

  B1_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Antibody Capture'[1:113,])

  ##B1_US_data[[i]]$'Antibody Capture'[c(1:10),]
  ##SEE ONLY 4:6 is being used
  B1_US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B1_US_data[[i]]$'Custom'[c(1:5),])
  B1_US_SeuratObj[[i]] <- RenameCells(B1_US_SeuratObj[[i]], new.names = paste(substr(colnames(B1_US_SeuratObj[[i]]), start = 1, stop = 17),B1LanesUnsorted[[i]], sep = ""))
  B1_US_SeuratObj[[i]]$Batch  <- rep("B1UnSort", length(colnames(B1_US_SeuratObj[[i]])))
  B1_US_SeuratObj[[i]]$Lane  <- rep(B1LanesUnsorted[[i]], length(colnames(B1_US_SeuratObj[[i]])))

}

######################################
##BATCH 2
for(i in 1:length(B2LanesUnsorted)){
  B2_US_data[[i]] = Read10X_h5(paste(rawdir2,"multi_config_",
                                     B2LanesUnsorted[[i]],
                                     "/outs/per_sample_outs/",
                                     "multi_config_",B2LanesUnsorted[[i]],
                                     "/count/sample_feature_bc_matrix.h5", sep=""))
}

##rownames(B2_US_data[[i]]$'Antibody Capture')
##rownames(B2_US_data[[i]]$Custom)
for(i in 1:length(B2LanesUnsorted)){
  B2_US_SeuratObj[[i]] <- CreateSeuratObject(counts = B2_US_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 0)
  
  B2_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B2_US_data[[i]]$'Antibody Capture'[1:113,])
  
  ##B2_US_data[[i]]$'Antibody Capture'[c(1:113),]
  ##SEE ONLY 4:6 or 1:6is being used
  B2_US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B2_US_data[[i]]$'Custom'[c(6:10),])
  B2_US_SeuratObj[[i]] <- RenameCells(B2_US_SeuratObj[[i]], new.names = paste(substr(colnames(B2_US_SeuratObj[[i]]), start = 1, stop = 17),B2LanesUnsorted[[i]], sep = ""))
  B2_US_SeuratObj[[i]]$Batch  <- rep("B2UnSort", length(colnames(B2_US_SeuratObj[[i]])))
  B2_US_SeuratObj[[i]]$Lane  <- rep(B2LanesUnsorted[[i]], length(colnames(B2_US_SeuratObj[[i]])))
  
}
######################################

######################################
##BATCH 3

for(i in 1:length(B3LanesUnsorted)){
  B3_US_data[[i]] = Read10X_h5(paste(rawdir3,"multi_config_",
                                     B3LanesUnsorted[[i]],
                                     "/outs/per_sample_outs/",
                                     "multi_config_",B3LanesUnsorted[[i]],
                                     "/count/sample_feature_bc_matrix.h5", sep=""))
}

##rownames(B1_US_data[[i]]$'Antibody Capture')
##rownames(B1_US_data[[i]]$Custom)
for(i in 1:length(B3LanesUnsorted)){
  B3_US_SeuratObj[[i]] <- CreateSeuratObject(counts = B3_US_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 0)
  
  B3_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B3_US_data[[i]]$'Antibody Capture'[1:113,])
  
  ##B3_US_data[[i]]$'Antibody Capture'[c(1:113),]
  ##SEE ONLY 4:6 or 1:6is being used
  #B3_US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B3_US_data[[i]]$'Custom'[c(7,8,9,10,11,12),])
  B3_US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B3_US_data[[i]]$'Custom'[c(1:6),])
  B3_US_SeuratObj[[i]] <- RenameCells(B3_US_SeuratObj[[i]], new.names = paste(substr(colnames(B3_US_SeuratObj[[i]]), start = 1, stop = 17),B3LanesUnsorted[[i]], sep = ""))
  B3_US_SeuratObj[[i]]$Batch  <- rep("B3UnSort", length(colnames(B3_US_SeuratObj[[i]])))
  B3_US_SeuratObj[[i]]$Lane  <- rep(B3LanesUnsorted[[i]], length(colnames(B3_US_SeuratObj[[i]])))
  
}
######################################

######################################
##BATCH 4

for(i in 1:length(B4LanesUnsorted)){
  B4_US_data[[i]] = Read10X_h5(paste(rawdir4,"multi_config_",
                                     B4LanesUnsorted[[i]],
                                     "/outs/per_sample_outs/",
                                     "multi_config_",B4LanesUnsorted[[i]],
                                     "/count/sample_feature_bc_matrix.h5", sep=""))
}

##rownames(B1_US_data[[i]]$'Antibody Capture')
##rownames(B1_US_data[[i]]$Custom)
for(i in 1:length(B4LanesUnsorted)){
  B4_US_SeuratObj[[i]] <- CreateSeuratObject(counts = B4_US_data[[i]]$'Gene Expression', assay = "RNA", min.feature = 0)
  
  B4_US_SeuratObj[[i]][["CITE"]] <- CreateAssayObject(counts = B4_US_data[[i]]$'Antibody Capture'[1:113,])
  
  ##B4_US_data[[i]]$'Antibody Capture'[c(1:113),]
  ##SEE ONLY 4:6 or 1:6is being used
  B4_US_SeuratObj[[i]][["HTO"]] <- CreateAssayObject(counts = B4_US_data[[i]]$'Custom'[c(7,8,9,10,11,12),])
  B4_US_SeuratObj[[i]] <- RenameCells(B4_US_SeuratObj[[i]], new.names = paste(substr(colnames(B4_US_SeuratObj[[i]]), start = 1, stop = 17),B4LanesUnsorted[[i]], sep = ""))
  B4_US_SeuratObj[[i]]$Batch  <- rep("B4UnSort", length(colnames(B4_US_SeuratObj[[i]])))
  B4_US_SeuratObj[[i]]$Lane  <- rep(B4LanesUnsorted[[i]], length(colnames(B4_US_SeuratObj[[i]])))
  
}
######################################





names(B1_US_data) = names(B1_US_SeuratObj) = B1LanesUnsorted
names(B2_US_data) = names(B2_US_SeuratObj) = B2LanesUnsorted


names(B3_US_data) = names(B3_US_SeuratObj) = B3LanesUnsorted
names(B4_US_data) = names(B4_US_SeuratObj) = B4LanesUnsorted


B1_US_merge = merge(B1_US_SeuratObj[[1]], B1_US_SeuratObj[2:length(B1LanesUnsorted)])
B2_US_merge = merge(B2_US_SeuratObj[[1]], B2_US_SeuratObj[2:length(B2LanesUnsorted)])
B3_US_merge = merge(B3_US_SeuratObj[[1]], B3_US_SeuratObj[2:length(B2LanesUnsorted)])
B4_US_merge = merge(B4_US_SeuratObj[[1]], B4_US_SeuratObj[2:length(B4LanesUnsorted)])
##gets killed if you dont get enough h_vmem=160G

##DOES THIS PLOT LOOKS RIGHT ??
pdf("Violinplot1.pdf")
VlnPlot(B2_US_merge, c("nCount_CITE", "nCount_HTO", "nCount_RNA"), group.by = "Lane", pt.size=0)
dev.off()


B1_US_merge <- subset(B1_US_merge, subset = nCount_CITE < 30000)
B1_US_merge <- subset(B1_US_merge, subset = nCount_HTO < 15001)
B1_US_merge <- NormalizeData(B1_US_merge, assay = "HTO", normalization.method = "CLR", margin = 2)
B1_US_merge <- ScaleData(B1_US_merge, assay = "HTO", model.use = "linear")
B1_US_merge = HTODemux(B1_US_merge, positive.quantile = 0.85)
Idents(B1_US_merge) <- "hash.ID"


B2_US_merge <- subset(B2_US_merge, subset = nCount_CITE < 30000)
B2_US_merge <- subset(B2_US_merge, subset = nCount_HTO < 15001)
B2_US_merge <- NormalizeData(B2_US_merge, assay = "HTO", normalization.method = "CLR", margin = 2)
B2_US_merge <- ScaleData(B2_US_merge, assay = "HTO", model.use = "linear")
B2_US_merge = HTODemux(B2_US_merge, positive.quantile = 0.85)
Idents(B2_US_merge) <- "hash.ID"


B3_US_merge <- subset(B3_US_merge, subset = nCount_CITE < 30000)
B3_US_merge <- subset(B3_US_merge, subset = nCount_HTO < 15001)
B3_US_merge <- NormalizeData(B3_US_merge, assay = "HTO", normalization.method = "CLR", margin = 2)
B3_US_merge <- ScaleData(B3_US_merge, assay = "HTO", model.use = "linear")
B3_US_merge = HTODemux(B3_US_merge, positive.quantile = 0.85)
Idents(B3_US_merge) <- "hash.ID"


B4_US_merge <- subset(B4_US_merge, subset = nCount_CITE < 30000)
B4_US_merge <- subset(B4_US_merge, subset = nCount_HTO < 15001)
B4_US_merge <- NormalizeData(B4_US_merge, assay = "HTO", normalization.method = "CLR", margin = 2)
B4_US_merge <- ScaleData(B4_US_merge, assay = "HTO", model.use = "linear")
B4_US_merge = HTODemux(B4_US_merge, positive.quantile = 0.85)
Idents(B4_US_merge) <- "hash.ID"




##################################################


B1_US_demuxbestList = list()
for(i in 1:length(B1LanesUnsorted)){
  B1_US_demuxbestList[[i]] = read.table(paste(rawdir,"DEMUX/multi_config_",B1LanesUnsorted[[i]],".best", sep = ""), sep = "\t", header = TRUE)
}

for(i in 1:length(B1LanesUnsorted)){
  B1_US_demuxbestList[[i]]$NewBarcode = paste(substr(B1_US_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B1LanesUnsorted[[i]], sep = "")
}

B1_US_demuxbestdf <- plyr::ldply(B1_US_demuxbestList, data.frame)
length(which(colnames(B1_US_merge) %in% B1_US_demuxbestdf$NewBarcode))
length(setdiff(colnames(B1_US_merge), B1_US_demuxbestdf$NewBarcode))
rownames(B1_US_demuxbestdf) <- B1_US_demuxbestdf$NewBarcode
B1_US_merge <- subset(B1_US_merge, cells =  B1_US_demuxbestdf$NewBarcode)
B1_US_merge <- AddMetaData(B1_US_merge, metadata = B1_US_demuxbestdf[colnames(B1_US_merge),])




B2_US_demuxbestList = list()
for(i in 1:length(B2LanesUnsorted)){
  B2_US_demuxbestList[[i]] = read.table(paste(rawdir2,"DEMUX/multi_config_",B2LanesUnsorted[[i]],".best", sep = ""), sep = "\t", header = TRUE)
}

for(i in 1:length(B2LanesUnsorted)){
  B2_US_demuxbestList[[i]]$NewBarcode = paste(substr(B2_US_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B2LanesUnsorted[[i]], sep = "")
}

B2_US_demuxbestdf <- plyr::ldply(B2_US_demuxbestList, data.frame)
length(which(colnames(B2_US_merge) %in% B2_US_demuxbestdf$NewBarcode))
length(setdiff(colnames(B2_US_merge), B2_US_demuxbestdf$NewBarcode))
rownames(B2_US_demuxbestdf) <- B2_US_demuxbestdf$NewBarcode
B2_US_merge <- subset(B2_US_merge, cells =  B2_US_demuxbestdf$NewBarcode)
B2_US_merge <- AddMetaData(B2_US_merge, metadata = B2_US_demuxbestdf[colnames(B2_US_merge),])






B3_US_demuxbestList = list()
for(i in 1:length(B3LanesUnsorted)){
  B3_US_demuxbestList[[i]] = read.table(paste(rawdir3,"DEMUX/multi_config_",B3LanesUnsorted[[i]],".best", sep = ""), sep = "\t", header = TRUE)
}

for(i in 1:length(B3LanesUnsorted)){
  B3_US_demuxbestList[[i]]$NewBarcode = paste(substr(B3_US_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B3LanesUnsorted[[i]], sep = "")
}

B3_US_demuxbestdf <- plyr::ldply(B3_US_demuxbestList, data.frame)
length(which(colnames(B3_US_merge) %in% B3_US_demuxbestdf$NewBarcode))
length(setdiff(colnames(B3_US_merge), B3_US_demuxbestdf$NewBarcode))
rownames(B3_US_demuxbestdf) <- B3_US_demuxbestdf$NewBarcode
B3_US_merge <- subset(B3_US_merge, cells =  B3_US_demuxbestdf$NewBarcode)
B3_US_merge <- AddMetaData(B3_US_merge, metadata = B3_US_demuxbestdf[colnames(B3_US_merge),])




B4_US_demuxbestList = list()
for(i in 1:length(B4LanesUnsorted)){
  B4_US_demuxbestList[[i]] = read.table(paste(rawdir4,"DEMUX/multi_config_",B4LanesUnsorted[[i]],".best", sep = ""), sep = "\t", header = TRUE)
}

for(i in 1:length(B4LanesUnsorted)){
  B4_US_demuxbestList[[i]]$NewBarcode = paste(substr(B4_US_demuxbestList[[i]]$BARCODE, start = 1, stop = 17),B4LanesUnsorted[[i]], sep = "")
}

B4_US_demuxbestdf <- plyr::ldply(B4_US_demuxbestList, data.frame)
length(which(colnames(B4_US_merge) %in% B4_US_demuxbestdf$NewBarcode))
length(setdiff(colnames(B4_US_merge), B4_US_demuxbestdf$NewBarcode))
rownames(B4_US_demuxbestdf) <- B4_US_demuxbestdf$NewBarcode
B4_US_merge <- subset(B4_US_merge, cells =  B4_US_demuxbestdf$NewBarcode)
B4_US_merge <- AddMetaData(B4_US_merge, metadata = B4_US_demuxbestdf[colnames(B4_US_merge),])



########################################################################################





##############################################

SNG.US <- list(subset(B1_US_merge, subset = DROPLET.TYPE == "SNG"),
               subset(B2_US_merge, subset = DROPLET.TYPE == "SNG"),
               subset(B3_US_merge, subset = DROPLET.TYPE == "SNG"),
               subset(B4_US_merge, subset = DROPLET.TYPE == "SNG")
               )


NEG.US <- list(subset(B1_US_merge, subset = DROPLET.TYPE == "AMB"),
               subset(B2_US_merge, subset = DROPLET.TYPE == "AMB"),
               subset(B3_US_merge, subset = DROPLET.TYPE == "AMB"),
               subset(B4_US_merge, subset = DROPLET.TYPE == "AMB")
)

names(SNG.US) <-  names(NEG.US) <- c("B1", "B2", "B3","B4")


mito.genes = grep(pattern = "^MT-", x = rownames(NEG.US[[1]]), value = TRUE)
(mean(NEG.US[[1]]$nFeature_RNA) + 3*sd(NEG.US[[1]]$nFeature_RNA))
(mean(NEG.US[[2]]$nFeature_RNA) + 3*sd(NEG.US[[2]]$nFeature_RNA))
(mean(NEG.US[[3]]$nFeature_RNA) + 3*sd(NEG.US[[3]]$nFeature_RNA))
(mean(NEG.US[[4]]$nFeature_RNA) + 3*sd(NEG.US[[4]]$nFeature_RNA))

for(i in 1:length(SNG.US)){
  SNG.US[[i]] <- AddMetaData(object = SNG.US[[i]], metadata = Matrix::colSums(SNG.US[[i]][mito.genes,])/Matrix::colSums(SNG.US[[i]]), col.name = "percent.mito")
  }

pdf("Violinplot2.pdf")
VlnPlot(SNG.US[[1]], features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mito"), pt.size = 0, ncol = 2)
VlnPlot(SNG.US[[2]], features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mito"), pt.size = 0, ncol = 2)
VlnPlot(SNG.US[[3]], features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mito"), pt.size = 0, ncol = 2)
VlnPlot(SNG.US[[4]], features= c("nCount_RNA","nCount_HTO","nCount_CITE", "nFeature_RNA","percent.mito"), pt.size = 0, ncol = 2)
dev.off()

pdf("Hist1.pdf")
hist(NEG.US[[1]]$nFeature_RNA, breaks = 100)
for(i in 1:length(NEG.US)){
  NEG.US[[i]] = subset(NEG.US[[i]], subset = nFeature_RNA < quantile(NEG.US[[i]]$nFeature_RNA, probs=0.95) )
}
dev.off()

pdf("Hist2.pdf")
hist(SNG.US[[1]]$nFeature_RNA, breaks = 100)
for(i in 1:length(SNG.US)){
  SNG.US[[i]] = subset(SNG.US[[i]], subset = nFeature_RNA < quantile(SNG.US[[i]]$nFeature_RNA, probs=0.95) )
}
dev.off()


for(i in 1:length(SNG.US)){
	SNG.US[[i]] = subset(SNG.US[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mito < 0.30 & nCount_CITE < 20000 )
}

hist(NEG.US[[1]]$nFeature_RNA, breaks = 100)
for(i in 1:length(NEG.US)){
	NEG.US[[i]] = subset(NEG.US[[i]], subset = nFeature_RNA < quantile(NEG.US[[i]]$nFeature_RNA, probs=0.95) )
}



##############
##############
##############
############## FAILS HERE
##############
#RatIgG1kiso,RatIgG1kiso,R2,5PNNNNNNNNNN(BC),ATCAGATGCCCTCAT,Antibody Capture
#RatIgG2akiso,RatIgG2akiso,R2,5PNNNNNNNNNN(BC),AAGTCAGGTTCGTTT,Antibody Capture
#ArmenianHamsterIgGiso,ArmenianHamsterIgGiso,R2,5PNNNNNNNNNN(BC),CCTGTCATTAAGACT,Antibody Capture
#RatIgG2bkIso,RatIgG2bkIso,R2,5PNNNNNNNNNN(BC),GATTCTTGACGACCT,Antibody Capture
#IgG2bkiso,IgG2bkiso,R2,5PNNNNNNNNNN(BC),ATATGTATCACGCGA,Antibody Capture
#IgG2akiso,IgG2akiso,R2,5PNNNNNNNNNN(BC),CTCCTACCTAAACTG,Antibody Capture
#IgG1kiso,IgG1kiso,R2,5PNNNNNNNNNN(BC),GCCGGACGACATTAA,Antibody Capture

#isotype.control.name.vec = c("RatIgG1kiso", "RatIgG2akiso", "ArmenianHamsterIgGiso",  "RatIgG2bkIso", "IgG2akiso", "IgG1kiso", "IgG2bkiso")

is.array(GetAssayData(SNG.US[[1]], assay = "CITE", slot = "counts")[c(rownames(SNG.US[[i]][["CITE"]][67:70,])),])
##WORKS GetAssayData(SNG.US[[i]], assay = "CITE", slot = "counts")[c(rownames(SNG.US[[i]][["CITE"]][67:70,])),]

isotype.control.name.vec =c(rownames(SNG.US[[i]][["CITE"]][67:70,]))
#isotype.control.name.vec = c("MouseIgG1kappaisotype-PROT","MouseIgG2akappaisotype-PROT",
#"MouseIgG2bkIsotype-PROT","RatIgG2bkIsotype-PROT")

#DOES IT NEED TO BE ORDERED WHY NOT WORKED
#isotype.control.name.vec = c("MouseIgG2bkIsotype_PROT","MouseIgG2akappaisotype_PROT","MouseIgG1kappaisotype_PROT","RatIgG2bkIsotype_PROT")
##rownames(B1_US_data[[i]]$'Antibody Capture')
##index 67:70

for(i in 1:length(SNG.US)){
   SNG.US[[i]] <- AddMetaData(object = SNG.US[[i]], metadata = colMeans(as.matrix(GetAssayData(SNG.US[[i]], assay = "CITE", slot = "counts")[isotype.control.name.vec,])), col.name = "isotype.mean")
}

#for(i in 1:length(SNG.US)){
#   SNG.US[[i]] <- AddMetaData(object = SNG.US[[i]], metadata = colMeans(
#as.matrix(GetAssayData(SNG.US[[i]], assay = "CITE", slot = "counts")[67:70,]), col.name = "isotype.mean"))
#}


for(i in 1:length(NEG.US)){
  NEG.US[[i]] <- AddMetaData(object = NEG.US[[i]], metadata = colMeans(as.matrix(GetAssayData(NEG.US[[i]], assay = "CITE", slot = "counts")[isotype.control.name.vec,])), col.name = "isotype.mean")
}

#for(i in 1:length(NEG.US)){
#  NEG.US[[i]] <- AddMetaData(object = NEG.US[[i]], metadata = colMeans(GetAssayData(NEG.US[[i]], assay = "CITE", slot = "counts")[67:70,]), col#.name = "isotype.mean")
#}



for(i in 1:length(SNG.US)){
  SNG.US[[i]] = subset(SNG.US[[i]], subset = isotype.mean < quantile(SNG.US[[i]]$isotype.mean, probs=0.99) )
}
for(i in 1:length(NEG.US)){
  NEG.US[[i]] = subset(NEG.US[[i]], subset = isotype.mean < quantile(NEG.US[[i]]$isotype.mean, probs=0.99) )
}

#filter out CITEseq outliers cells

for(i in 1:length(SNG.US)){
  SNG.US[[i]] <- AddMetaData(object = SNG.US[[i]], metadata = log1p(SNG.US[[i]]$nCount_CITE), col.name = "nCount_CITE.log1p")
}
for(i in 1:length(NEG.US)){
  NEG.US[[i]] <- AddMetaData(object = NEG.US[[i]], metadata = log1p(NEG.US[[i]]$nCount_CITE), col.name = "nCount_CITE.log1p")
}

for(i in 1:length(SNG.US)){
  SNG.US[[i]] = subset(SNG.US[[i]], subset = nCount_CITE.log1p < (median(SNG.US[[i]]$nCount_CITE.log1p) + 3*mad(SNG.US[[i]]$nCount_CITE.log1p)))
  SNG.US[[i]] = subset(SNG.US[[i]], subset = nCount_CITE.log1p >  (median(SNG.US[[i]]$nCount_CITE.log1p) - 3*mad(SNG.US[[i]]$nCount_CITE.log1p)))
}
for(i in 1:length(NEG.US)){
  NEG.US[[i]] = subset(NEG.US[[i]], subset = nCount_CITE.log1p < (median(NEG.US[[i]]$nCount_CITE.log1p) + 3*mad(NEG.US[[i]]$nCount_CITE.log1p)))
  NEG.US[[i]] = subset(NEG.US[[i]], subset = nCount_CITE.log1p >  (median(NEG.US[[i]]$nCount_CITE.log1p) - 3*mad(NEG.US[[i]]$nCount_CITE.log1p)))
}





#for(i in 1:length(SNG.US)){
#  norm.dsb.temp = DSBNormalizeProtein(cell_protein_matrix = rbind(as.matrix(GetAssayData(SNG.US[[i]][["HTO"]], slot = "counts")),
#                                                                  as.matrix(GetAssayData(SNG.US[[i]][["CITE"]], slot = "counts"))),
#                                      empty_drop_matrix = rbind(as.matrix(GetAssayData(NEG.US[[i]][["HTO"]], slot = "counts")),
#                                                                as.matrix(GetAssayData(NEG.US[[i]][["CITE"]], slot = "counts"))),
#                                      define.pseudocount = TRUE, pseudocount.use = 10, denoise.counts = TRUE,
#                                      isotype.control.name.vec = isotype.control.name.vec)
#  SNG.US[[i]][["HTO"]] <- SetAssayData(SNG.US[[i]][["HTO"]], slot = "data", new.data = norm.dsb.temp[1:6,])
#  SNG.US[[i]][["CITE"]] <- SetAssayData(SNG.US[[i]][["CITE"]], slot = "data", new.data = norm.dsb.temp[1:113,])
#}


i=1
 norm.dsb.temp1 = DSBNormalizeProtein(cell_protein_matrix = rbind(as.matrix(GetAssayData(SNG.US[[i]][["HTO"]], slot = "counts")),
                                                                  as.matrix(GetAssayData(SNG.US[[i]][["CITE"]], slot = "counts"))),
                                      empty_drop_matrix = rbind(as.matrix(GetAssayData(NEG.US[[i]][["HTO"]], slot = "counts")),
                                                                as.matrix(GetAssayData(NEG.US[[i]][["CITE"]], slot = "counts"))),
                                      define.pseudocount = TRUE, pseudocount.use = 10, denoise.counts = TRUE,
                                      isotype.control.name.vec = isotype.control.name.vec)
  SNG.US[[i]][["HTO"]] <- SetAssayData(SNG.US[[i]][["HTO"]], slot = "data", new.data = norm.dsb.temp1[1:5,])
  SNG.US[[i]][["CITE"]] <- SetAssayData(SNG.US[[i]][["CITE"]], slot = "data", new.data = norm.dsb.temp1[6:118,])
  

i=2
 norm.dsb.temp2 = DSBNormalizeProtein(cell_protein_matrix = rbind(as.matrix(GetAssayData(SNG.US[[i]][["HTO"]], slot = "counts")),
                                                                  as.matrix(GetAssayData(SNG.US[[i]][["CITE"]], slot = "counts"))),
                                      empty_drop_matrix = rbind(as.matrix(GetAssayData(NEG.US[[i]][["HTO"]], slot = "counts")),
                                                                as.matrix(GetAssayData(NEG.US[[i]][["CITE"]], slot = "counts"))),
                                      define.pseudocount = TRUE, pseudocount.use = 10, denoise.counts = TRUE,
                                      isotype.control.name.vec = isotype.control.name.vec)
  SNG.US[[i]][["HTO"]] <- SetAssayData(SNG.US[[i]][["HTO"]], slot = "data", new.data = norm.dsb.temp2[1:5,])
  SNG.US[[i]][["CITE"]] <- SetAssayData(SNG.US[[i]][["CITE"]], slot = "data", new.data = norm.dsb.temp2[6:118,])
  

  
  i=3
  norm.dsb.temp2 = DSBNormalizeProtein(cell_protein_matrix = rbind(as.matrix(GetAssayData(SNG.US[[i]][["HTO"]], slot = "counts")),
                                                                   as.matrix(GetAssayData(SNG.US[[i]][["CITE"]], slot = "counts"))),
                                       empty_drop_matrix = rbind(as.matrix(GetAssayData(NEG.US[[i]][["HTO"]], slot = "counts")),
                                                                 as.matrix(GetAssayData(NEG.US[[i]][["CITE"]], slot = "counts"))),
                                       define.pseudocount = TRUE, pseudocount.use = 10, denoise.counts = TRUE,
                                       isotype.control.name.vec = isotype.control.name.vec)
  SNG.US[[i]][["HTO"]] <- SetAssayData(SNG.US[[i]][["HTO"]], slot = "data", new.data = norm.dsb.temp2[1:6,])
  SNG.US[[i]][["CITE"]] <- SetAssayData(SNG.US[[i]][["CITE"]], slot = "data", new.data = norm.dsb.temp2[7:119,])
  
  i=4
  norm.dsb.temp2 = DSBNormalizeProtein(cell_protein_matrix = rbind(as.matrix(GetAssayData(SNG.US[[i]][["HTO"]], slot = "counts")),
                                                                   as.matrix(GetAssayData(SNG.US[[i]][["CITE"]], slot = "counts"))),
                                       empty_drop_matrix = rbind(as.matrix(GetAssayData(NEG.US[[i]][["HTO"]], slot = "counts")),
                                                                 as.matrix(GetAssayData(NEG.US[[i]][["CITE"]], slot = "counts"))),
                                       define.pseudocount = TRUE, pseudocount.use = 10, denoise.counts = TRUE,
                                       isotype.control.name.vec = isotype.control.name.vec)
  SNG.US[[i]][["HTO"]] <- SetAssayData(SNG.US[[i]][["HTO"]], slot = "data", new.data = norm.dsb.temp2[1:6,])
  SNG.US[[i]][["CITE"]] <- SetAssayData(SNG.US[[i]][["CITE"]], slot = "data", new.data = norm.dsb.temp2[7:119,])
  
  
##NOT GOING TO DO MANUAL GATING AS ITS OPTIONAL I THINK
#manual gating of HTOs
dir.create("BatchCsvFiles")

for(i in 1:length(SNG.US)){
  temp.table = t(rbind(as.matrix(GetAssayData(SNG.US[[i]], assay = "CITE", slot = "data")), as.matrix(GetAssayData(SNG.US[[i]], assay = "HTO", slot = "data"))))
  rownames(temp.table) <- 1:length(rownames(temp.table))
  write.table(temp.table ,
              file = paste("BatchCsvFiles/US.", names(SNG.US[i]), ".csv", sep = ""), quote = FALSE, row.names=TRUE, col.names=NA, sep=",")
}

#save objects for clustering
dir.create("SeuratObjects")

saveRDS(SNG.US, file = "SeuratObjects/SNG.US.forClust.rds")









