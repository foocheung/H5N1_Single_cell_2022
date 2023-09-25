library("Seurat") #load Seurat 4
library("dplyr")
library("matrixStats")
library('tidyverse')
library('genefilter')
library("viridis")
library('pheatmap')
library('parallelDist')

SNG.US <- readRDS(file = "SeuratObjects/SNG.US.ADTclust.rds")


dir.create("ClusterVis")

for(i in 1:length(SNG.US)){
  Idents(SNG.US[[i]]) <- "adtfilt_snn_res.1"
}  

#ADTmarkers of predictedcelltypes
ADTMrks <- list()
for(i in 1:length(SNG.US)){
  Idents(SNG.US[[i]]) <- "predicted.celltype.l2"
  ADTMrks[[i]] <- FindAllMarkers(SNG.US[[i]], assay = "CITE", only.pos = TRUE,
                                 max.cells.per.ident = 1000)
}

celltypeADTMrks <- list()
for(i in 1:length(SNG.US)){
  celltypeADTMrks[[i]] = ADTMrks[[i]] %>% group_by(cluster) 
}
##celltypeADTMrksSelected = intersect(celltypeADTMrks[[1]]$gene, celltypeADTMrks[[2]]$gene)
celltypeADTMrksSelected = intersect(celltypeADTMrks[[1]]$gene, intersect(celltypeADTMrks[[2]]$gene, intersect(celltypeADTMrks[[3]]$gene,celltypeADTMrks[[4]]$gene)))

for(i in 1:length(SNG.US)){
  SNG.US[[i]] <- RunUMAP(SNG.US[[i]], assay = "CITE", features = rownames(SNG.US[[i]][["CITE"]][which(rownames(SNG.US[[i]][["CITE"]]) %in% celltypeADTMrksSelected)]), n_neighbors=50L, min_dist=1)
  ggsave(paste0("ClusterVis/", names(SNG.US[i]), ".adtfilt.umap.png"), width = 5, height = 5,
         plot = AugmentPlot(DimPlot(SNG.US[[i]], pt.size = 0.5, reduction = "umap", group.by = "adtfilt_snn_res.1", label=TRUE)) + NoLegend()
  )
}


quantile_breaks <- function(xs, n = 100) {
  breaks <- quantile(xs, probs = seq(0.5, 1, length.out = n))
  breaks[!duplicated(breaks)]
}




for(i in 1:length(SNG.US)){
  SNG.US[[i]]$adtfilt_snn_res.1.char <- as.character(SNG.US[[i]]$adtfilt_snn_res.1) 
  Idents(SNG.US[[i]]) <- "adtfilt_snn_res.1.char"
  SNG.US[[i]] <- SetAssayData(SNG.US[[i]], assay = "CITE", slot = "scale.data", new.data = GetAssayData(SNG.US[[i]], assay = "CITE", slot = "data"))
  aver.temp = AverageExpression(SNG.US[[i]], assays = "CITE", slot = "scale.data", return.seurat=T)
  mat_breaks <- quantile_breaks(GetAssayData(aver.temp, assay = "CITE", slot = "scale.data"), n = 101)
  protstoplot = rownames(SNG.US[[i]][["CITE"]][which(rownames(SNG.US[[i]][["CITE"]]) %in% celltypeADTMrksSelected)])
  
  ggsave(filename = paste0("ClusterVis/", names(SNG.US[i]),".adtfilt.averExp_res.1_CITE.pdf"), width = 10, height = 18,
         plot = pheatmap(GetAssayData(aver.temp, assay = "CITE", slot = "scale.data")[protstoplot,], scale = "none", border_color=NA, inferno(length(mat_breaks) - 1), breaks = mat_breaks))
  
}



## visualization
for(i in 1:length(SNG.US)){
  ggsave(filename = paste0("ClusterVis/", names(SNG.US[i]),"Umap_Unsorted.SeuratPredictedFromCITEref.pdf") , 
         AugmentPlot(DimPlot(SNG.US[[i]], reduction = "umap", group.by = "predicted.celltype.l2", label = TRUE)) + NoLegend(), width = 15, height = 10)
  write.csv(table(SNG.US[[i]]$predicted.celltype.l1, Idents(SNG.US[[i]])), paste0("ClusterVis/", names(SNG.US[i]),"UnsortedADTclustvsPredictedFromSeuratCITERefl1.csv"))
  write.csv(table(SNG.US[[i]]$predicted.celltype.l2, Idents(SNG.US[[i]])), paste0("ClusterVis/", names(SNG.US[i]),"UnsortedADTvsPredictedFromSeuratCITERefl2.csv"))
  pheatmap(prop.table(x = table(SNG.US[[i]]$predicted.celltype.l1, Idents(SNG.US[[i]])), margin =  2)  %>% '*'(100) %>% round(2),
           filename = paste0("ClusterVis/", names(SNG.US[i]),".UnsortedADTclustvsPredictedFromSeuratCITERefl1.pdf"), width = 15, height = 10,
           cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE)
  pheatmap(prop.table(x = table(SNG.US[[i]]$predicted.celltype.l2, Idents(SNG.US[[i]])), margin =  2)  %>% '*'(100) %>% round(2),
           filename = paste0("ClusterVis/", names(SNG.US[i]),".UnsortedADTclustvsPredictedFromSeuratCITERefl2.pdf"), width = 15, height = 10,
           cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE)
  #  pheatmap(prop.table(x = table(SNG.US[[i]]$ManualGatePop, Idents(SNG.US[[i]])), margin =  2)  %>% '*'(100) %>% round(2),
  #         filename = paste0("ClusterVis/", names(SNG.US[i]),".UnsortedADTclustvsManualGatePop.pdf"), width = 30, height = 15,
  #         cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE)
  
}





##TO DO





#name the clusters. bring back in the clusternames here from csv files
B1.ADTclustNames = read.csv("SampleMeta/B1clustername.csv", header=TRUE)
B2.ADTclustNames = read.csv("SampleMeta/B2clustername.csv", header=TRUE)
B3.ADTclustNames = read.csv("SampleMeta/B3clustername.csv", header=TRUE)
B4.ADTclustNames = read.csv("SampleMeta/B4clustername.csv", header=TRUE)
#B5.ADTclustNames = read.csv("SampleMeta/B5.ADTclustNames.csv", header=TRUE)

SNG.US[[1]]$celltype <- plyr::mapvalues(x = SNG.US[[1]]$adtfilt_snn_res.1, from = B1.ADTclustNames$ADTclust, to = as.character(B1.ADTclustNames$Celltype))
SNG.US[[1]]$coursecelltype <- plyr::mapvalues(x = SNG.US[[1]]$adtfilt_snn_res.1, from = B1.ADTclustNames$ADTclust, to = as.character(B1.ADTclustNames$CourseCelltype))
SNG.US[[2]]$celltype <- plyr::mapvalues(x = SNG.US[[2]]$adtfilt_snn_res.1, from = B2.ADTclustNames$ADTclust, to = as.character(B2.ADTclustNames$Celltype))
SNG.US[[2]]$coursecelltype <- plyr::mapvalues(x = SNG.US[[2]]$adtfilt_snn_res.1, from = B2.ADTclustNames$ADTclust, to = as.character(B2.ADTclustNames$CourseCelltype))
SNG.US[[3]]$celltype <- plyr::mapvalues(x = SNG.US[[3]]$adtfilt_snn_res.1, from = B3.ADTclustNames$ADTclust, to = as.character(B3.ADTclustNames$Celltype))
SNG.US[[3]]$coursecelltype <- plyr::mapvalues(x = SNG.US[[3]]$adtfilt_snn_res.1, from = B3.ADTclustNames$ADTclust, to = as.character(B3.ADTclustNames$CourseCelltype))
SNG.US[[4]]$celltype <- plyr::mapvalues(x = SNG.US[[4]]$adtfilt_snn_res.1, from = B4.ADTclustNames$ADTclust, to = as.character(B4.ADTclustNames$Celltype))
SNG.US[[4]]$coursecelltype <- plyr::mapvalues(x = SNG.US[[4]]$adtfilt_snn_res.1, from = B4.ADTclustNames$ADTclust, to = as.character(B4.ADTclustNames$CourseCelltype))
#SNG.US[[5]]$celltype <- plyr::mapvalues(x = SNG.US[[5]]$adtfilt_snn_res.1, from = B5.ADTclustNames$ADTclust, to = as.character(B5.ADTclustNames$Celltype))
#SNG.US[[5]]$coursecelltype <- plyr::mapvalues(x = SNG.US[[5]]$adtfilt_snn_res.1, from = B5.ADTclustNames$ADTclust, to = as.character(B5.ADTclustNames$CourseCelltype))

merge.US <- merge(SNG.US[[1]], SNG.US[2:length(SNG.US)])




##testSNG <-merge(SNG.USorig$B1[1:24132],SNG.USorig$B2[1:23756], add.cell.ids = c("batch1", "batch2"), project = "H5Merge")

#add metadata

#merge.US$Donor[which(is.na(merge.US$Donor))] = "CHI014"
#merge.US$Sample = paste(merge.US$Batch, merge.US$Donor,
#                         merge.US$GatedHashCall, sep = "_")
merge.US$Donor = sapply(strsplit(as.character(merge.US$BEST.GUESS),split = ","),'[',1)
#merge.US$Donor = sapply(strsplit(as.character(merge.US$Donor),split = "_"),'[',3)
merge.US$Sample<-paste(merge.US$Batch,"_",merge.US$Donor,"_", merge.US$HTO_classification, sep = "")

samplemetadata.US = read.csv("SampleMeta/h5cite_sheet_sample_sheet.csv", header = TRUE, colClasses = c(rep("character", 10)))
merge.US$Batch_Sample = merge.US$Sample
Idents(merge.US) <- "Sample"
merge.US$Subject <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Subject)
## The following `from` values were not present in `x`: 
###31_9980966277_R07C01_HTO-6, 31_9980966277_R07C01_HTO-7, 31_9980966277_R07C01_HTO-8, 31_9980966277_R07C01_HTO-9, 31_9980966277_R07C01_HTO-10, 31_9980966277_R07C01_HTO-12, 31_9980966277_R07C01_HTO-13
merge.US$Subject <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Subject)
merge.US$Age <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Age)
merge.US$Gender <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Sex)
merge.US$Group <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Group)
merge.US$Adjuvant <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Adjuvant)
#merge.US$MatchedPair <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$CITEseq.pairs.Matched)
merge.US$Timepoint <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$TimePoint)
#merge.US$X1st.PCR..aka.d1..for.reference. <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$X1st.PCR..aka.d1..for.reference.)
#merge.US$ra1 <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$ra1)
#merge.US$eth <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$eth)
#merge.US$Time1 <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Time1)
#merge.US$Time1_Sgene <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Time1_Sgene)
#merge.US$Time1_Ngene <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Time1_Ngene)
#merge.US$Time1_ORF1ab <- plyr::mapvalues(x = Idents(merge.US), from = samplemetadata.US$Sample, to = samplemetadata.US$Time1_ORF1ab)
merge.US$Donor.Timepoint = paste(merge.US$Donor, merge.US$Timepoint, sep="_")

# table(merge.US$Donor.Timepoint, merge.US$SelectionGroup) 
# table(merge.US$Donor.Timepoint, merge.US$Age) 
# table(merge.US$Donor.Timepoint, merge.US$Gender) 
# table(merge.US$Donor.Timepoint, merge.US$Timepoint)
# table(merge.US$Donor.Timepoint, merge.US$Symptomatic.or.asymptomatic ) 
# table(merge.US$Batch, merge.US$celltype ) 
# table(merge.US$Batch, merge.US$coursecelltype ) 
# table(subset(merge.US, Donor == "CHI014")$Batch, subset(merge.US, Donor == "CHI014")$celltype ) 


merge.US <- DietSeurat(merge.US, scale.data=TRUE)

saveRDS(merge.US, file = "SeuratObjects/merge.US.rds")


