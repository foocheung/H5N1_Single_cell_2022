#this is tested using R 4.0.3 on a high-performance comupting node with 16 cores and at least 384 gb or ram. (Cicero Conda env)
library("Seurat") #load Seurat 4
library(statmod)
library(magrittr)
library(limma)
library(edgeR)
library(tidyverse)
library(ggsci)
library("dplyr")
library("matrixStats")
library('tidyverse')
library('genefilter')
library("viridis")
library('pheatmap')
library('readr')
library(Biobase)
library(tidyseurat)
library(fgsea)
source("./Rfns/Full_DE_Workflow_MPM_AJMmod_20210330.R")

remerged.US.QCd <- readRDS("SeuratObjects/remerged.US.WCT.QCd.ForDE.Rds")


# make output directories
dir.create(path = "DE_wctcelltype_QCd", showWarnings=FALSE)
dir.create(path = "DE_wctcelltype_QCd/figure_output", showWarnings=FALSE)
dir.create(path = "DE_wctcelltype_QCd/gsea_output", showWarnings=FALSE)
dir.create(path = "DE_wctcelltype_QCd/model_output", showWarnings=FALSE)
dir.create(path = "DE_wctcelltype_QCd/module_score_out", showWarnings=FALSE)
dir.create(path = "DE_wctcelltype_QCd/color_schemes", showWarnings=FALSE)
dir.create(path = "DE_wctcelltype_QCd/object_output", showWarnings=FALSE)


###########################################################

remerged.US.QCd_orig<-remerged.US.QCd 
remerged.US.QCd <-remerged.US.QCd  %>% filter(grepl("D",Timepoint)) 

remerged.US.QCd$Batch_Donor_Timepoint <-paste(remerged.US.QCd$Batch, 
                                              remerged.US.QCd$Donor,
                                              remerged.US.QCd$Timepoint, sep="_")



remerge.US.DeltaDelta <- remerged.US.QCd
#format metadata for pipeline
Idents(remerge.US.DeltaDelta) <- "wctcelltype"

#remove low abundance samples and cell types

celltypes = unique(remerge.US.DeltaDelta$wctcelltype)
celltypes

celltypes_mytable.DeltaDelta = GetSubjectCelltypeTable(Seurat.Object = remerge.US.DeltaDelta, celltype_column = "wctcelltype", sample_column = "Batch_Donor_Timepoint")
celltypes_mytable.DeltaDelta

celltypes = c("CD16 Mono", "CD14 Mono", "CD8 TEM", "gdT", "CD4 TCM", "B memory", 
              "B naive", "B intermediate",   "CD4 Naive", 
              "CD4 CTL", "NK", "NK_CD56bright", "CD8 Naive", "CD8 TCM", "MAIT", 
              "pDC", "cDC2") ##, "cDC1")


## Create aggregate dataset 1 -- Pseudobulk counts at the sample level.
celltypes_sl.DeltaDelta = MakePseudobulkList(seurat.object = remerge.US.DeltaDelta, celltype_column = "wctcelltype",
                                             sample_column = "Batch_Donor_Timepoint", vector.of.celltypes = celltypes)


## Set up model and contrast matrix


names(celltypes_sl.DeltaDelta)
celltypes_sl.DeltaDelta[[5]]$RNA[1:5,]

remerge.US.DeltaDelta$Subject<-gsub( "-","_", remerge.US.DeltaDelta$Subject)
remerge.US.DeltaDelta$Group <-paste(remerge.US.DeltaDelta$Group, 
                                        remerge.US.DeltaDelta$Timepoint, sep="_")

MakeMetaTableFromSeurat(seurat_object = remerge.US.DeltaDelta, sample_column = "Batch_Donor_Timepoint",
                        variable_column = c("Subject"),                              
                        aggregate.data.list = celltypes_sl.DeltaDelta,
                        celltypes.vector = celltypes)

celltypes_sxmd.DeltaDelta = MakeMetaTableFromSeurat(seurat_object = remerge.US.DeltaDelta, sample_column = "Batch_Donor_Timepoint",
                                                    variable_column = "Group",                              
                                                    aggregate.data.list = celltypes_sl.DeltaDelta,
                                                    celltypes.vector = celltypes)
celltypes_sxmd.DeltaDelta
celltypes_sxmd.DeltaDeltaSubject = MakeMetaTableFromSeurat(seurat_object = remerge.US.DeltaDelta, sample_column = "Batch_Donor_Timepoint",
                                                           variable_column = "Subject",                              
                                                           aggregate.data.list = celltypes_sl.DeltaDelta,
                                                           celltypes.vector = celltypes)
celltypes_sxmd.DeltaDeltaSubject

celltypes_sxmd.DeltaDeltaBatch = MakeMetaTableFromSeurat(seurat_object = remerge.US.DeltaDelta, sample_column = "Batch_Donor_Timepoint",
                                                         variable_column = "Batch",                              
                                                         aggregate.data.list = celltypes_sl.DeltaDelta,
                                                         celltypes.vector = celltypes)
celltypes_sxmd.DeltaDeltaBatch

celltypes_sxmd.DeltaDeltaTime = MakeMetaTableFromSeurat(seurat_object = remerge.US.DeltaDelta, sample_column = "Batch_Donor_Timepoint",
                                                        variable_column = "Timepoint",                              
                                                        aggregate.data.list = celltypes_sl.DeltaDelta,
                                                        celltypes.vector = celltypes)
celltypes_sxmd.DeltaDeltaTime


###MUST RUN PREVIOUS SECTION AGAIN TO GET THIS WORKING
###

for(i in 1:length(celltypes_sxmd.DeltaDelta)){
  rownames(celltypes_sxmd.DeltaDelta[[i]]) <- NULL
  rownames(celltypes_sxmd.DeltaDeltaBatch[[i]]) <- NULL
  rownames(celltypes_sxmd.DeltaDeltaSubject[[i]]) <- NULL
  celltypes_sxmd.DeltaDelta[[i]] =
    celltypes_sxmd.DeltaDelta[[i]] %>%
    select(sample,Group) %>%
    column_to_rownames("sample")
  
  
  celltypes_sxmd.DeltaDeltaBatch[[i]] =
    celltypes_sxmd.DeltaDeltaBatch[[i]] %>%
    select(sample,Batch) %>%
    column_to_rownames("sample")	
  celltypes_sxmd.DeltaDelta[[i]]$Batch <- celltypes_sxmd.DeltaDeltaBatch[[i]][rownames(celltypes_sxmd.DeltaDelta[[i]]), 1]
  
  celltypes_sxmd.DeltaDeltaSubject[[i]] =
    celltypes_sxmd.DeltaDeltaSubject[[i]] %>%
    select(sample,Subject) %>%
    column_to_rownames("sample")	
  celltypes_sxmd.DeltaDelta[[i]]$Subject <- celltypes_sxmd.DeltaDeltaSubject[[i]][rownames(celltypes_sxmd.DeltaDelta[[i]]), 1]
  
  
}



form.DeltaDelta <- ~ 0 + Group + Batch  + Subject
#form.DeltaDelta <- ~ 0 + Group + Batch # + Subject
celltypes_design.DeltaDelta <-list()
for(i in 1:length(celltypes_sxmd.DeltaDelta)){
  celltypes_design.DeltaDelta[[i]]<- model.matrix(form.DeltaDelta, celltypes_sxmd.DeltaDelta[[i]])
}


cont.matrix.DeltaDelta <-list()
for(i in 1:length(celltypes_sxmd.DeltaDelta)){
  
  cont.matrix.DeltaDelta[[i]] <- makeContrasts( 
    
    
    g1=GroupNEG_D22-GroupNEG_D21,
    g2=GroupNEG_D21_2h-GroupNEG_D21,
    g3=GroupNEG_D21_4h-GroupNEG_D21,
    g4=GroupNEG_D21_12h-GroupNEG_D21,
    
    g5=GroupPOS_D22-GroupPOS_D21,  
    g6=GroupPOS_D21_2h-GroupPOS_D21,  
    g7=GroupPOS_D21_4h-GroupPOS_D21,  
    g8=GroupPOS_D21_12h-GroupPOS_D21,  
    
    
    #   g3=(GroupPOS_D1-GroupPOS_D0) - (GroupNEG_D1-GroupNEG_D0),
    
    #   g4=GroupPOS_D100-GroupPOS_D21,  
    
    #    g5=GroupNEG_D100-GroupNEG_D21,
    
    #  g6=GroupPOS_D100-GroupPOS_D21,  
    
    #  g7=GroupNEG_D100-GroupNEG_D21,
    
    
    levels=celltypes_design.DeltaDelta[[i]])
}

write.csv(prop.table(table(remerged.US.QCd$Batch_Donor_Timepoint, remerged.US.QCd$wctcelltype), margin=1), file = "WCTcelltypeFreqs.csv")

celltypes_freq <- read.csv(file = "WCTcelltypeFreqs.csv", row.names=1)

celltypes_freq.eset <- Biobase::ExpressionSet(assayData = as.matrix(celltypes_freq))
celltypes_freq.fit <- lmFit(t(Biobase::exprs(celltypes_freq.eset))[,rownames(celltypes_design.DeltaDelta[[1]])], 
                            design= celltypes_design.DeltaDelta[[1]], 
                            contrasts = cont.matrix.DeltaDelta[[1]] )

celltypes_freq.ebfit <- eBayes(celltypes_freq.fit)



#topTable(celltypes_freq.ebfit, coef=1, number = 40) #"preVsPostInAsym" contrast
#topTable(celltypes_freq.ebfit, coef=2, number = 40) #"preVsPostInSym" contrast
#topTable(celltypes_freq.ebfit, coef=3, number = 40) #"Diff" contrast


# ### New function "Get Sample Weights"
GetSampleWeights = function(subject.celltype.table){
  
  tab = subject.celltype.table
  celltype_sums = rowSums(tab)
  weights = apply(tab, 2, function(x){ x / celltype_sums  } )
  weights = weights %>% t %>% as.data.frame()
  return(weights)
}

celltypes_cell_weights.DeltaDelta = GetSampleWeights(subject.celltype.table = celltypes_mytable.DeltaDelta$table[celltypes,])
celltypes_cell_weights.DeltaDelta






celltypes_genes.use1.DeltaDelta = list()
celltypes_dflist.DeltaDelta = list()
celltypes_efitlist.DeltaDelta = list()

#
#for (i in c(1:8,10:length(celltypes_sl.DeltaDelta))) {
for (i in 1:length(celltypes_sl.DeltaDelta)) {
  celltypes_dflist.DeltaDelta[[i]] = DGEList(celltypes_sl.DeltaDelta[[i]]$RNA)
  celltypes_dflist.DeltaDelta[[i]] = calcNormFactors(celltypes_dflist.DeltaDelta[[i]], method = "RLE")
  celltypes_genes.use1.DeltaDelta[[i]] = filterByExpr(celltypes_dflist.DeltaDelta[[i]], min.count = 1, min.prop=0.1)
  celltypes_dflist.DeltaDelta[[i]] = celltypes_dflist.DeltaDelta[[i]][celltypes_genes.use1.DeltaDelta[[i]], keep.lib.sizes=FALSE]
}
lapply(celltypes_dflist.DeltaDelta, FUN = function(x) print(dim(x)))

celltypes_v.DeltaDelta <- list()
celltypes_fit.DeltaDelta <- list()
celltypes_lcpm.DeltaDelta <- list()
celltypes_efit.DeltaDelta <- list()
celltypes_cfit.DeltaDelta <- list()

for (i in 1:length(celltypes_sl.DeltaDelta)) {
  #for (i in c(1:8,10:length(celltypes_sl.DeltaDelta))) {
  celltypes_v.DeltaDelta[[i]] <- voom(celltypes_dflist.DeltaDelta[[i]], celltypes_design.DeltaDelta[[i]], plot=FALSE)
  celltypes_cellweight_i.DeltaDelta = celltypes_cell_weights.DeltaDelta[colnames(celltypes_dflist.DeltaDelta[[i]]), i]
  celltypes_new_weight.DeltaDelta = apply(celltypes_v.DeltaDelta[[i]]$weights, 1, function(x) x * celltypes_cellweight_i.DeltaDelta) %>% t
  celltypes_fit.DeltaDelta[[i]] <- lmFit(celltypes_v.DeltaDelta[[i]], celltypes_design.DeltaDelta[[i]], weights = celltypes_new_weight.DeltaDelta)
  celltypes_cfit.DeltaDelta[[i]] = contrasts.fit(fit = celltypes_fit.DeltaDelta[[i]], contrasts = cont.matrix.DeltaDelta[[i]])
  celltypes_efit.DeltaDelta[[i]] <- eBayes(celltypes_cfit.DeltaDelta[[i]])
  celltypes_lcpm.DeltaDelta[[i]] <- cpm(celltypes_dflist.DeltaDelta[[i]], log=TRUE)
  
}

names(celltypes_v.DeltaDelta) <- names(celltypes_fit.DeltaDelta) <- names(celltypes_lcpm.DeltaDelta) <- names(celltypes_efit.DeltaDelta) <- names(celltypes_cfit.DeltaDelta) <- celltypes



DeltaDelta_1 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta, 
                              coefficient.number = 1, contrast.name = "Diff")


DeltaDelta_2 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta, 
                              coefficient.number = 2, contrast.name = "Diff")

DeltaDelta_3 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta, 
                              coefficient.number = 3, contrast.name = "Diff")


DeltaDelta_4 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta, 
                              coefficient.number = 4, contrast.name = "Diff")
DeltaDelta_5 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta, 
                              coefficient.number = 5, contrast.name = "Diff")

DeltaDelta_6 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta, 
                              coefficient.number = 6, contrast.name = "Diff")

DeltaDelta_7 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta, 
                              coefficient.number = 7, contrast.name = "Diff")


DeltaDelta_8 = GetRankResults(limma.fit.object.list = celltypes_efit.DeltaDelta,
                              coefficient.number = 8, contrast.name = "Diff")


btm = gmtPathways("./Genesets/patterns.gmt")


f1<-RunFgseaOnRankList(rank.list.celltype =                                                                                     
                         DeltaDelta_1, 
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)

f2<-RunFgseaOnRankList(rank.list.celltype =                                                                                     
                         DeltaDelta_2, 
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)

f3<-RunFgseaOnRankList(rank.list.celltype =                                                                                     
                         DeltaDelta_3, 
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)

f4<-RunFgseaOnRankList(rank.list.celltype =                                                                                     
                         DeltaDelta_4, 
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)

f5<-RunFgseaOnRankList(rank.list.celltype =                                                                                     
                         DeltaDelta_5, 
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)

f6<-RunFgseaOnRankList(rank.list.celltype =                                                                                     
                         DeltaDelta_6, 
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)

f7<-RunFgseaOnRankList(rank.list.celltype =                                                                                     
                         DeltaDelta_7, 
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)

f8<-RunFgseaOnRankList(rank.list.celltype =
                         DeltaDelta_8,
                       pathways = btm, positive.enrich.only = FALSE, celltypes.vector = celltypes)


celltypes_score.DeltaDelta1 = lapply(f1, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
 # filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})

celltypes_score.DeltaDelta2 = lapply(f2, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
 # filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})
celltypes_score.DeltaDelta3 = lapply(f3, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
#  filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})


celltypes_score.DeltaDelta4 = lapply(f4, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
 # filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})


celltypes_score.DeltaDelta5 = lapply(f5, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
#  filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})

celltypes_score.DeltaDelta6 = lapply(f6, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
#  filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})


celltypes_score.DeltaDelta7 = lapply(f7, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
#  filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})



celltypes_score.DeltaDelta8 = lapply(f8, function(x){ x = x %>%
  select(pathway, padj, NES,ES, leadingEdge,celltype) %>% 
#  filter(padj < 0.1)  %>%
  mutate(n_logp = -log10(padj))
})

system("mkdir DE_wctcelltype_QCd/")
system("mkdir DE_wctcelltype_QCd/figure_output")
system("mkdir DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta")

##################################3
#merge the list into a single data frame
do.call(rbind, celltypes_score.DeltaDelta1) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_neg_adj_D22-D21.tsv')


#########################################
do.call(rbind, celltypes_score.DeltaDelta2) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_neg_adj_D21_2h-D21.tsv')


#########################################

do.call(rbind, celltypes_score.DeltaDelta3) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_neg_adj_D21_4h-D21.tsv')

#########################################

do.call(rbind, celltypes_score.DeltaDelta4) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_neg_adj_D21_12h-D21.tsv')


#########################################

do.call(rbind, celltypes_score.DeltaDelta5) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_pos_adj_D22-D21.tsv')

#########################################

do.call(rbind, celltypes_score.DeltaDelta6) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_pos_adj_D21_2h-D21.tsv')

#########################################

do.call(rbind, celltypes_score.DeltaDelta7) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_pos_adj_D21_4h-D21.tsv')


#########################################

do.call(rbind, celltypes_score.DeltaDelta8) %>% 
  rowwise() %>% 
  mutate_if(is.list, ~paste(unlist(.), collapse = ',')) %>% 
  write_tsv('DE_wctcelltype_QCd/figure_output/module_heatmaps.DeltaDelta/le_patterns_pos_adj_D21_12h-D21.tsv')

