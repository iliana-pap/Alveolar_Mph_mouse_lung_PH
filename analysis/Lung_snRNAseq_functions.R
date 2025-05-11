
preprocessing <- function(seurat_obj, resolution){
  
  seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(object = seurat_obj, npcs = 30, verbose = FALSE,seed.use = 8734)
  seurat_obj <- RunTSNE(object = seurat_obj, reduction = "pca", dims = 1:20, seed.use = 8734, check_duplicates = FALSE)
  seurat_obj <- FindNeighbors(object = seurat_obj, reduction = "pca", dims = 1:20)
  for(k in 1:length(resolution)){
    seurat_obj <- FindClusters(object = seurat_obj, resolution = resolution[k], random.seed = 8734)
  }
  seurat_obj <- RunUMAP(object = seurat_obj, reduction = "pca", dims = 1:20, seed.use = 8734)
  return(seurat_obj)
}


avgHeatmapTotalHeatmap <- function(seurat, selGenes, colVecIdent, colVecCond=NULL,
                                   ordVec=NULL, gapVecR=NULL, gapVecC=NULL,cc=F,
                                   cr=FALSE, condCol=FALSE,path, analysis_type, len){
  
  selGenes <- selGenes$gene
  ## assay data
  clusterAssigned <- as.data.frame(Idents(seurat)) %>%
    dplyr::mutate(cell=rownames(.))
  colnames(clusterAssigned)[1] <- "ident"
  seuratDat <- GetAssayData(seurat)
  
  ## genes of interest
  genes <- data.frame(gene=rownames(seurat)) %>% 
    mutate(geneID=gsub("^.*?\\.", "", gene)) %>% filter(geneID %in% selGenes)
  
  ## matrix with averaged cnts per ident
  logNormExpres <- as.data.frame(t(as.matrix(
    seuratDat[which(rownames(seuratDat) %in% genes$gene),])))
  
  logNormExpres <- logNormExpres %>% dplyr::mutate(cell=rownames(.)) %>%
    dplyr::left_join(.,clusterAssigned, by=c("cell")) %>%
    dplyr::select(-cell) %>% dplyr::group_by(ident)  %>%
    dplyr::summarise_all(mean)
  logNormExpresMa <- logNormExpres %>% dplyr::select(-ident) %>% as.matrix()
  rownames(logNormExpresMa) <- logNormExpres$ident
  logNormExpresMa <- t(logNormExpresMa)
  rownames(logNormExpresMa) <- gsub("^.*?\\.","",rownames(logNormExpresMa))
  
  ## remove genes if they are all the same in all groups
  ind <- apply(logNormExpresMa, 1, sd) == 0
  logNormExpresMa <- logNormExpresMa[!ind,]
  genes <- genes[!ind,]
  
  ## color columns according to cluster
  annotation_col <- as.data.frame(gsub("(^.*?_)","",
                                       colnames(logNormExpresMa)))%>%
    dplyr::mutate(label=gsub("(_.*$)","",colnames(logNormExpresMa)))
  colnames(annotation_col)[1] <- "col1"
  annotation_col <- annotation_col %>%
    dplyr::mutate(cond = gsub("(^[0-9]_?)","",col1)) %>%
    dplyr::select(cond, label)
  rownames(annotation_col) <- colnames(logNormExpresMa) 
  
  cell_order <- names(colVecIdent)
  annotation_col <- annotation_col %>% arrange(factor(label, levels = cell_order ))
  
  ann_colors = list(
    cond = colVecCond,
    label=colVecIdent)
  
  if(is.null(ann_colors$cond)){
    annotation_col$cond <- NULL
  }
  
  ## adjust order
  logNormExpresMa <- logNormExpresMa[selGenes,]
  if(is.null(ordVec)){
    ordVec <- levels(seurat)
  }

  # Reorder ordVec
  ordVec <- factor(ordVec, levels = cell_order)
  ordVec <- as.character(sort(ordVec))
  
  logNormExpresMa <- logNormExpresMa[,ordVec]
  
  if (len == 1){
    logNormExpresMa <-logNormExpresMa[selGenes[1],, drop=FALSE]
  }
  p<- pheatmap(logNormExpresMa, scale="row" ,treeheight_row = 0, cluster_rows = cr, 
               cluster_cols = cc,
               color = colorRampPalette(c("#2166AC" , "white","#C16622FF"))(50),
               annotation_col = annotation_col, cellwidth=14.5, cellheight=15,
               annotation_colors = ann_colors, gaps_row = gapVecR, gaps_col = gapVecC,
               treeheight_col = 0, border_color = FALSE)
  print(p)
}

avgHeatmap2 <- function(Selgenes,logNormExpresMa1, logNormExpresMa2,colVecIdent, colVecCond=dep_vec,gapVecR=NULL, gapVecC=NULL,cc=FALSE,
                        cr=FALSE, condCol=FALSE, type = type, folder_name = folder_name){
  ann_colors = list(
    condition = colVecCond,
    celltype=colVecIdent)
  
  celltype_label <-names(colVecIdent)
  colnames(logNormExpresMa1) <- paste(colnames(logNormExpresMa1), "Control", sep = "_")
  colnames(logNormExpresMa2) <- paste(colnames(logNormExpresMa2), "PH-cHFpEF", sep = "_")
  
  # If there is a gene that is included in one condition and not in another add new row with zero where it is missing
  if ((dim(logNormExpresMa1)[1] != dim(logNormExpresMa2)[1])) {
    gene_to_be_added1 <- setdiff(rownames(logNormExpresMa2),rownames(logNormExpresMa1))
    rn <-rownames(logNormExpresMa1)
    logNormExpresMa1<- rbind(logNormExpresMa1, matrix(0,length(gene_to_be_added1) ,length(colnames(logNormExpresMa1))))
    rownames(logNormExpresMa1) <- c(rn,c(gene_to_be_added1))
    logNormExpresMa1 <- logNormExpresMa1[match(rownames(logNormExpresMa2),rownames(logNormExpresMa1)), ]
  }
  
  logNormExpresMa <- cbind(logNormExpresMa1,logNormExpresMa2)
  
  annotation_col_updated = data.frame(
    celltype = factor(rep(celltype_label , rep(2, length(celltype_label)))),
    condition = rep(c('Control', 'PH-cHFpEF') , length(celltype_label))
  )
  
  annot_df <- annotation_col_updated
  annot_df$name_col <- paste0(annot_df$celltype, '_',annot_df$condition)
  logNormExpresMa_debug <- logNormExpresMa
  
  # Reorder matrix columns based on the annotation
  logNormExpresMa<-logNormExpresMa[,annot_df$name_col]
  
  rownames(annotation_col_updated) <- colnames(logNormExpresMa)
  
  p <- pheatmap(logNormExpresMa, scale="row" ,treeheight_row = 0, cluster_rows = cr, treeheight_col = 0,
                cluster_cols = F,
                color = colorRampPalette(c(c("#749B58FF" , "white","#67001F")))(50),
                annotation_col = annotation_col_updated, cellwidth=15, cellheight=10,
                annotation_colors = ann_colors, gaps_row = gapVecR, gaps_col = gapVecC,show_colnames=FALSE)
  p
}

compute_avg_gene_expression_matrix <- function(seurat, selGenes,ordVec = NULL){
  
  genes <- data.frame(gene=rownames(seurat)) %>%
    mutate(geneID=gsub("^.*?\\.", "", gene))
  
  Selgenes_df <- filter(genes, geneID %in% Selgenes)
  colnames(Selgenes_df) <- c('geneID','gene')
  Selgenes_df <-Selgenes_df[match(Selgenes, Selgenes_df$gene),]
  Idents(seurat) <- seurat$label
  ## assay data
  clusterAssigned <- as.data.frame(Idents(seurat)) %>%
    dplyr::mutate(cell=rownames(.))
  colnames(clusterAssigned)[1] <- "ident"
  seuratDat <- GetAssayData(seurat)
  genes <- data.frame(gene=rownames(seurat)) %>%
    mutate(geneID=gsub("^.*\\.", "", gene)) %>% filter(geneID %in% selGenes)
  
  ## matrix with averaged cnts per ident
  logNormExpres <- as.data.frame(t(as.matrix(
    seuratDat[which(rownames(seuratDat) %in% genes$gene),])))
  logNormExpres <- logNormExpres %>% dplyr::mutate(cell=rownames(.)) %>%
    dplyr::left_join(.,clusterAssigned, by=c("cell")) %>%
    dplyr::select(-cell) %>% dplyr::group_by(ident) %>%
    dplyr::summarise_all(mean)
  logNormExpresMa <- logNormExpres %>% dplyr::select(-ident) %>% as.matrix()
  rownames(logNormExpresMa) <- logNormExpres$ident
  logNormExpresMa <- t(logNormExpresMa)
  rownames(logNormExpresMa) <- gsub("^.*?\\.","",rownames(logNormExpresMa))
  ## remove genes if they are all the same in all groups
  ind <- apply(logNormExpresMa, 1, sd) == 0
  logNormExpresMa <- logNormExpresMa[!ind,]
  genes <- genes[!ind,]
  gene_to_be_added <-setdiff(selGenes, rownames(logNormExpresMa))
  rn <-rownames(logNormExpresMa)
  logNormExpresMa <- rbind(logNormExpresMa, matrix(0,length(gene_to_be_added) ,length(colnames(logNormExpresMa))))
  rownames(logNormExpresMa) <- c(rn,c(gene_to_be_added))
  logNormExpresMa <- logNormExpresMa[rownames(logNormExpresMa), ]
  
  ## adjust order
  logNormExpresMa <- logNormExpresMa[selGenes,]
  if(is.null(ordVec)){
    ordVec <- levels(seurat)
  }
  return(logNormExpresMa)
}



# Functions for pathway analysis
customize_parameters <- function(Vec,DEmarkers,organism,datatype,disease_phase, output_path){
  
  # Initialize dictionary
  enrichcl_dict <- dict()
  
  # Contains aggregated enrichment analysis S4 objects after all subanalyses
  enrichcl_list <- list()
  
  analysis_res = Enrichment_Analysis(Vec[1],DEmarkers, "GO",datatype, organism,enrichcl_dict,enrichcl_list,disease_phase,output_path)
  analysis_res = Enrichment_Analysis(Vec[2],DEmarkers, "GO",datatype, organism, analysis_res$enrichcl_dict,analysis_res$enrichcl_list,disease_phase,output_path)
  
  param <- list("qscore_min" = min(analysis_res$enrichcl_dict$qscore_min),
                "qscore_max" = max(analysis_res$enrichcl_dict$qscore_max),
                "GeneRatioNum_min" = min(analysis_res$enrichcl_dict$GeneRatioNum_min),
                "GeneRatioNum_max" = max(analysis_res$enrichcl_dict$GeneRatioNum_max),
                "Count_min" = min(analysis_res$enrichcl_dict$Count_min),
                "Count_max" = max(analysis_res$enrichcl_dict$Count_max))
  
  return(list("param" = param, "enrichcl_list" = analysis_res$enrichcl_list))
  
}

MapGenetoENSEMBL<- function(DEGmarkers, organismDB,output_path){
  
  mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
  G_list<- getBM(attributes=c("ensembl_gene_id", "mgi_symbol"),
                 filters = 'mgi_symbol',
                 values=DEGmarkers$gene,
                 mart=mart, uniqueRows=T)
  
  
  DEGmarkers <- left_join(DEGmarkers, G_list,
                          by= c("gene"="mgi_symbol"))
  
  umappped_genes <- filter(DEGmarkers, is.na(ensembl_gene_id))$gene
  DEGmarkers[DEGmarkers$gene == "2900097C17Rik", "ensembl_gene_id"] <- "ENSMUSG00000102869"
  DEGmarkers[DEGmarkers$gene == "Fam129a", "ensembl_gene_id"] <- "ENSMUSG00000026483"
  DEGmarkers[DEGmarkers$gene == "Fam49a", "ensembl_gene_id"] <- "ENSMUSG00000020589"
  DEGmarkers[DEGmarkers$gene == "Impad1", "ensembl_gene_id"] <- "ENSMUSG00000066324"
  DEGmarkers[DEGmarkers$gene == "Ndufb1-ps", "ensembl_gene_id"] <- "ENSMUSG00000113902"
  DEGmarkers[DEGmarkers$gene == "CT010467.1", "ensembl_gene_id"] <- "ENSMUSG00000045999"
  DEGmarkers[DEGmarkers$gene == "Plac9b", "ensembl_gene_id"] <- "ENSMUSG00000095304"
  
  
  DEGmarkers <- rename(DEGmarkers,c( 'ensembl_gene_id'=  'ENS'))
  
  return(DEGmarkers)
}


MapENSEMBLtoENTREZID <- function(DEgenes, organismDB){
  entrezids<-bitr(DEgenes$ENS, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organismDB)
  mapped_ids = entrezids[!duplicated(entrezids[c("ENSEMBL")]),]
  names(mapped_ids)[names(mapped_ids) == 'ENSEMBL'] <- 'ENS'
  DEgenes <-merge(DEgenes,mapped_ids, by = "ENS")
  
  return(DEgenes)
}

Enrichment_Analysis <- function(anal_type,DEmarkers, analysisType, analysisData, organismDB, enrichcl_dict,enrichcl_list, enrichcl_obj_list, disease_phase,output_path) {
  
  # Filter markers on anal_type
  DEmarkers <- subset(DEmarkers, cluster== anal_type)
  
  Ensembl_id_presence <- DEmarkers%>%mutate(gene=map(., ~str_detect(.x, 'ENSMUSG'))%>%pmap_lgl(any))
  
  if (any(sum(Ensembl_id_presence$gene))){
    ENS <- vapply(strsplit(DEmarkers$gene,"\\."), `[`, 1, FUN.VALUE=character(1))
    SYMBOL <- vapply(strsplit(DEmarkers$gene,"\\."), `[`, 2, FUN.VALUE=character(1))
    DEmarkers <- cbind(DEmarkers, ENS,SYMBOL ) %>% slice_min(., p_val_adj, n=250)
    
  }
  else{
    
    DEmarkers <- MapGenetoENSEMBL(DEmarkers,organismDB,output_path) %>% slice_min(., p_val_adj, n=250)
  }
  
  
  if (analysisType == "GO" & analysisData == 'ENSEMBL' ) {
    # Perform GO Enrichment Analysis
    enrichObj <- enrichGO(gene  = unique(DEmarkers$ENS),
                          OrgDb         = organismDB,
                          keyType       = analysisData,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.05, #pvalueCutoff  = 0.05
                          qvalueCutoff  = 0.05) #qvalueCutoff  = 0.05
  }
  else if (analysisType == "GO" & analysisData == 'SYMBOL' ) {
    # Perform GO Enrichment Analysis
    enrichObj <- enrichGO(gene  = unique(DEmarkers$SYMBOL),
                          OrgDb         = organismDB,
                          keyType       = analysisData,
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable = FALSE)
  }
  else if (analysisType == "KEGG" & analysisData == 'ncbi-geneid') {
    DEmarkers <- MapENSEMBLtoENTREZID(DEmarkers,organismDB)
    
    # Perform KEGG pathways Enrichwment Analysis
    enrichObj <- enrichKEGG(gene= unique(DEmarkers$ENTREZID),
                            organism= "mouse",
                            keyType       = analysisData,
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05)
  }
  
  # Add qscore, GeneRatio, analysis, color as S4 Class elements
  sc1 <- setClass('qscore', contains = 'enrichResult', slots = c(qscore = 'numeric'))
  enrichPlusqscore <- as(enrichObj, 'qscore')
  enrichPlusqscore@result$qscore <- unlist(lapply(enrichPlusqscore@result$p.adjust,function(x) -log(x,base=10)))
  
  sc2 <- setClass('analysis', contains = 'qscore', slots = c(grp = 'numeric'))
  enrichclPlusanal <- as(enrichPlusqscore, 'analysis')
  enrichclPlusanal@result$analysis <- c(replicate(length(enrichclPlusanal@result$Count), anal_type))
  
  sc3 <- setClass('GeneRatioNum', contains = 'analysis', slots = c(GeneRatioNum = 'numeric'))
  enrichclPlusgration <- as(enrichclPlusanal, 'GeneRatioNum')
  enrichclPlusgration@result$GeneRatioNum <- parse_ratio(c(enrichclPlusgration@result$GeneRatio))
  
  sc4 <- setClass('Color', contains = 'GeneRatioNum', slots = c(Color = 'numeric'))
  enrichcl <- as(enrichclPlusgration, 'Color')

  color <- "orange"
  
  enrichcl@result$Color <- c(replicate(length(enrichcl@result$Count), color))
  
  enrichcl_dict[["qscore_min"]] <- c(enrichcl_dict[["qscore_min"]],min(enrichcl@result$qscore))
  enrichcl_dict[["qscore_max"]] <- c(enrichcl_dict[["qscore_max"]],max(enrichcl@result$qscore))
  
  enrichcl_dict[["GeneRatioNum_min"]] <- c(enrichcl_dict[["GeneRatioNum_min"]],min(enrichcl@result$GeneRatioNum))
  enrichcl_dict[["GeneRatioNum_max"]] <- c(enrichcl_dict[["GeneRatioNum_max"]],max(enrichcl@result$GeneRatioNum))
  
  enrichcl_dict[["Count_min"]] <- c(enrichcl_dict[["Count_min"]],min(enrichcl@result$Count))
  enrichcl_dict[["Count_max"]] <- c(enrichcl_dict[["Count_max"]],max(enrichcl@result$Count))
  
  enrichcl_list <- c(enrichcl_list,enrichcl)
  
  enrichcl_obj_list <- c(enrichcl_obj_list,enrichObj)
  
  print("End Enrichment analysis ...")
  return(list("enrichcl_dict" = enrichcl_dict, "enrichcl_list" = enrichcl_list, "enrichcl_obj_list" = enrichcl_obj_list))
}


