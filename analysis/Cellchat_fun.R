
library(tibble)
library(tidyr)
library(circlize)
library(magrittr)

LoadCellChatDB <- function(cellchat,type){
  if(type == "human"){
    CellChatDB <- CellChatDB.human
  }else{
    CellChatDB <- CellChatDB.mouse
  }
  showDatabaseCategory(CellChatDB)
  # use all CellChatDB for cell-cell communication analysis
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  return(cellchat)
}

DataPreprocessing <- function(cellchat){
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network
  cellchat <- projectData(cellchat, PPI.human)
  return(cellchat)
}

CellToCellCommunicationInference <- function(cellchat){
  # Compute the communication probability and infer cellular communication network
  cellchat <- computeCommunProb(cellchat,population.size = FALSE)
  # Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat )
  # Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  
  # Compute network centrality scores - slot 'netP' : inferred intercellular communication network of signaling pathways
  cellchat <- netAnalysis_computeCentrality(cellchat)
  
  return(cellchat)
}

CellToCellCommunicationAnalysis <- function(cellchat,type){
  # Preprocessing of Expression data
  cellchat <- LoadCellChatDB(cellchat,type)
  print("CellChatDB Loaded!")
  
  # Set the ligand-receptor interaction database
  cellchat <- DataPreprocessing(cellchat)
  print("CellChatDB preprocessed!")
  
  # Inference of cell-cell communication network
  cellchat <- CellToCellCommunicationInference(cellchat)
  
  return(cellchat)
}

computeCommunProbPathway<- function (object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05)
{
  if (is.null(net)) {
    net <- object@net
  }
  if (is.null(pairLR.use)) {
    pairLR.use <- object@LR$LRsig
  }
  prob <- net$prob
  prob[net$pval > thresh] <- 0
  pathways <- unique(pairLR.use$pathway_name)
  group <- factor(pairLR.use$pathway_name, levels = pathways)
  prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum),
                         c(2, 3, 1))
  pathways.sig <- pathways[apply(prob.pathways, 3, sum) !=
                             0]
  prob.pathways.sig <- prob.pathways[, , pathways.sig]
  idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE,
              index.return = TRUE)$ix
  pathways.sig <- pathways.sig[idx]
  prob.pathways.sig <- prob.pathways.sig[, , idx]
  if (is.null(object)) {
    netP = list(pathways = pathways.sig, prob = prob.pathways.sig)
    return(netP)
  }
  else {
    object@netP$pathways <- pathways.sig
    object@netP$prob <- prob.pathways.sig
    return(object)
  }
}

MergeCellchat <- function(name1,name2,control.no.cellchat,TCRM.no.cellchat){

  object.list <- list(data1 = control.no.cellchat, data2 = TCRM.no.cellchat)
  names(object.list) <- c(name1, name2)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  return(list(object.list = object.list, cellchat = cellchat ))
}
