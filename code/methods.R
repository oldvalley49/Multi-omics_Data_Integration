.libPaths("/users/tfurutan/R/4.3")
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(rliger)
library(bindSC)
library(irlba)
library(umap)
#install additional libraries
library(assertthat)

getSubsample <- function(seurat.obj, k) {
  
  # calculate number of cells in original dataset
  num.cells <- length(Cells(seurat.obj))
  # calculate hte number of cells for the 1/kth sample
  subsample.size <- num.cells %/% k
  pbmc.subsample <- pbmc[, sample(colnames(pbmc), size =subsample.size, replace=False)]
  print(paste0("Subsample size: ", length(Cells(pbmc.subsample))))
  return(pbmc.subsample)
  
}

getFragmentsENCODE <- function(atac) {
  orig.ident <- unique(atac@meta.data[["orig.ident"]])
  fragments <- list()
  file.path <- "/scMulti/ENCODE/cellranger_human/"
  extract_barcode <- function(x) {
    substr(x, 13, nchar(x))
  }
  for (origin in orig.ident) {
    file.path.origin <- paste0(file.path, origin, "/", "atac_fragments.tsv.gz")
    cells.orig <- colnames(subset(atac, subset = orig.ident == origin))
    cells.orig <- sapply(cells.orig, extract_barcode)
    print(origin)
    print(file.path.origin)
    print(length(cells.orig))
    fragments[[origin]] <- CreateFragmentObject(path = file.path.origin, cells = cells.orig)
  }
  Fragments(atac) <- fragments
  return(atac)
}

getFragmentsPBMC <- function(atac) {
  frag <- CreateFragmentObject(path = "/users/tfurutan/data/fragments.tsv.gz")
  Fragments(atac) <- frag
  return(atac)
}

QC <- function(seurat.obj){
  seurat.obj <- subset(
    x = seurat.obj,
    subset = nCount_ATAC < 7e4 &
      nCount_ATAC > 5e3 &
      nCount_RNA < 25000 &
      nCount_RNA > 1000 &
      percent.mt < 20
  )
  return(seurat.obj)
}

workflowSeurat <- function(seurat.obj, variable.num) {
  
  DefaultAssay(seurat.obj) <- "RNA"
  #RNA workflow
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj, nfeatures = variable.num)
  seurat.obj <- ScaleData(seurat.obj)
  seurat.obj <- RunPCA(seurat.obj)
  
  DefaultAssay(seurat.obj) <- "ATAC"
  #ATAC workflow
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  Annotation(seurat.obj) <- annotations
  seurat.obj <- RunTFIDF(seurat.obj)
  seurat.obj <- FindTopFeatures(seurat.obj, min.cutoff = "q0")
  seurat.obj <- RunSVD(seurat.obj)
  
  DefaultAssay(seurat.obj) <- "RNA"
  return(seurat.obj)
}

loadData <- function(tissue){
  if (tissue == "pbmc"){
    multiome.data <- readRDS("/users/tfurutan/data/pbmc_seurat.rds")
  }
  else {
    multiome.data <- readRDS("/dcs04/hongkai/data/zfu17/scMulti/ENCODE/proc_data/tissue_integration_list.rds")
    multiome.data <- multiome.data[[tissue]]
    multiome.data@assays[["RNA"]] <- as(object = multiome.data@assays[["RNA"]], Class = "Assay5")
  }
  return(multiome.data)
}

runSeurat <- function(multiome.data, tissue, var.num) {
  #QC
  multiome.data <- QC(multiome.data)
  
  #split
  DefaultAssay(multiome.data) <- "RNA"
  rna <- DietSeurat(multiome.data, assays = "RNA")
  DefaultAssay(multiome.data) <- "ATAC"
  atac <- DietSeurat(multiome.data, assays = "ATAC")
  
  rna[["orig.assay"]] <- "RNA"
  atac[["orig.assay"]] <- "ATAC"
  
  #Get fragments
  if (tissue == "pbmc"){
    atac <- getFragmentsPBMC(atac)
  }
  else {
    atac <- getFragmentsENCODE(atac)
  }
  #perform standard workflow
  rna <- NormalizeData(rna)
  rna <- FindVariableFeatures(rna, nfeatures = var.num)
  rna <- ScaleData(rna)
  rna <- RunPCA(rna)
  
  #ATAC workflow
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  Annotation(atac) <- annotations
  atac <- RunTFIDF(atac)
  atac <- FindTopFeatures(atac, min.cutoff = "q0")
  atac <- RunSVD(atac)
  
  use.features = VariableFeatures(rna)
  #impute expression from ATAC
  gene.activities <- GeneActivity(atac, features = use.features)
  # add gene activities as a new assay
  atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
  # normalize gene activities
  DefaultAssay(atac) <- "ACTIVITY"
  atac <- NormalizeData(atac)
  atac <- ScaleData(atac, features = rownames(atac))
  rna <- DietSeurat(rna, features = use.features)
  
  transfer.anchors <- FindTransferAnchors(reference = rna, query = atac, features = use.features, reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
  
  refdata <- GetAssayData(rna, assay = "RNA", slot = "data")[use.features, ]
  
  imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata, weight.reduction = atac[["lsi"]],
                             dims = 2:30)
  atac[["RNA"]] <- imputation
  
  coembed <- merge(x = rna, y = atac, add.cell.ids = c("RNA", "ATAC"))
  #remove unneeded objects
  rm(atac)
  rm(rna)
  coembed <- ScaleData(coembed, features = use.features, do.scale = FALSE)
  coembed <- RunPCA(coembed, features = use.features, verbose = FALSE)
  #output cell embedding as csv file
  write.csv(Embeddings(coembed, reduction = "pca"), file = paste0('/users/tfurutan/integration/var/coembed/', tissue, "/", 'seurat_v3_v', var.num, '.csv'), row.names = TRUE)
  return(length(use.features))
  #make sure to return cor.num at some point
  
}

runLiger <- function(multiome.data, tissue, var.num){
  multiome.data <- QC(multiome.data)

  DefaultAssay(multiome.data) <- "RNA"
  rna <- DietSeurat(multiome.data, assays = "RNA")
  DefaultAssay(multiome.data) <- "ATAC"
  atac <- DietSeurat(multiome.data, assays = "ATAC")

  rna[["orig.assay"]] <- "RNA"
  atac[["orig.assay"]] <- "ATAC"

  if (tissue == "pbmc"){
    atac <- getFragmentsPBMC(atac)
  }
  else {
    atac <- getFragmentsENCODE(atac)
  }
  
  rna <- NormalizeData(rna)
  rna <- FindVariableFeatures(rna, nfeatures = var.num)
  use.features <- VariableFeatures(rna)
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  Annotation(atac) <- annotations
  
  gene.activities <- GeneActivity(atac, features = use.features)
  atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
  #set unique names
  colnames(rna) = paste0("RNA_", colnames(rna))
  colnames(atac) = paste0("ATAC_", colnames(atac))
  
  rna_counts <- rna[["RNA"]]$counts
  
  liger.obj <- createLiger(list(atac = atac@assays[["ACTIVITY"]]@counts, rna = rna))
  rm(atac)
  #integrate
  liger.obj <- normalize(liger.obj)
  liger.obj <- selectGenes(liger.obj, useDatasets = "rna", thresh=0)
  liger.obj <- scaleNotCenter(liger.obj)
  
  liger.obj <- runIntegration(liger.obj, k = 20)
  liger.obj <- quantileNorm(liger.obj)
  liger.obj <- runCluster(liger.obj, nNeighbors = 30, resolution = 0.2)
  
  #liger.obj <- runUMAP(liger.obj, nNeighbors = 30, minDist = 0.3)
  #dir.create(paste0("/users/tfurutan/coembed/", tissue))
  #dir.create(paste0("/users/tfurutan/rna_counts/", tissue))
  
  write.csv(liger.obj@H.norm, file = paste0('/users/tfurutan/integration/var/coembed/', tissue, "/", 'liger_v', var.num, '.csv'), row.names = TRUE)
}

runbindSC <- function(multiome.data, tissue, var.num){
  multiome.data <- QC(multiome.data)

  DefaultAssay(multiome.data) <- "RNA"
  rna <- DietSeurat(multiome.data, assays = "RNA")
  DefaultAssay(multiome.data) <- "ATAC"
  atac <- DietSeurat(multiome.data, assays = "ATAC")

  rna[["orig.assay"]] <- "RNA"
  atac[["orig.assay"]] <- "ATAC"

  if (tissue == "pbmc"){
    atac <- getFragmentsPBMC(atac)
  }
  else {
    atac <- getFragmentsENCODE(atac)
  }
  
  rna <- NormalizeData(rna)
  rna <- FindVariableFeatures(rna, nfeatures = var.num)
  rna <- ScaleData(rna)
  rna <- RunPCA(rna)
  rna <- FindNeighbors(rna, dims = 1:20, reduction="pca")
  rna <- FindClusters(rna, resolution = 0.5)
  
  use.features <- VariableFeatures(rna)
  
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  Annotation(atac) <- annotations
  
  atac <- RunTFIDF(atac)
  atac <- FindTopFeatures(atac, min.cutoff = "q0")
  atac <- RunSVD(atac)
  atac <- FindNeighbors(atac, dims = 1:20, reduction = "lsi")
  atac <- FindClusters(atac, resolution = 0.5)
  
  gene.activities <- GeneActivity(atac, features = use.features)
  atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
  
  #set unique names
  colnames(rna) = paste0("RNA_", colnames(rna))
  colnames(atac) = paste0("ATAC_", colnames(atac))
  
  DefaultAssay(atac) <- "ACTIVITY"
  atac <- NormalizeData(atac)
  atac <- ScaleData(atac, features = rownames(atac))
  
  genes.use <- intersect(use.features, rownames(atac))
  X <- rna[["RNA"]]$counts[genes.use, ]
  Y <- atac[["ATAC"]]$counts
  Z0 <- atac[["ACTIVITY"]]$counts[genes.use, ]
  
  
  K <- 30
  out <- irlba(Matrix::crossprod(X,Z0), nu =K, nv =K) #check if this sparse, return sparse?
  rownames(out$u) <- colnames(X)
  rownames(out$v) <- colnames(Z0)
  colnames(out$u) <- paste("C",seq(1,K,1), sep="")
  colnames(out$v) <- paste("C",seq(1,K,1), sep="")
  result <- list()
  result$X <- out$u
  result$Z0 <- out$v
  
  x <- result$X
  z0 <- result$Z0
  y  <- atac@reductions[["lsi"]]@cell.embeddings[, 1:30]
  
  
  res <- BiCCA( X = t(x) ,
                Y = t(y), 
                Z0 =t(z0), 
                X.clst = rna$seurat_clusters,
                Y.clst = atac$seurat_clusters,
                alpha = 0.5, 
                lambda = 0.5,
                K = 15,
                temp.path  = "out",
                num.iteration = 50,
                tolerance = 0.01,
                save = TRUE,
                parameter.optimize = FALSE,
                block.size = 0)
  
  write.csv(rbind(res$u, res$r), file = paste0('/users/tfurutan/integration/var/coembed/', tissue, "/", 'bindsc_v', var.num, '.csv'), row.names = TRUE)
}
