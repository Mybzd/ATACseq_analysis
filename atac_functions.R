library(Seurat)
library(ggplot2)
library(patchwork)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(future)

#==============================================================================##
# CREATE MULTIOME OBJECTS
#==============================================================================##

# 1.seuOB = list of h5 matrices paths
# 2.frag.path = list of fragment paths

CreateSeuraMultiome=function(h5_mat, metadata, frag.path, min.cells.rna, min.cells.atac, annotation, outdir, save, name){

message("Reading in h5 data...")
  counts <- lapply(h5_mat, Read10X_h5)

  # add in RNA
message("Creating Seurat Object with RNA...")
  seuOBs <- lapply(counts, function(x){CreateSeuratObject(x$'Gene Expression', min.cells=min.cells.rna, metadata=metadata)})
  names(seuOBs)=lapply(h5_mat, function(x){strsplit(x,"/")[[1]][2]})

  # add in ATAC
message("Adding in ATAC...")
  for (i in 1:length(seuOBs)){
   grangeCounts <- list()
   grangeUse <- list()
   atacCounts <- list()
   grangeCounts[[i]] <- StringToGRanges(rownames(counts[[i]]$Peaks), sep = c(":", "-"))
   grangeUse[[i]] <- seqnames(grangeCounts[[i]]) %in% standardChromosomes(grangeCounts[[i]])
   atacCounts[[i]] <- counts[[i]]$Peaks[as.vector(grangeUse[[i]]),]  

   seuOBs[[i]][["ATAC"]]= CreateChromatinAssay(
     counts = atacCounts[[i]],
     sep = c(":", "-"),
     fragments = frag.path[[i]],
     min.cells = min.cells.atac,
     annotation = annotation)
     }

  if (save==TRUE) {saveRDS(seuOBs, paste0(outdir,name,"_multiome_object.rds"))}
  return(seuOBs)
}


#==============================================================================##
# CALL PEAKS
#==============================================================================##

# 1.seuOB = list of seuOBS
# 2.frag.path = list of fragment paths

CallPeaksMACs=function(seuOBs, frag.path, macs2.path, annotation, outdir, save, name){

message("Reading in data...")

for (i in 1:length(seuOBs)){
  DefaultAssay(seuOBs[[i]])<-"ATAC"
  }

message("Calling peaks...")
  peaks <- lapply(seuOBs, CallPeaks, macs2.path = macs2.path)
  # around 15min of runtime

  # remove peaks on non standard chromosomes and regions
message("MACS done...")
  peaks <- lapply(peaks, keepStandardChromosomes, pruning.mode = "coarse")
message("pt.2...")
  peaks <- lapply(peaks, subsetByOverlaps, ranges=blacklist_hg38_unified, invert = TRUE)

message("Saving...")
if (save==TRUE) {saveRDS(peaks, paste0(outdir,name,"_multiome_Peaks.rds"))}
}


# reduce peaks
#message("Combining peaks...")

# combinedPeaks=reduce(x=unlist(peaks[[1]],peaks[[2]]) etc..) - (not working using Reduce)

#==============================================================================##
# QUANTIFY PEAKS
#==============================================================================##

QuantifyPeaks=function(seuOBs, frag.path, combinedPeaks, annotation, outdir, save, name){

message("Quantifying counts...")
  macsCounts=list()
   for (i in 1:length(seuOBs)){
    DefaultAssay(seuOBs[[i]])<-"ATAC"
    macsCounts[[i]] <- FeatureMatrix(
      fragments = Fragments(seuOBs[[i]]),
      features = combinedPeaks,
      cells = colnames(seuOBs[[i]])
       )

   seuOBs[[i]][["peaks"]] <- CreateChromatinAssay(
    counts = macsCounts[[i]],
    fragments = frag.path[[i]],
    annotation = annotation
    )
   } 
  if (save==TRUE) {saveRDS(seuOBs, paste0(outdir,name,"_multiome_object_PeaksMX.rds"))}
  return(seuOBs)
}

#==============================================================================##
# NORMALISE
#==============================================================================##

NormaliseATAC=function(seuOBs, save, name){

   applyNorm=function(x){
    DefaultAssay(x) <- "peaks"
    x <- FindTopFeatures(x, min.cutoff = "q0")
    x <- RunTFIDF(x)
    x <- RunSVD(x)
     }

  seuObs=lapply(seuObs, applyNorm)
  
  if (save==TRUE) {saveRDS(seuOBs, paste0((outdir,name,"multiome_SEUobject_normalised.rds"))}
  return(seuOBs)  
}


#==============================================================================##
# SUBSET OBS
#==============================================================================##

SubsetCells=function(seuOBs, meta, save, name){

seuOBsubed=setNames(vector(mode="list", length=length(seuOBs)),names(seuOBs))

pools=lapply(names(seuObs), function(x){strsplit(x,"_")[[1]][1]})

    for (i in 1:length(seuObs)){
      pool=pools[[i]] 

s1=levels(as.factor(meta[which(meta$Pool==pool),]$Sample))[1]
s2=levels(as.factor(meta[which(meta$Pool==pool),]$Sample))[2]

n1=levels(as.factor(meta[which(meta$Sample==n1),]$sample_name))[1]
n2=levels(as.factor(meta[which(meta$Sample==n2),]$sample_name))[2]

cellsUse1=meta$Cell[which(meta$Sample==n1)]
cellsUse2=meta$Cell[which(meta$Sample==n2)]

SEUob1=subset(seuOBs[[i]], cells=cellsUse1)
SEUob2=subset(seuOBs[[i]], cells=cellsUse2)

	renamedCOLS1=vector()
	renamedCOLS2=vector()

   for (k in 1:length(colnames(SEUob1@assays$RNA@counts))){
     renamedCOLS1[k]=paste0(s1,"_",pool,"_",colnames(SEUob1@assays$RNA@counts)[k])}

   for (k in 1:length(colnames(SEUob2@assays$RNA@counts))){
     renamedCOLS2[k]=paste0(s2,"_",pool,"_",colnames(SEUob2@assays$RNA@counts)[k])}

   SEUob1=RenameCells(SEUob1, new.names=renamedCOLS1)
   SEUob2=RenameCells(SEUob2, new.names=renamedCOLS2)
   
  seuOBsubed[[i]]=list(SEUob1,SEUob2)
}

if (save==TRUE) {saveRDS(seuOBsubed, paste0(outdir,name,"_multiome_object_subsetted.rds"))}
return(seuOBsubed)
    }


#==============================================================================##
# MERGE
#==============================================================================##

# also not sure how to do this in a function
#combined <- merge(
#  x = seuObs[[1]],
#  y = list(seuObs[[2]],seuObs[[3]],seuObs[[4]],seuObs[[5]],seuObs[[6]],seuObs[[7]],seuObs[[8]],seuObs[[9]],seuObs[[10]],
#          seuObs[[11]],seuObs[[12]],seuObs[[13]],seuObs[[14]],seuObs[[15]],seuObs[[16]],seuObs[[17]],seuObs[[18]])
#)

#==============================================================================##
# NUCLEOSOME FILTERING
#==============================================================================##

NucleosomeTSSfiltering=function(combinedOB, nN,nT, save, name){

DefaultAssay(combinedOB) <- "ATAC"

combinedOB <- NucleosomeSignal(object = combinedOB)

# compute TSS enrichment score per cell
combinedOB <- TSSEnrichment(object =combinedOB, fast=FALSE)


combinedOB <- subset(
  x = combinedOB,
  subset =nucleosome_signal < nN & TSS.enrichment > nT
)

if (save==TRUE) {saveRDS(combinedOB, paste0(outdir,name,"_multiome_combinedObject_TTSnucleosome.rds"))}

return(combinedOB)
    }


#==============================================================================##
# NORMALISE COMBINED SEURAT
#==============================================================================##

NormaliseCombined=function(combinedOB, save, name){

DefaultAssay(combinedOB) <- "peaks"

combinedOB <- FindTopFeatures(combinedOB, min.cutoff = "q0")
combinedOB <- RunTFIDF(combinedOB)
combinedOB <- RunSVD(
  combinedOB,
  reduction.key = 'LSI_',
  reduction.name = 'lsi',
  irlba.work = 400
)
combinedOB <- RunUMAP(combinedOB, dims = 2:30, reduction = 'lsi')

if (save==TRUE) {saveRDS(combinedOB, paste0(outdir,name,"_multiome_combinedObject_normalised.rds"))}

return(combinedOB)

    }

#==============================================================================##
# INTEGRATE SAMPLES
#==============================================================================##

IntegrateObs=function(combinedOB,seuObs, save, name){

integration.anchors <- FindIntegrationAnchors(
  object.list = seuObs,
  anchor.features = rownames(seuObs[[1]]),
  reduction = "rlsi",
  dims = 2:30
)

if (save==TRUE) {saveRDS(combinedOB, paste0(outdir,name,"_multiome_combinedObject_integration.anchors.rds"))}

DefaultAssay(combinedOB) <- "peaks"
# integrate LSI embeddings
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = combinedOB[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

if (save==TRUE) {saveRDS(integrated, paste0(outdir,name,"_multiome_combinedObject_integrated.rds"))}

return(integrated)

    }

#==============================================================================##
# NORMALISE INTEGRATED
#==============================================================================##

normaliseIntegrated=function(integrated, save, name){

integrated=NormalizeData(integrated, normalization.method = "LogNormalize",verbose = FALSE)
integrated=RunPCA(integrated)
integrated=FindNeighbors(integrated,dims=1:50)
integrated=FindClusters(integrated,pc.use=1:50)
integrated=RunUMAP(integrated, pc.use=1:50, reduction = "umap")

if (save==TRUE) {saveRDS(integrated, paste0(outdir,name,"_multiome_combinedObject_integrated_normalised.rds"))}

return(integrated)

    }

#==============================================================================##
# NORMALISE INTEGRATED
#==============================================================================##
normaliseMultiome=function(integrated, save, name){

DefaultAssay(integrated)<-"ATAC"
integrated <- FindTopFeatures(integrated, min.cutoff = "q0")
integrated <- RunTFIDF(integrated)
integrated <- RunSVD(integrated)
integrated <- RunUMAP(integrated, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

DefaultAssay(integrated) <- "RNA"
integrated <- SCTransform(integrated)
integrated <- RunPCA(integrated)

if (save==TRUE) {saveRDS(integrated, paste0(outdir,name,"_multiome_combinedObject_integrated_normalisedMULTIOMEob.rds"))}

return(integrated)

    }

#==============================================================================##
# TRANSFER ANCHORS 
#==============================================================================##

transferAnchors=function(integrated,reference, save, name){

DefaultAssay(integrated) <- "SCT"

# transfer cell type labels from reference to query
transfer_anchors <- FindTransferAnchors(
  reference = reference,
  query = integrated,
  normalization.method = "SCT",
  reference.reduction = "pca",
  recompute.residuals = FALSE,
  dims = 1:50
)

predictions <- TransferData(
  anchorset = transfer_anchors, 
  refdata = reference$celltype,
  weight.reduction = integrated[['pca']],
  dims = 1:50
)

integrated <- AddMetaData(
  object = integrated,
  metadata = predictions
)

if (save==TRUE) {
saveRDS(transfer_anchors, paste0(outdir,name,"_multiome_combinedObject_integrated_transferAnchors.rds"))
saveRDS(predictions, paste0(outdir,name,"_multiome_combinedObject_integrated_predictions.rds")) 
saveRDS(integrated, paste0(outdir,name,"_multiome_combinedObject_integrated_labelsAdded.rds")) 
}

    }

#==============================================================================##
# MULTIOMODAL NEIGHBORS
#==============================================================================##

multimodalNeighbors=function(integrated, save, name){

# build a joint neighbor graph using both assays
integrated <- FindMultiModalNeighbors(
  object = integrated,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# build a joint UMAP visualization
integrated <- RunUMAP(
  object = integrated,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

if (save==TRUE) {saveRDS(integrated, paste0(outdir,name,"_multiome_combinedObject_integrated_multimodalUMAP.rds"))}

    }

