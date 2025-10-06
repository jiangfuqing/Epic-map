#take Dox-rep3 as example
#load package
pkgs <- c("Signac","Seurat","ArchR","Rsamtools","GenomeInfoDb","EnsDb.Hsapiens.v75","EnsDb.Mmusculus.v79","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Mmusculus.UCSC.mm10","future","ggplot2","dplyr","reshape2","stringr")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

# accept args from pbs file
args <- commandArgs(TRUE)

sampleId <- args[1]

plan("multicore", workers = 46)
options(future.globals.maxSize = 50 * 1024 ^ 3) # for 5 Gb RAM

set.seed(1234)

# load the RNA and ATAC data
counts <- Read10X_h5(paste0(sampleId,"-raw_peak_bc_matrix.h5"))
fragpath <- paste0(sampleId,"-filtered_fragments.tsv.gz")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
    
# create ATAC assay and add it to the object
atac <- CreateChromatinAssay(counts = counts, sep = c(":", "-"), fragments = fragpath, annotation = annotation)
# Create seurat object
df <- CreateSeuratObject(counts = atac, assay = "ATAC")

# filter non-tissue grids use position_{sampleId}.txt file which generated from photoshop and matlab soft.
location <- read.table(paste0("position_",sampleId,".txt"), sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])[-1]
x <- as.data.frame(x)
colnames(x) <- "array"
y <- read.csv("~/barcode/Epic_barcode_filter.csv",header = T) 
z <- merge(x,y,by="array")

df <- subset(df, cells=z$barcode)

# The set of peaks identified using Cellranger often merges distinct peaks that are close together. This can create a problem for certain analyses, particularly motif enrichment analysis and peak-to-gene linkage.
# call peaks using MACS2
peaks <- CallPeaks(df, macs2.path = "~/anaconda3/envs/python39/bin/macs2") # effective.genome.size=1.87e+09,

# remove peaks on nonstandard chromosomes and in genomic blacklist regions, mouse:blacklist_mm10, human:blacklist_hg38_unified.
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(df),
  features = peaks,
  cells = colnames(df)
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
df[["peaks"]] <- CreateChromatinAssay(counts = macs2_counts,fragments = fragpath,annotation = annotation)

# read spatial barcode file.
barcode <- read.csv("~/barcode/Epic_barcode.csv",header = T,row.names = 1)

# read and filter per_barcode_metrics.csv file which from cellranger
metrics <- read.csv(paste0(sampleId,"-singlecell.csv"),header = T,row.names = 1)

metrics <- merge(metrics,barcode,by="row.names")
# used when your input is filter_fragments.tsv
# metrics <- metrics[metrics$Row.names %in% location,]

rownames(metrics) <- metrics$Row.names

metrics <- transform(metrics, atac_percent_mito=(metrics$mitochondrial/metrics$total)*100, FRiP=(metrics$peak_region_fragments/metrics$passed_filters)*100, percent_TSS_fragments=(metrics$TSS_fragments/metrics$passed_filters)*100, log10_unique_fragments=log10(metrics$passed_filters), Sample=sampleId)

metrics <- metrics[,c("FRiP", "log10_unique_fragments", "atac_percent_mito", "percent_TSS_fragments", "Sample", "array_col", "array_row" )]

df <- AddMetaData(object = df, metadata = metrics)

#perform ATAC reduction and cluster analysis
#The set of peaks identified using Cellranger often merges distinct peaks that are close together. This can create a problem for certain analyses, particularly motif enrichment analysis and peak-to-gene linkage. So we use macs2 to identify a more accurate set of peaks.

DefaultAssay(df) <- "peaks"
for (i in c("q5", "q15")){
df <- FindTopFeatures(df, min.cutoff = i) %>% RunTFIDF() %>% RunSVD()
df <- FindNeighbors(df, reduction = 'lsi', dims = 2:30) %>% FindClusters(resolution = c(0.5,0.8,1), verbose = FALSE, algorithm = 3) %>% RunUMAP(reduction = 'lsi', dims = 2:30, reduction.name = paste0("umap.",i))
df@meta.data[paste0(i,"_snn_res.0.5")] <- df@meta.data["peaks_snn_res.0.5"]
df@meta.data[paste0(i,"_snn_res.0.8")] <- df@meta.data["peaks_snn_res.0.8"]
df@meta.data[paste0(i,"_snn_res.1")] <- df@meta.data["peaks_snn_res.1"]
}
df@meta.data["peaks_snn_res.0.5"] <- df@meta.data["peaks_snn_res.0.8"] <- df@meta.data["peaks_snn_res.1"] <- NULL 

# compute gene accessibility using fragments counts and add to the Seurat
gene.activities <- GeneActivity(df)
df[['gene_activities']] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(df) <- "gene_activities"
df <- NormalizeData(df, assay = 'gene_activities', normalization.method = 'LogNormalize', scale.factor = median(df$nCount_gene_activities))
df <- ScaleData(df, features = rownames(df))
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000) %>% RunPCA(verbose = FALSE) %>% FindNeighbors(dims = 1:15) %>% FindClusters(resolution = c(0.5,0.8,1.0)) %>% RunUMAP(dims = 1:15, reduction.name = "umap.gene_activities")

# plot Cluster in each resolution
#imported_raster=OpenImageR::readImage(paste0(sampleId,"-HE.jpg"))
#g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

#merge clusters which have no meaningful for Dox-rep3
df$merge_cluster <- df$q5_snn_res.1
df$merge_cluster[df$merge_cluster %in% c("0","2")] <- "0"
df$merge_cluster[df$merge_cluster %in% "3"] <- "2"
df$merge_cluster[df$merge_cluster %in% "4"] <- "3"
df$merge_cluster <- as.character(as.numeric(as.character(df$merge_cluster))+1)

clusters <- "merge_cluster"

pdf(file = paste0(sampleId,"_Cluster_each.pdf"), width=4.6, height=4.0)
for (i in clusters) {
  p1 <- DimPlot(df, reduction = "umap.q5", group.by =i, label = T, pt.size = 1.5, label.size=4) + ggtitle(i) + theme(plot.title = element_text(hjust = 0.5,size = 15))
  plot(p1)
}
dev.off()

pdf(file = paste0(sampleId,"_Signac_Spatial_Clusters.pdf"), width=8.6, height=8.6)
for (i in clusters){
  p2 <- ggplot(df@meta.data, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df@meta.data[,i])) +
    #scale_color_gradientn(colours = c("black", "green")) + 
    #scale_color_gradientn(colours = c("blue","green", "red"),
    #                      oob = scales::squish) +
    ggtitle(i) +
    #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(size = 3.5,shape=19)+
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name=NULL, limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name=NULL,limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          #legend.title = element_text(colour="black", size=15, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  plot(p2)
}
dev.off()

#save raw seurat object
saveRDS(df, file = paste0(sampleId,"_Signac.rds"))