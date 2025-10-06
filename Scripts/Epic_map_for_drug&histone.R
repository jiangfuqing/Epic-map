#take Dox_H3K27ac-rep2 as example
#load package
pkgs <- c("Signac","Seurat","ArchR","Rsamtools","GenomeInfoDb","EnsDb.Hsapiens.v75","EnsDb.Mmusculus.v79","BSgenome.Hsapiens.UCSC.hg38","COSG",
          "BSgenome.Mmusculus.UCSC.mm10","future","ggplot2","dplyr","reshape2","stringr","clusterProfiler","org.Hs.eg.db","org.Mm.eg.db")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

# accept args from pbs file
args <- commandArgs(TRUE)

colnames1 <- args[1]
colnames2 <- args[2]
modality1 <- args[3]
modality2 <- args[4]

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))

color_list<-c("1"="#fcb462","2"="#9500ff","3"="#fb4f2b","4"="#4993b6","5"="#1e25ce","6"="#197a35","7"="#fd8000","8"="#d73027","9"="#8e0054","10"="#bb207c")

type <- "peaks" #"gene_activities" or "peaks"
colnameslist <- c(colnames1, colnames2)
modalitylist <- c(modality1, modality2)

drug <- paste0(colnameslist[1],"/fragments.tsv.gz")
histone <- paste0(colnameslist[2],"/fragments.tsv.gz")
fragpath <- list(drug,histone)
names(fragpath) <- modalitylist

counts1 <- Read10X_h5(paste0(colnameslist[1],"/raw_peak_bc_matrix.h5")) 
counts2 <- Read10X_h5(paste0(colnameslist[2],"/raw_peak_bc_matrix.h5")) 

atac1 <- CreateChromatinAssay(counts = counts1, sep = c(":", "-"), fragments = fragpath[[1]], annotation = annotation)
atac2 <- CreateChromatinAssay(counts = counts2, sep = c(":", "-"), fragments = fragpath[[2]], annotation = annotation)

df <- CreateSeuratObject(counts = atac1, assay = paste0(modalitylist[1],"_ATAC"))
df[[paste0(modalitylist[2],"_ATAC")]] <- atac2

# filter non-tissue grids use position_{sampleId}.txt file which generated from photoshop and matlab soft.
location <- read.table(paste0("position_",colnameslist[1],".txt"), sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])[-1]
x <- as.data.frame(x)
colnames(x) <- "array"
y <- read.csv("~/barcode/Epic_barcode_filter.csv",header = T) 
z <- merge(x,y,by="array")

df <- subset(df, cells=z$barcode)

# read spatial barcode file.
barcode <- read.csv("~/barcode/Epic_barcode.csv",header = T,row.names = 1)
df <- AddMetaData(df, metadata = barcode)

# call peaks using MACS2
for (i in modalitylist){
peaks <- CallPeaks(df[[paste0(i,"_ATAC")]], macs2.path = "~/anaconda3/envs/python37/bin/macs2") # effective.genome.size=1.87e+09,
# remove peaks on nonstandard chromosomes and in genomic blacklist regions, mouse:blacklist_mm10, human:blacklist_hg38_unified.
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(df[[paste0(i,"_ATAC")]]),
  features = peaks,
  cells = colnames(df[[paste0(i,"_ATAC")]])
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
df[[paste0(i,"_peaks")]] <- CreateChromatinAssay(counts = macs2_counts, fragments = fragpath[[i]], annotation = annotation)

DefaultAssay(df) <- paste0(i,'_peaks')
df <- FindTopFeatures(df) %>% RunTFIDF() %>%  RunSVD(reduction.name = paste0(i,'_lsi'))
df <- FindNeighbors(df, reduction = paste0(i,'_lsi'), dims = 2:30) %>% FindClusters(resolution = c(0.5,0.8,1.0,1.2), verbose = FALSE, algorithm = 3) %>% 
    RunUMAP(reduction = paste0(i,'_lsi'), dims = 2:30, reduction.name = paste0("umap.",paste0(i,'_peaks')))

# compute gene scores using fragments counts and add to Seurat.
df[[paste0(i,'_gene_activities')]] <- CreateAssayObject(counts = GeneActivity(df, assay=paste0(i,'_peaks')))
DefaultAssay(df) <- paste0(i,'_gene_activities')
df <- NormalizeData(df, normalization.method = 'LogNormalize', scale.factor = median(df@meta.data[,paste0('nCount_',i,'_gene_activities')]))
df <- ScaleData(df, features = rownames(df))
df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000) %>% RunPCA(verbose = FALSE, reduction.name = paste0(i,'_pca'))
df <- FindNeighbors(df, reduction = paste0(i,'_pca'), dims = 1:15) %>% FindClusters(resolution = c(0.5,0.8,1.0,1.2)) %>% 
    RunUMAP(reduction = paste0(i,'_pca'), dims = 1:15, reduction.name = paste0("umap.",paste0(i,'_gene_activities')))
}

#combined  analysis
if (type=="peaks"){
df <- FindMultiModalNeighbors(
  df, reduction.list = list(paste0(modalitylist[1],'_lsi'), paste0(modalitylist[2],'_lsi')), 
  dims.list = list(2:30, 2:30), modality.weight.name = "histone.weight", k.nn = 10
)} else if (type=="gene_activities"){
df <- FindMultiModalNeighbors(
  df, reduction.list = list(paste0(modalitylist[1],'_pca'), paste0(modalitylist[2],'_pca')), #l2.norm = FALSE
  dims.list = list(1:15, 1:15), modality.weight.name = "histone.weight", k.nn = 10
)
}

df <- FindClusters(df, graph.name = "wsnn", algorithm = 3, resolution = c(0.5,0.8,1.0,1.2), verbose = FALSE)
df <- RunUMAP(df, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key = "wnnUMAP_")

#merge some clusters which have no meaningful
#Drug
df$merge_cluster <- df$DOX_peaks_snn_res.1.2
df$merge_cluster[df$merge_cluster %in% c("0","3","4")] <- "0"
df$merge_cluster[df$merge_cluster %in% c("5")] <- "3"
df$merge_cluster[df$merge_cluster %in% c("6")] <- "4"
df$merge_cluster_DOX <- as.character(as.numeric(as.character(df$merge_cluster))+1)

#Histone
df$merge_cluster <- df$H3K27ac_peaks_snn_res.0.5
df$merge_cluster_H3K27ac <- as.character(as.numeric(as.character(df$merge_cluster))+1)

#Combined
df$merge_cluster <- df$wsnn_res.0.8
df$merge_cluster_Combined <- as.character(as.numeric(as.character(df$merge_cluster))+1)

df$merge_cluster <- NULL

#plot Cluster in each resolution
#imported_raster=OpenImageR::readImage(paste0(sampleId,"-HE.jpg"))
#g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

pdf(file = paste0(colnameslist[1],"&",colnameslist[2],"_UMAP_cluster.pdf"), width=8.6, height=8.6)
DimPlot(df, reduction = "umap.DOX_peaks", group.by = "merge_cluster_DOX", label = T, pt.size = 1.5, label.size=4) + 
    ggtitle("merge_cluster_DOX") + theme(plot.title = element_text(hjust = 0.5,size = 15)) + scale_color_manual(values = color_list)
DimPlot(df, reduction = "umap.H3K27me3_peaks", group.by = "merge_cluster_H3K27ac", label = T, pt.size = 1.5, label.size=4) + 
    ggtitle("merge_cluster_H3K27ac") + theme(plot.title = element_text(hjust = 0.5,size = 15)) + scale_color_manual(values = color_list)
DimPlot(df, reduction = "umap.wnn", group.by = "merge_cluster_Combined", label = T, pt.size = 1.5, label.size=4) + 
    ggtitle("merge_cluster_Combined") + theme(plot.title = element_text(hjust = 0.5,size = 15)) + scale_color_manual(values = color_list)
dev.off()

indexs <- grep("merge", colnames(df@meta.data), value = TRUE)
pdf(file = paste0(colnameslist[1],"&",colnameslist[2],"_Spatial_cluster.pdf"), width=8.6, height=8.6)
for (i in indexs){
  p5 <- ggplot(df@meta.data, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df@meta.data[,i])) +
    #scale_color_gradientn(colours = c("black", "green")) + 
    #scale_color_gradientn(colours = c("blue","green", "red"),
    #                      oob = scales::squish) +
    scale_color_manual(values = color_list) + 
    ggtitle(i) +
    #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(size = 3.5,shape=19)+
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          #legend.title = element_text(colour="black", size=15, face="bold"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  plot(p5)
}
dev.off()

saveRDS(df, file = paste0(colnameslist[1],"&",colnameslist[2],"_WNN_Seurat.rds"))
