#take JQ1-rep3 as example
#load package
pkgs <- c("ArchR", "Rsamtools", "presto", "Seurat", "grid", "patchwork", "dplyr", "ggplot2", "parallel")
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

# accept args from pbs file
args <- commandArgs(TRUE)

sampleId <- args[1]

color_list<-c("1"="#fcb462","2"="#9500ff","3"="#fb4f2b","4"="#4993b6","5"="#1e25ce","6"="#197a35","7"="#fd8000","8"="#d73027","9"="#8e0054","10"="#bb207c")

#threads default to 1 in windows
addArchRThreads(threads = 50)

set.seed(1234)

#####  data processing #####
addArchRGenome("mm10")
inputFiles <- paste0(sampleId,'-filtered_fragments.tsv.bgz')

set.seed(123)
## Create ArchRProject
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = paste0(sampleId,"_filter"),
  minTSS = 0,
  minFrags = 0,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  force = TRUE,
  TileMatParams = list(tileSize = 5000)
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = paste0(sampleId,"_filter"),
  copyArrows = TRUE
)

## Data normalization and dimensionality reduction
for (i in c(25000, 20000, 15000)){
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix", 
  name = paste0("LSI_", i),
  iterations = 2, 
  clusterParams = list(
    resolution = c(2), 
    sampleCells = 10000, 
    n.start = 10
  ),
  varFeatures = i,
  dimsToUse = 1:30,
  force = TRUE
)

for (j in c(0.3, 0.5, 0.8)){
proj <- addClusters(
  input = proj, 
  method = "seurat",
  reducedDims = paste0("LSI_", i),
  name = paste0("Clusters_ArchR_",j,"_",i),
  resolution = j,
  prefix = "",
  force = TRUE
)
}

proj <- addUMAP(
  ArchRProj = proj, 
  reducedDims = paste0("LSI_", i),
  name = paste0("UMAP_", i), 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)}

#read the spatial barcode
barcode <- read.csv("~/barcode/Epic_barcode.csv",header = T,row.names = 1)

#read and filter singlecell.csv file which from cellranger

##### data filtering #####
#filter spot not in tissue, generate position.txt file from matlab.
location <- read.table(paste0("position_",sampleId,".txt"), sep =",", header = FALSE, dec =".", stringsAsFactors = F)
x <- as.character(location[1,])
x = x[-1]
x <- as.data.frame(x)
colnames(x) <- "array"
y <- read.csv("~/barcode/Epic_barcode_filter.csv",header = T) 
z <- merge(x,y,by="array")
z <- z[,-1]

barcode <- barcode[z,]

#read and filter singlecell.csv file which from cellranger
singelcell <- read.csv(paste0(sampleId,"-singlecell.csv"),header = T,row.names = 1)

singelcell <- merge(singelcell,barcode,by="row.names")

rownames(singelcell) <- paste0(sampleId,"_filter#",singelcell$Row.names)

singelcell<-transform(singelcell, percent_mito=(singelcell$mitochondrial/singelcell$total)*100, 
                      FRiP=(singelcell$peak_region_fragment/singelcell$passed_filters)*100, 
                      percent_TSS_fragments=(singelcell$TSS_fragments/singelcell$passed_filters)*100,
                      total_fragments=singelcell$passed_filters,Sample=sampleId)

QC_index <- c("array_col","array_row","percent_mito","FRiP","percent_TSS_fragments","total_fragments")

for (i in QC_index) {
  proj <- addCellColData(ArchRProj = proj, data = singelcell[,i], cells = rownames(singelcell),name = i, force = T)
}

# plot atac spatial feature
#imported_raster=OpenImageR::readImage(paste0(sampleId,"-HE.jpg"))
#g <- rasterGrob(imported_raster, width=unit(1,"npc"), height=unit(1,"npc"), interpolate = FALSE)

pdf(file = paste0(sampleId,"_Spatial_QC_ArchR.pdf"), width=8.6, height=8.6)

df <- as.data.frame(proj@cellColData)
df <- transform(df, log10_unique_fragments=log10(df$nFrags))

for (i in c("TSSEnrichment","nFrags","log10_unique_fragments","total_fragments","FRiP","percent_mito","percent_TSS_fragments")){
  
  xlab=mean(df[,i])
  ylab=paste0("Filter_",i)
  
  p1 <- ggplot(df,aes(x=Sample,y=df[,i],fill=Sample)) +
    geom_violin(alpha=0.8,width=1) + 
    guides(fill=F)+xlab(paste0("mean=",round(xlab,2))) + 
    ylab(ylab) + 
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"), 
          axis.text=element_text(size=20), axis.title=element_text(size=20,face="bold"), 
          panel.background = element_blank())
  
  p2 <- ggplot(df, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df[,i])) + scale_color_gradientn(colours = ArchRPalettes$solarExtra) +
    ggtitle(ylab) +
    #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(size = 1,shape=19)+
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  plot(p1|p2)
}
dev.off()

# plot atac spatial cluster
proj$merge_cluster <- proj$Clusters_ArchR_0.8_20000

clusters <- "merge_cluster"
pdf(file = paste0(sampleId,"_Spatial_Cluster_ArchR.pdf"), width=8.6, height=8.6)

for (i in clusters){
  p3 <- plotEmbedding(proj, name = i, embedding = "UMAP_20000", rastr = FALSE, plotAs = "points", pal = color_list[1:length(unique(df[,i]))], size = 0.5, labelAsFactors=F, labelMeans=F)

  p4 <- ggplot(df, aes(x = as.numeric(array_col), y = as.numeric(array_row), color=df[,i])) +
    scale_color_manual(values = color_list[sort(as.numeric(unique(df[,i])))],
                       labels=paste0("",sort(as.numeric(unique(df[,i]))))) +
    ggtitle(i) +
    #annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    geom_point(size = 1,shape=19) +
    expand_limits(x = 0, y = 0) +
    scale_x_continuous(name="X", limits = c(NA, NA), expand = expansion(mult = c(-0.013, -0.013))) +
    scale_y_reverse(name="Y", limits = c(NA, NA), expand = expansion(mult = c(-0.013, 0.008))) +
    coord_equal(xlim=c(0,51),ylim=c(51,1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20,face="bold"),
          legend.text=element_text(size=20),
          legend.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  
  plot(p3|p4)}

dev.off()

saveRDS(proj,paste0(sampleId,"_ArchR.rds"))
