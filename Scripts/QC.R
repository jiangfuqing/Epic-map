pkgs <- c("ArchR","tidyverse","Matrix","patchwork","grid","magrittr","viridis","igraph","ggalluvial") 
suppressWarnings(suppressPackageStartupMessages(lapply(pkgs,library,character.only=TRUE)))

addArchRGenome("mm10")
addArchRThreads(threads = 50)

#proj <- readRDS(paste0(omics,"_proj_all.rds"))
omics <- "Single" # "Double", "Single" 

sampleId <- c("DOX-rep1","DOX-rep2","DOX-rep3","DOX-20um-rep1","DOX-20um-rep2","JQ1-rep1","JQ1-rep2","JQ1-rep3","BRD4","THZ1-rep1","THZ1-rep2","THZ1-rep3")
input_ATAC <- c()
for (i in sampleId){input_ATAC[i] <- paste0(i,"/",i,'_filtered_fragments.tsv.bgz')}

ArrowFiles <- createArrowFiles(
  inputFiles = input_ATAC,
  sampleNames = sampleId,
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
  outputDirectory = "Single",
  #copyArrows = TRUE
)

metrics <- list()
for(i in sampleId){metrics[[i]] = read.csv(paste0(i,"/singlecell.csv"), header = T, row.names = 1)
                   metrics[[i]] <- metrics[[i]][-1,]
                   rownames(metrics[[i]]) <- paste0(i,"#",rownames(metrics[[i]]))}

singlecell <- bind_rows(metrics)

singlecell <- singlecell[rownames(proj@cellColData),]

singlecell <- transform(singlecell, percent_mito=(singlecell$mitochondrial/singlecell$total)*100, 
                      FRiP=(singlecell$peak_region_fragment/singlecell$passed_filters)*100, 
                      percent_TSS_fragments=(singlecell$TSS_fragments/singlecell$passed_filters)*100,
                      total_fragments=singlecell$passed_filters)
                      
QC_index <- c("percent_mito","FRiP","percent_TSS_fragments","total_fragments")

for (i in QC_index) {
  proj <- addCellColData(ArchRProj = proj, data = singlecell[,i], cells = rownames(singlecell), name = i, force = T)
}

p1 <- plotTSSEnrichment(proj)
p2 <- plotFragmentSizes(proj)
plotPDF(p1, p2, name = paste0(omics,"_All_QC1"), addDOC = FALSE, width = 4, height = 4)


pdf(paste0(omics,"_All_QC2.pdf"), width= 5.6, height = 4.6)

colors <- rainbow(length(sampleId))
data <- data.frame(proj@cellColData)
data <- transform(data, log10_unique_fragments=log10(data$nFrags))

#reorder the sample
data$Sample <- factor(data$Sample, levels = sampleId)

for (i in c("TSSEnrichment","log10_unique_fragments","FRiP","percent_mito","percent_TSS_fragments")){
p <- ggplot(data, mapping = aes(x = Sample, y = !!sym(i))) +
  geom_violin(aes(fill = Sample), trim = FALSE) +
  geom_point(size = 0) +
  #ylim(0, 10) +
  geom_boxplot(width = 0.2, linewidth = 0, outlier.shape = NA) +
  scale_fill_manual(values = colors) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot(p)
}
dev.off()