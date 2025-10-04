# Epic-map data analysis pipeline
1. Download and install cellranger-atac software from 10x genomics and replace their default barcode file:
a. Enter cellranger/lib/python/atac/barcodes/;
b. Replace default barcode file "737K-cratac-v1.txt.gz" with new custom barcode (in Barcode, also named as "737K-cratac-v1.txt.gz").

2. Generate ATAC fragments file and gene expression matrix file from MISAR-seq data:
bash Epic-map.sh

3. After fill tissue region with white and non-tissue with black in photoshop or other image processing software manually:
python Figure_filter.py

4. Analysis Epic-map data:
#for Epic-map in small-molecule drug:
Rscript Epic_map_for_drug_Signac.R $sampleId or Rscript Epic_map_for_drug_ArchR.R $sampleId
#for Epic-map in both small-molecule drug and histone modification:
Rscript Epic_map_for_drug&histone.R $sampleId1 $sampleId2 $sampleId1 $modality1 $modality2
