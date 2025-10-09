# Spatial-Epigenome & Chemical drug binding map (Epic-map)
![image](https://github.com/jiangfuqing/Epic-map/blob/main/Epic-map.jpg)

# Epic-map data analysis pipeline
1. Download and install cellranger-atac software from 10x genomics and replace their default barcode file:
a. Enter cellranger-atac-2.0.0/lib/python/atac/barcodes/;
b. Replace default barcode file "737K-cratac-v1.txt.gz" with new custom barcode (in Barcode, also named as "737K-cratac-v1.txt.gz").

2. Generate fragments file and peak_bc_matrix file from Epic-map data:
bash Epic-map.sh

3. Fill tissue region with white and non-tissue with black in photoshop manually and saved as $sampleId-PS.jpg:
python Figure_filter.py -i $sampleId

4. Analysis Epic-map data:
#for Epic-map in small-molecule drug:
Rscript Epic_map_for_drug_Signac.R $sampleId or Rscript Epic_map_for_drug_ArchR.R $sampleId
#for Epic-map in both small-molecule drug and histone modification:
Rscript Epic_map_for_drug&histone.R $sampleId1 $sampleId2 $sampleId1 $modality1 $modality2
