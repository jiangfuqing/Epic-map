#prepara files as iStar recommended: he-raw.jpg, locs-raw.tsv, pixel-size-raw.txt, radius-raw.txt. Set pixel-size-raw=0.46,  radius-raw=54
set -e
device="cuda" 
pixel_size=0.5  
n_genes=1000  

echo $pixel_size > ${prefix} pixel-size.txt
python rescale.py ${prefix} --image
python preprocess.py ${prefix} --image

python extract_features.py ${prefix} --device=${device}

python get_mask.py ${prefix}embeddings-hist.pickle ${prefix}mask-small.png

python select_genes.py --n-top=${n_genes} "${prefix}cnts.tsv" "${prefix}gene-names.txt"

python rescale.py ${prefix} --locs --radius

python impute.py ${prefix} --epochs=400 --device=${device}

python plot_imputed.py ${prefix}

python cluster.py --filter-size=8 --min-cluster-size=20 --n-clusters=10 --mask=${prefix}mask-small.png ${prefix}embeddings-gene.pickle ${prefix}clusters-gene/

python aggregate_imputed.py ${prefix}
python reorganize_imputed.py ${prefix}
python differential.py ${prefix}
python plot_spots.py ${prefix}