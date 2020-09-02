for i in $(ls /data4/heshuai/RAW_data/1-SingleCell/3-HCA/2-Mapping | grep cDNA)
do
sed "s/tissue_type/$i/g" /data4/heshuai/RAW_data/1-SingleCell/3-HCA/3-analysis/0-Raw_Seurat_cluster/Single_tissue_clustering.lsf|bsub
done
