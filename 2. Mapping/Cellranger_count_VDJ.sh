##----------------------- Cellranger vdj --------------####
# /data/home/heshuai/software/cellranger-3.1.0/cellranger-cs/3.1.0/bin
# cellranger --version (3.1.0)
# Copyright (c) 2019 10x Genomics, Inc.  All rights reserved.

#BSUB -n 16
#BSUB -o %J.our
#BSUB -e %J.err
#BSUB -R span[hosts=1]
##BSUB -q smp

cellranger vdj \
--id=sample_name \
--localcores=16  \
--reference=/data/home/heshuai/reference_data/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0 \
--fastqs=/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY_SS1_SCC3_BCR_reversed/AHY725CCXY/sample_name_1,\
/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY_SS1_SCC3_BCR_reversed/AHY725CCXY/sample_name_2,\
/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY_SS1_SCC3_BCR_reversed/AHY725CCXY/sample_name_3,\
/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY_SS1_SCC3_BCR_reversed/AHY725CCXY/sample_name_4 \
--sample=sample_name
