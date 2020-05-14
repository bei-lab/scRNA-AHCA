##-------------------- Cellranger count----------------###
# /data/home/heshuai/software/cellranger-3.0.1/cellranger-cs/3.0.1/bin
# cellranger --version (3.0.1)
# Copyright (c) 2019 10x Genomics, Inc.  All rights reserved.

#BSUB -n 16
#BSUB -o %J.our
#BSUB -e %J.err
#BSUB -R span[hosts=1]
##BSUB -q smp

cellranger count \
--force-cells=13000 \
--id=sample_name \
--localcores=16  \
--transcriptome=/data/home/hanbw/database/refdata-cellranger-GRCh38_and_EBV_5p  \
--fastqs=/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY/AHY725CCXY/sample_name_1,\
/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY/AHY725CCXY/sample_name_2,\
/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY/AHY725CCXY/sample_name_3,\
/data4/heshuai/RAW_data/1-SingleCell/2-NKTCL-SingcleCell/NKTCL_SC1811_YGY_lung/AHY725CCXY/AHY725CCXY/sample_name_4 \
--sample=sample_name

