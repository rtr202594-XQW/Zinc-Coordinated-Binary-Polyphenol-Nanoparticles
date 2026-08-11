#!/bin/bash

#########################################################
# Mouse RNA-seq analysis pipeline
# RSEM + Bowtie2
# Genome: Mus musculus GRCm39
#########################################################


############################
############################

PROJECT="/home/lxx/mouse_RNAseq"
RAW="${PROJECT}/01_raw_data"
REF="${PROJECT}/02_reference/Mus_musculus.GRCm39.dna.toplevel.fa"
GTF="${PROJECT}/02_reference/Mus_musculus.GRCm39.114.gtf"
INDEX_DIR="${PROJECT}/03_rsem_index"
QUANT_DIR="${PROJECT}/04_rsem_quant"
MATRIX_DIR="${PROJECT}/05_matrix"
THREADS=32
############################
# 样本列表
############################
SAMPLES=(
TA-1
TA-2
TA-3
TA-4
TB-1
TB-2
TB-3
TB-4
TG-1
TG-2
TG-3
TG-4
)



############################
############################
mkdir -p ${INDEX_DIR}
mkdir -p ${QUANT_DIR}
mkdir -p ${MATRIX_DIR}

echo "======================================"
echo "Mouse RNA-seq pipeline start"
echo "Time: $(date)"
echo "======================================"



############################
# 建立RSEM索引
############################

RSEM_INDEX="${INDEX_DIR}/mm39_mouse"



if [ ! -f "${RSEM_INDEX}.grp" ]

then


echo "======================================"
echo "Building RSEM reference index"
echo "======================================"


rsem-prepare-reference \
--gtf ${GTF} \
--bowtie2 \
${REF} \
${RSEM_INDEX}


else

echo "RSEM index already exists, skip"

fi





############################
# RSEM表达定量
############################


echo "======================================"
echo "Start RSEM quantification"
echo "======================================"


for SAMPLE in ${SAMPLES[@]}

do


echo "--------------------------------------"
echo "Processing ${SAMPLE}"
echo "--------------------------------------"



R1="${RAW}/${SAMPLE}_1.fq.gz"

R2="${RAW}/${SAMPLE}_2.fq.gz"



if [ ! -f ${R1} ]

then

echo "Cannot find ${R1}"
exit 1

fi



if [ ! -f ${R2} ]

then

echo "Cannot find ${R2}"
exit 1

fi



rsem-calculate-expression \
-p ${THREADS} \
--paired-end \
--bowtie2 \
${R1} \
${R2} \
${RSEM_INDEX} \
${QUANT_DIR}/${SAMPLE}



echo "${SAMPLE} finished"

done




############################
# 合并表达矩阵
############################


echo "======================================"
echo "Generate expression matrix"
echo "======================================"


find ${QUANT_DIR} \
-name "*.genes.results" \
> ${MATRIX_DIR}/quant_files.txt



perl /home/lxx/miniconda3/envs/rna-seq/bin/abundance_estimates_to_matrix.pl \
--est_method RSEM \
--gene_trans_map none \
--cross_sample_norm none \
--out_prefix ${MATRIX_DIR}/mouse_RNAseq \
--quant_files ${MATRIX_DIR}/quant_files.txt



echo "======================================"
echo "ALL DONE"
echo "Finished: $(date)"
echo "======================================"