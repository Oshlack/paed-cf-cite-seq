#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --mem=15G
#SBATCH --array=1-2
#SBATCH --time=48:00:00
#SBATCH --job-name=C133_cellsnp-lite.sh
#SBATCH --output=cellsnp-lite_%A_%a.out
#SBATCH --partition=regular

# Run cellsnp-lite on capture-specific BAMs
# Peter Hickey
# 2021-05-04

# Setup ------------------------------------------------------------------------

module load cellsnp-lite/1.2.0

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

# Project specific variables ---------------------------------------------------

PROJECT="C133_Neeland"
PROJECT_ROOT="/stornext/Projects/score/Analyses/${PROJECT}"
SEQDIR=${PROJECT_ROOT}/extdata/210312_A01221_0027_BHJ7TKDSXY
CELLRANGERDIR=${SEQDIR}/CellRanger
OUTDIR=${PROJECT_ROOT}/data/cellsnp-lite
mkdir ${OUTDIR}
REGIONSVCF="/stornext/Projects/score/Indexes/cellSNP/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.bgz"

# Run cellsnp-lite  ------------------------------------------------------------

cellsnp-lite --samFile ${CELLRANGERDIR}/C133_${SLURM_ARRAY_TASK_ID}/outs/possorted_genome_bam.bam \
             --outDir ${OUTDIR}/C133_${SLURM_ARRAY_TASK_ID} \
             --regionsVCF ${REGIONSVCF} \
             --barcodeFile ${PROJECT_ROOT}/data/emptyDrops/C133_${SLURM_ARRAY_TASK_ID}.barcodes.txt \
             --nproc 20 \
             --minMAF 0.1 \
             --minCOUNT 20 \
             --gzip
