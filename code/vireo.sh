#!/bin/bash
#SBATCH --cpus-per-task=20
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --array=1-2
#SBATCH --time=48:00:00
#SBATCH --job-name=C133_vireo.sh
#SBATCH --output=vireo_%A_%a.out
#SBATCH --partition=regular

# Run vireo on capture-specific BAMs
# Peter Hickey
# 2021-05-05

# Setup ------------------------------------------------------------------------

module load vireoSNP/0.5.6

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

# Project specific variables ---------------------------------------------------

PROJECT="C133_Neeland"
PROJECT_ROOT="/stornext/Projects/score/Analyses/${PROJECT}"
CELLSNPDIR=${PROJECT_ROOT}/data/cellsnp-lite
OUTDIR=${PROJECT_ROOT}/data/vireo
mkdir ${OUTDIR}

# Run vireo  -------------------------------------------------------------------

vireo --cellData=${CELLSNPDIR}/C133_${SLURM_ARRAY_TASK_ID} \
      --nDonor=8 \
      --outDir=${OUTDIR}/C133_${SLURM_ARRAY_TASK_ID} \
      --nInit=200 \
      --nproc=20
