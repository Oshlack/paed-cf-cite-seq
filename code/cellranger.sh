#!/bin/bash
#SBATCH --partition=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --job-name=C133_cellranger.sh

# Run CellRanger
# Peter Hickey
# 2020-03-23

# NOTE: This was sample was captured over 2 GEM wells (C133_1, C133_2) and
#       sequenced over 1 Illumina run:
#
#       1. 210312_A01221_0027_BHJ7TKDSXY: GEX, ADT, and HTO of C133_1, C133_2
#          @ MCRI on NovaSeq.

# Setup ------------------------------------------------------------------------

module load cellranger/5.0.0
module load bcl2fastq/2.19.1

# Project specific variables ---------------------------------------------------

PROJECT="C133_Neeland"
PROJECT_ROOT="/stornext/Projects/score/Analyses/${PROJECT}"
SEQDIRS=( ${PROJECT_ROOT}/extdata/210312_A01221_0027_BHJ7TKDSXY )
# NOTE: $RUNDIRS not required because data have already been demultiplexed.

LIBRARIESCSV=${PROJECT_ROOT}/data/sample_sheets/libraries.csv
AGGRCSV=${PROJECT_ROOT}/data/sample_sheets/${PROJECT}.libraries.csv

FEATUREREFCSV=${PROJECT_ROOT}/data/sample_sheets/features.csv

TRANSCRIPTOME="/stornext/Projects/score/Indexes/refdata-gex-GRCh38-2020-A"

CAPTURES=( C133_1 C133_2 )

# General variables ------------------------------------------------------------

OUTDIR=${PROJECT_ROOT}/output/CellRanger
mkdir -p ${OUTDIR}
# NOTE: All `cellranger count` and `cellranger aggr` output goes into directory
#       for the first (and only) sequencing run.
CELLRANGERDIR=${SEQDIRS[0]}/CellRanger
mkdir -p ${CELLRANGERDIR}

# Generate FASTQs with cellranger mkfastq  -------------------------------------

# NOTE: 210312_A01221_0027_BHJ7TKDSXY has already been demultiplexed from the
#       BCLs to FASTQs by MCRI.

# Single-library analyses with cellranger count --------------------------------

# NOTE: CellRanger writes a bunch of output to $(pwd) so have to move to the
#       directory where you want the output to go. Sigh.
cd ${CELLRANGERDIR} || exit
for CAPTURE in "${CAPTURES[@]}"
do
  # NOTE: Extract the relevant FASTQ files from the libraries.csv for this
  #       particular capture.
  echo ${CAPTURE}
  SAMPLE=$(echo $CAPTURE | cut -d '_' -f 2)
  head -n 1 ${LIBRARIESCSV} > ${CAPTURE}.libraries.csv
  grep "C133-MN-GEX-${SAMPLE}" ${LIBRARIESCSV} >> ${CAPTURE}.libraries.csv
  grep "C133-MN-ADT-${SAMPLE}" ${LIBRARIESCSV} >> ${CAPTURE}.libraries.csv
  grep "C133-MN-HTO-${SAMPLE}" ${LIBRARIESCSV} >> ${CAPTURE}.libraries.csv

  # NOTE: ADT+HTO libraries in L002-ChromiumscRNAV3-ATB were sequenced with a
  #       75bp R1, so need to truncate to 28bp.
  cellranger count --id=${CAPTURE} \
                   --transcriptome=${TRANSCRIPTOME} \
                   --libraries=${CAPTURE}.libraries.csv \
                   --feature-ref=${FEATUREREFCSV} \
                   --localcores=20 \
                   --localmem=100 \
                   --r1-length=28

  cp ${CAPTURE}/outs/web_summary.html \
     ${OUTDIR}/${CAPTURE}.web_summary.html
  cp ${CAPTURE}/outs/cloupe.cloupe ${OUTDIR}/${CAPTURE}.cloupe
done

# Aggregating Multiple GEM Wells with cellranger aggr --------------------------

# NOTE: CellRanger writes a bunch of output to $(pwd) so have to move to the
#       directory where you want the output to go. Sigh.
cd ${CELLRANGERDIR} || exit

# No normalization
cellranger aggr --id=${PROJECT} \
                --csv=${AGGRCSV} \
                --normalize=none
cp ${PROJECT}/outs/web_summary.html \
   ${OUTDIR}/${PROJECT}.web_summary.html
cp ${PROJECT}/outs/count/cloupe.cloupe ${OUTDIR}/${PROJECT}.cloupe

# Read-depth normalization
cellranger aggr --id=${PROJECT}_read_depth_normalized \
                --csv=${AGGRCSV} \
                --normalize=mapped
cp ${PROJECT}_read_depth_normalized/outs/web_summary.html \
   ${OUTDIR}/${PROJECT}_read_depth_normalized.web_summary.html
cp ${PROJECT}_read_depth_normalized/outs/count/cloupe.cloupe ${OUTDIR}/${PROJECT}_read_depth_normalized.cloupe
