#!/usr/bin/env bash
# Run FastQC (and MultiQC on the output).
# Peter Hickey
# 2021-03-25

module load fastqc

OUTDIR="../output/FastQC"
mkdir -p ${OUTDIR}

fastqc -o ${OUTDIR} \
       --threads 10 \
       ../extdata/210312_A01221_0027_BHJ7TKDSXY/L002-ChromiumscRNAV31-DI/C133*fastq.gz \
       ../extdata/210312_A01221_0027_BHJ7TKDSXY/L002-ChromiumscRNAV3-ATB/C133*fastq.gz

multiqc --title C133_Neeland \
        --outdir ${OUTDIR} \
        --no-megaqc-upload \
        ${OUTDIR}
