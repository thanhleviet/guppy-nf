#!/usr/bin/env nextflow
RUN="/data/Rachel_Gilroy_Meta_Run2and3_050922/Run_3/20220905_0916_1-E9-H9_PAK44116_d6bade27/fastq_pass"
nextflow run /qib/instruments/n91636/software/guppy-nf/main.nf \
--fast5 "$RUN" \
--barcode_kit "EXP-NBD196" \
--deplex \
--min_qscore 10 \
--outdir "/data/staging/Rachel_Gilroy_Meta_Run2and3_050922" \
--map "/data/Rachel_Gilroy_Meta_Run2and3_050922/Run_3/20220905_0916_1-E9-H9_PAK44116_d6bade27/fastq_pass/map" \
-w "/data/staging/work/Rachel_Gilroy_Meta_Run2and3_050922" \
-resume
 
# Enrich_bacHiT
# SQK-RPB004
# EXP-NBD196
# coronahit
# EXP-NBD114 Expansion