#!/usr/bin/env nextflow
RUN="/qib/instruments/n121571/Run_11_Non_Adaptive_22EPA077CP_090822/Run_11_Non_Adaptive_22EPA077CP_090822/20220809_0956_MN41075_FAT84700_28fedb04"
nextflow run /qib/instruments/n91636/software/guppy-nf/main.nf \
--fast5 "$RUN" \
--sup_model \
--basecall \
--cuda "1,2,3" \
--chunks_per_runner 600 \
--min_qscore 10 \
--outdir "/data/staging/Run_11_Non_Adaptive" \
--map "${RUN}/map" \
-w "/data/staging/work/Run_11_Non_Adaptive" \
-resume

# Enrich_bacHiT
# SQK-RPB004
# EXP-NBD196
# coronahit