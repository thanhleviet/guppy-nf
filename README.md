## Command

```
RUN="/data/runs";
nextflow run /qib/instruments/n91636/software/guppy-nf/main.nf \
--basecall \
--hac_prom \
--fast5 "$RUN" \
--outdir /data/staging/ \
--map "${RUN}/map" \
-w "/data/staging/work" \
-resume
```