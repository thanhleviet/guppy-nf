#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process guppy_basecaller {
    tag "${fast5.baseName}"
    
    cache true

    publishDir "${params.outdir}",
        mode: "move",
        overwrite: true,
        saveAs: { filename -> filename }

    input:
        path(fast5)

    output:
        path("fastq"), emit: fastq
    script:
        model = params.hac_prom ? "dna_r9.4.1_450bps_hac_prom.cfg " : "dna_r9.4.1_450bps_hac.cfg"
        cuda = params.full_cuda ? "'cuda:all'" : "cuda:0,1"
        """
        guppy_basecaller --input_path $fast5 \\
            -r \\
            --save_path fastq \\
            --records_per_fastq 0 \\
            --compress_fastq \\
            -c ${model} -x ${cuda} \\
            --num_callers ${params.num_callers}
        """
}


process guppy_barcoder {
    publishDir "${params.outdir}/named_samples",
        mode: "move",
        overwrite: true,
        saveAs: { filename -> filename } 
    
    tag {"Running"}
    
    cache true

    cpus 48

    input:
    path("fastq")
    output:
    path("*.fastq.gz")
    script:
    barcode_kit=params.barcode_kit ? "-d  $params.barcoding_folder --barcode_kits '$params.barcode_kit'": ""
    trim_barcodes=params.trim_barcodes ? "--trim_barcodes" : ""
    output=params.trim_barcodes ? "deplexed-barcodes-trimmed" : "deplexed" 
    
    """
    guppy_barcoder -i fastq \\
    -r \\
    -s $output \\
    -q 0 \\
    $barcode_kit \\
    $trim_barcodes \\
    -t ${task.cpus}
    
    # Combine fastq
    cd $output
    for f in barcode*;do echo \$f;cat \$f/*.fastq | pigz -p ${task.cpus} - > ../\${f}.fastq.gz;done
    cd ..
    brename -e -p '(barcode[0-9]{2})' -r '{kv}' -k ${params.map}
    cat <(zcat barcode*.gz) <(cat ${output}/unclassified/*.*) | pigz -p ${task.cpus} - > unclassified.fastq.gz
    rm barcode*.fastq.gz
    """
}


fast5 = Channel.fromPath(params.fast5, type: 'dir', checkIfExists: true)

workflow {
    if (params.basecall) {
        guppy_basecaller(fast5)
        fastq = guppy_basecaller.out.fastq
    } else {
        fastq = fast5
    }
  
  guppy_barcoder(fastq)
}