#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process guppy_basecaller {
    tag "${fast5.baseName}"
    
    cache true

    publishDir "${params.outdir}",
        mode: "copy",
        overwrite: true,
        saveAs: { filename -> filename }

    input:
        path(fast5)

    output:
        path("fastq/pass"), emit: fastq
        tuple val(fast5.baseName), path("fastq/sequencing_summary.txt"), emit: sequencing_summary
    script:
        model = params.sup_model ? "dna_r9.4.1_450bps_sup.cfg" : "dna_r9.4.1_450bps_hac.cfg"
        cuda = params.full_cuda ? "'cuda:all'" : "cuda:0"
        """
        guppy_basecaller --input_path $fast5 \\
            -r \\
            --min_qscore ${params.min_qscore} \\
            --save_path fastq \\
            --records_per_fastq 0 \\
            --chunks_per_runner ${params.chunks_per_runner} \\
            --compress_fastq \\
            -c ${model} -x ${cuda} \\
            --num_callers ${params.num_callers}
        """
}


process guppy_barcoder {
    publishDir "${params.outdir}/named_samples/${output}",
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
    path("*.stats")
    path("${output}/barcoding_summary.txt")
    
    script:
    barcode_kit=params.barcode_kit ? "-d  $params.barcoding_folder --barcode_kits '$params.barcode_kit'": ""
    trim_barcodes=params.trim_barcodes ? "--trim_barcodes" : ""
    output=params.trim_barcodes ? "deplexed-barcodes-trimmed" : "deplexed" 
    barcoder_options = params.barcoder_options ? "$params.barcoder_options" : ""
    """
    guppy_barcoder -i fastq \\
    -r \\
    -s $output \\
    -q 0 \\
    $barcoder_options \\
    $barcode_kit \\
    $trim_barcodes \\
    -t ${task.cpus}
    
    # Combine fastq
    cd $output
    for f in barcode*;do echo \$f;cat \$f/*.fastq | pigz -p ${task.cpus} - > ../\${f}.fastq.gz;done
    cd ..
    if [[ -f '${params.map}' ]]; then
        brename -e -p '(barcode[0-9]{2})' -r '{kv}' -k ${params.map}
        cat <(zcat barcode*.gz) <(cat ${output}/unclassified/*.*) | pigz -p ${task.cpus} - > unclassified.fastq.gz
        rm barcode*.fastq.gz 2> /dev/null
    else
        cat ${output}/unclassified/*.* | pigz -p ${task.cpus} - > unclassified.fastq.gz
    fi
    seqkit stats *.fastq.gz > reads.stats
    """
}

// --require_barcodes_both_ends \\

process STATS_PYCOQC {
    publishDir "${params.outdir}", mode: "copy"
    
    tag {"🏃🏃🏃🏃🏃🏃"}
    
    conda '/qib/instruments/n91636/software/miniconda3/envs/nanopore'

    cpus 8

    input:
    tuple val(run_name), path("sequencing_summary.txt")
    
    output:
    path("${run_name}.html")
    
    script:
    """
    pycoQC -f sequencing_summary.txt -o ${run_name}.html --min_pass_qual ${params.min_qscore}
    """
}


process MINIONQC {
    
    publishDir "${params.outdir}", mode: "copy" 
    
    tag {"🏃🏃🏃🏃🏃🏃"}
    
    conda '/qib/instruments/n91636/software/miniconda3/envs/nanopore/envs/minionqc'

    cpus 8

    input:
        tuple val(run_name), path("sequencing_summary.txt")
    output:
        path("${run_name}_minionqc")
    script:
    """
    MinIONQC.R -i sequencing_summary.txt -o ${run_name}_minionqc -p ${task.cpus} -q ${params.min_qscore}
    """
}

fast5 = Channel.fromPath(params.fast5, type: 'dir', checkIfExists: true)

workflow {
    if (params.basecall) {
        guppy_basecaller(fast5)
        fastq = guppy_basecaller.out.fastq
        STATS_PYCOQC(guppy_basecaller.out.sequencing_summary)
        MINIONQC(guppy_basecaller.out.sequencing_summary)
    } else {
        fastq = fast5
    }
  if (params.deplex) {
      guppy_barcoder(fastq)
  }
  
}