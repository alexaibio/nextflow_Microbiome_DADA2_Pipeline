#!/usr/bin/env nextflow


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run nf-core/mag --input '*_R{1,2}.fastq.gz' -profile docker
    nextflow run nf-core/mag --manifest 'manifest.tsv' -profile docker
    
    Mandatory arguments:
      --input [file]                  Path to input data (must be surrounded with quotes)
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated)
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more
      --FW_primer [str]             Forward primer sequence
	  --RV_primer [str]             Reverse primer sequence

    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// you can move it to config
params.reads = "$baseDir/DATA/raw/*.{R1,R2}.fastq"
params.outdir = "$baseDir/OUT/"






println """\
         D A D A 2 - N F   P I P E L I N E
         ===================================
         reads            : ${params.reads}
         outdir           : ${params.outdir}
         FW_primer        : ${params.FW_primer}
         RV_primer        : ${params.RV_primer}
         """
         .stripIndent()


//ch_read_pairs = Channel.fromFilePairs(params.reads, checkIfExists: true)
// we need two channels because we use it twice
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .into{ch_read_pairs1; ch_read_pairs2}



/*
* Trim each read-pair with cutadapt
* TODO : add all adapters!
*/

process trimming {
    tag " Cutadapt to ${pair_id} "  
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    tuple val(pair_id), path(reads) from ch_read_pairs1

    output:
    tuple val(pair_id), path("0_trimmed/*.*") into (ch_fastq_trimmed_manifest_1,  ch_fastq_trimmed_manifest_2)
    file "0_trimmed/*.*" into ch_fastq_trimmed_files
    //file "cutadapt_log_*.txt" into ch_fastq_cutadapt_log

    script:
    """
    mkdir -p 0_trimmed
    cutadapt -g ${params.FW_primer} -G ${params.RV_primer}  \
        -o 0_trimmed/${reads[0]} -p 0_trimmed/${reads[1]} ${reads[0]} ${reads[1]} > cutadapt_log_${pair_id}.txt
    """
}

//ch_fastq_trimmed_manifest.view { "::::::: Manifest: " + it }


process quality_fastqc{
    tag " Quality of ${pair_id} "  
    publishDir "${params.outdir}/1_fastQC", mode: 'symlink'

    input:
    tuple val(pair_id), path(reads) from ch_fastq_trimmed_manifest_1

    output:
    path "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -q ${reads}

    """
}




process filterAndTrim {
    tag " DADA trimming of ${pair_id} " 
    publishDir "${params.outdir}/2_dada2-FilterAndTrim", mode: "link"

    input:
    set val(pair_id), file(reads) from ch_fastq_trimmed_manifest_2

    output:
    set val(pair_id), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" into filteredReads
    file "*.R1.filtered.fastq.gz" into forReads
    file "*.R2.filtered.fastq.gz" into revReads
    file "*.trimmed.txt" into trimTracking

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2); packageVersion("dada2")

    out <- filterAndTrim(fwd = "${reads[0]}",
                        filt = paste0("${pair_id}", ".R1.filtered.fastq.gz"),
                        rev = "${reads[1]}",
                        filt.rev = paste0("${pair_id}", ".R2.filtered.fastq.gz"),
                        trimLeft = c(0,0),
                        truncLen = c(${params.truncFor},${params.truncRev}),
                        maxEE = c(${params.maxEEFor},${params.maxEERev}),
                        truncQ = ${params.truncQ},  
                        maxN = ${params.maxN},                      
                        compress = TRUE,
                        verbose = TRUE,
                        multithread = 2)

    write.csv(out, paste0("${pair_id}", ".trimmed.txt"))
    """
}


// ERROR calculation
