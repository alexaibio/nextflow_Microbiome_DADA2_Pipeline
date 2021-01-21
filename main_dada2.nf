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

params.rawfolder = "$baseDir/DATA/raw/"
params.reads = "$baseDir/DATA/raw/*.{R1,R2}.fastq"
params.outdir = "$baseDir/OUT/"
params.FW_primer = "CCAGCAGCYGCGGTAAN" // -g - 5’ adapter  (front), forwardPrimer
params.RV_primer = "CCGTCAATTCNTTTRAGT" //["CCGTCAATTCNTTTRAGT", "CCGTCAATTTCTTTGAGT", "CCGTCTATTCCTTTGANT"] // -G

println """\
         D A D A 2 - N F   P I P E L I N E
         ===================================
         reads            : ${params.reads}
         outdir           : ${params.outdir}
         FW_primer        : ${params.FW_primer}
         RV_primer        : ${params.RV_primer}
         """
         .stripIndent()


ch_read_pairs = Channel.fromFilePairs(params.reads, checkIfExists: true)
//ch_read_pairs.view()



/*
    * Trim each read-pair with cutadapt
*/

process trimming {
    tag "Cutadapt to ${pair_id}"  
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    set val(pair_id), file(reads) from ch_read_pairs

    output:
    set val(pair_id), file("0_trimmed/*.*") into ch_fastq_trimmed_manifest 
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

// head(fnFs) head(fnRs) head(sample.names)
// check if forward and reverse reads match

process quality_fastqc{
    tag "Quality of ${pair_id}"  
    publishDir "${params.outdir}/1_fastQC", mode: 'symlink'

    input:
    set val(pair_id), file(reads) from ch_fastq_trimmed_manifest

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -q ${reads}

    """
}


