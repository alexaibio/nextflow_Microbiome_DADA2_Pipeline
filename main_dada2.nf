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


println """\
         D A D A 2 - N F   P I P E L I N E - @Alex!!! 
         ===========================================
         reads folder     : ${params.reads}
         out folder       : ${params.outdir}

         FW_primers       : ${params.FW_primer}
         RV_primers       : ${params.RV_primer}

         Trunc-For/Rev    : ${params.truncFor}, ${params.truncRev}
         maxEE            : ${params.maxEEFor}, ${params.maxEERev}
         TODO             : add more...
         """
         .stripIndent()


// Create channels with reads: we need two channels because we use it twice. Obsolete in DSL2
Channel
    .fromFilePairs(params.reads, checkIfExists: true)
    .into{ch_read_pairs1; ch_read_pairs2}



// Trim each read-pair with cutadapt
process cut_primers {
    tag " Cutadapt to ${pair_id} "  
    publishDir "${params.outdir}/0_trimmed", mode: 'symlink'

    input:
    tuple val(pair_id), path(reads) from ch_read_pairs1

    output:
    tuple val(pair_id), path("0_trimmed/*.*") into (ch_fastq_trimmed_manifest_1,  ch_fastq_trimmed_manifest_2)
    file "0_trimmed/*.*" into ch_fastq_trimmed_files
    file "cutadapt_log_*.txt" into ch_fastq_cutadapt_log

    script:
    // Note: I moved -g into ${params.FW_primer} itself, a bad desicion, better use fasta file 
    """
    mkdir -p 0_trimmed
    cutadapt ${params.FW_primer} ${params.RV_primer}  \
        -o 0_trimmed/${reads[0]} -p 0_trimmed/${reads[1]} \
        ${reads[0]} ${reads[1]} > cutadapt_log_${pair_id}.txt
    """
}


// it is also possible to use multiqc, not fastqc! 
process initial_quality_fastqc{
    tag " Generate Quality reports for ${pair_id} "  
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


process initial_quality_multiqc {
    tag " MultiQC combined report " 
    publishDir "${params.outdir}/1_multiQC", mode:'copy'  // or symlink?
       
    input:
    path '*' from ch_fastq_cutadapt_log.mix(ch_fastqc_results).collect()
    
    output:
    path 'multiqc_report.html'
     
    script:
    """
    multiqc . 
    """
} 




// ******* START PROCESSING  **********

process filterAndTrim {
    cpus 4

    tag " DADA2: trimming/filtering of ${pair_id} " 
    publishDir "${params.outdir}/2_dada2-FilterAndTrim", mode: "link"

    input:
    set val(pair_id), file(reads) from ch_fastq_trimmed_manifest_2

    output:
    set val(pair_id), "*.R1.filtered.fastq.gz", "*.R2.filtered.fastq.gz" into filteredReads
    file "*.R1.filtered.fastq.gz" into forReads
    file "*.R2.filtered.fastq.gz" into revReads
    file "*.trim_report.csv" into trimTracking  //trimmed.txt

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
                        multithread = ${task.cpus})

    write.csv(out, paste0("${pair_id}", ".trim_report.csv"))
    """
}


// combine all individual trimming reports into one file for convinience
process generateFilteringReport {
    tag " Combining filtering report "
    publishDir "${params.outdir}/2_dada2-FilterAndTrim", mode: "link"

    // operator collects all the items emitted by a channel to a List and return as a sole emission
    input:
    file trimData from trimTracking.collect()

    output:
    file "all.trimmed.csv" into trimmedReadTracking

    script:
    """
    #!/usr/bin/env Rscript
    trimmedFiles <- list.files(path = '.', pattern = '*.trim_report.csv')
    sample.names <- sub('.trim_report.csv', '', trimmedFiles)
    trimmed <- do.call("rbind", lapply(trimmedFiles, function (x) as.data.frame(read.csv(x))))
    colnames(trimmed)[1] <- "Sequence"
    trimmed\$SampleID <- sample.names
    write.csv(trimmed, "all.trimmed.csv", row.names = FALSE)
    """
}



// ERROR calculation
process LearnErrorsForward {
    cpus 4

    tag " DADA2: error rate calculation for Forward reads "
    publishDir "${params.outdir}/3_dada2-LearnErrors", mode: "link"

    input:
    file fReads from forReads.collect()

    output:
    file "errorsF.RDS" into errorsFor
    file "R1.err.pdf"

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2);
    packageVersion("dada2")

    # Learn forward error rates
    filtFs <- strsplit('${fReads}', split = ' ')[[1]]
    errF <- learnErrors(filtFs, nbases=1e8, multithread=${task.cpus}, randomize=TRUE, verbose=1, MAX_CONSIST=20)
    
    pdf(file="R1.err.pdf")
    plotErrors(errF, nominalQ=TRUE)
    dev.off()

    saveRDS(errF, "errorsF.RDS")
    """
}

process LearnErrorsReverse {
    cpus 4

    tag " DADA2: error rate calculation for Reverse reads "
    publishDir "${params.outdir}/3_dada2-LearnErrors", mode: "link"

    input:
    file rReads from revReads.collect()

    output:
    file "errorsR.RDS" into errorsRev
    file "R2.err.pdf"

    script:
    """
    #!/usr/bin/env Rscript
    library(dada2);
    packageVersion("dada2")

    # Learn forward error rates
    filtRs <- strsplit('${rReads}', split = ' ')[[1]]
    errR <- learnErrors(filtRs, nbases=1e8, multithread=${task.cpus}, randomize=TRUE, verbose=1, MAX_CONSIST=20)
    
    pdf(file="R2.err.pdf")
    plotErrors(errR, nominalQ=TRUE)
    dev.off()

    saveRDS(errR, "errorsR.RDS")
    """
}


// Dereplication and sample inference

// for (sam in sample.names)

process LearnErrorsReverse {
    cpus 4

    tag " DADA2: Dereplication and sample inference "
    publishDir "${params.outdir}/4_dada2-sample_inference", mode: "link"

    input:

    output:


    script:
    """
    derepF <- derepFastq(filtFs[[sam]])
    derepR <- derepFastq(filtRs[[sam]])
    
    ### SAMPLE INFERENCE 
    dadaF <- dada(derepF, err=errF, pool=TRUE, multithread = TRUE)
    dadaR <- dada(derepR, err=errR, pool=TRUE, multithread = TRUE)
    
    ### Merge forward and reverse reads
    # NOTE if merger is zero, it means your reads aren't overlapping after truncation - check trimming
    # mergePairs requires 20 nts of overlap by default
    # https://github.com/benjjneb/dada2/issues/419
    merger <- dada2::mergePairs(dadaF, derepF, dadaR, derepR)
    
    if (length(merger$sequence)==0){
        print(" !!!!!  You have a PROBLEM !!!!!! ")
        print("  Forward and Revers reads are not overlapping during merging! Please check trimming parameters!")
        print(" justConcatenate=TRUE has been used for this sample")
        merger <- dada2::mergePairs(dadaF, derepF, dadaR, derepR, justConcatenate=TRUE)  
        }
    mergers[[sam]] <- merger
    
    # add number of merged sequences to a filter.log
    filter.log[sam,"merged"] <- length(merger$sequence)
    
    counter <- counter + 1
    cat("   ::SAMPLE #", counter, "Done...", sam, ":", length(merger$sequence),  " merged sequences.")

    
    
    """

}

