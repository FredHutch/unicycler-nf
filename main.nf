#!/usr/bin/env nextflow

// Default values
params.min_fasta_length = 100

def helpMessage() {
    log.info"""
    Usage:

    nextflow run fredhutch/unicycler-nf <ARGUMENTS>
    
    Arguments:
      --sample_sheet          CSV file listing samples to analyze
      --output_folder      Folder to place outputs

    Options:
      --mode               Either 'normal' (default), 'conservative', or 'bold'
      --min_fasta_length   Minimum contig length (default: 100)

    Sample Sheet:
      The sample_sheet is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `name` and a column `fastq`.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Read in the sample sheet
Channel.from(file(params.sample_sheet))
    .splitCsv(header: true, sep: ",")
    .map { sample ->
        tuple(sample["name"], file(sample["fastq"]))}
    .set{ input_ch }


// Run the assembly
process unicycler {
    container "quay.io/biocontainers/unicycler:0.4.7--py37hdbcaa40_1"
    cpus 16
    memory "128 GB"
    // errorStrategy "retry"
    publishDir "${output_folder}/"

    input:
    val threads from 16
    val min_fasta_length from params.min_fasta_length
    set val(genome_name), file(fastq) from input_ch

    output:
    file "${genome_name}/${genome_name}.gfa"
    file "${genome_name}/${genome_name}.fasta"
    file "${genome_name}/${genome_name}.log"

    afterScript "rm -rf *"

"""
set -e

mkdir ${genome_name}

unicycler \
    -l ${fastq} \
    -o ${genome_name} \
    --min_fasta_length ${min_fasta_length} \
    --keep 0 \
    -t ${threads}

mv ${genome_name}/assembly.gfa ${genome_name}/${genome_name}.gfa
mv ${genome_name}/assembly.fasta ${genome_name}/${genome_name}.fasta
mv ${genome_name}/unicycler.log ${genome_name}/${genome_name}.log
"""
}
