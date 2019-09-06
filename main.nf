#!/usr/bin/env nextflow

// Default values
params.min_fasta_length = 100
params.long_reads = false
params.short_reads = false
params.prokka = false

def helpMessage() {
    log.info"""
    Usage:

    nextflow run fredhutch/unicycler-nf <ARGUMENTS>
    
    Arguments:
      --sample_sheet       CSV file listing samples to analyze
      --output_folder      Folder to place outputs

    Options:
      --short_reads        Sample sheet contains short read data (`short_R1` and `short_R2`)
      --long_reads         Sample sheet contains long read data (`long_reads`)
      --min_fasta_length   Minimum contig length (default: 100)
      --prokka             Run the genome annotation tool Prokka
      --help               Display this message

    Sample Sheet:
      The sample_sheet is a CSV with a header indicating which samples correspond to which files.
      The file must contain the column `name`, and `long_reads`, `short_R1`, `short_R2` as appropriate.

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

if (!params.long_reads && !params.short_reads){
    log.info"""
###################################################################
Please specify either --long_reads or --short_reads as appropriate.
###################################################################
"""
    helpMessage()
    exit 0
}

// Run UniCycler with each of three different assembly modes
unicycler_modes = ['conservative', 'normal', 'bold']

// Use the --short_reads and --long_reads flags to determine what processes to run

if (params.short_reads){
    // Short reads only
    Channel.from(file(params.sample_sheet)
        .splitCsv(header: true, sep: ","))
        .map { sample ->
            tuple(sample["name"], file(sample["short_R1"]), file(sample["short_R2"]))}
        .set{ illumina_ch }

    // Run the assembly
    process unicyclerShortReadsOnly {
        container "quay.io/biocontainers/unicycler:0.4.7--py37hdbcaa40_1"
        cpus 16
        memory "128 GB"
        errorStrategy "retry"
        publishDir "${params.output_folder}/${genome_name}/${read_type}/${mode}/"

        input:
        val threads from 16
        val read_type from "short_reads"
        val min_fasta_length from params.min_fasta_length
        set val(genome_name), file(short_R1), file(short_R2) from illumina_ch
        each mode from unicycler_modes

        output:
        file "${genome_name}/${genome_name}.${read_type}.${mode}.gfa"
        file "${genome_name}/${genome_name}.${read_type}.${mode}.log"
        set val("${genome_name}"), file("${genome_name}/${genome_name}.${read_type}.${mode}.fasta.gz"), val("${mode}"), val("short_reads") into contigsShortReadsOnly
        

        afterScript "rm -rf *"

    """
    set -e

    mkdir ${genome_name}

    unicycler \
        -1 ${short_R1} \
        -2 ${short_R2} \
        -o ${genome_name} \
        --min_fasta_length ${min_fasta_length} \
        --keep 0 \
        --mode ${mode} \
        -t ${threads}

    mv ${genome_name}/assembly.gfa ${genome_name}/${genome_name}.${read_type}.${mode}.gfa
    mv ${genome_name}/assembly.fasta ${genome_name}/${genome_name}.${read_type}.${mode}.fasta
    mv ${genome_name}/unicycler.log ${genome_name}/${genome_name}.${read_type}.${mode}.log
    gzip ${genome_name}/${genome_name}.${read_type}.${mode}.fasta
    """
    }
}

if (params.long_reads){
    // Long reads only
    Channel.from(file(params.sample_sheet)
        .splitCsv(header: true, sep: ","))
        .map { sample ->
            tuple(sample["name"], file(sample["long_reads"]))}
        .set{ long_read_ch }

    // Run the assembly
    process unicyclerLongReadsOnly {
        container "quay.io/biocontainers/unicycler:0.4.7--py37hdbcaa40_1"
        cpus 16
        memory "128 GB"
        errorStrategy "retry"
        publishDir "${params.output_folder}/${genome_name}/${read_type}/${mode}/"

        input:
        val threads from 16
        val read_type from "long_reads"
        val min_fasta_length from params.min_fasta_length
        set val(genome_name), file(long_reads) from long_read_ch
        each mode from unicycler_modes

        output:
        file "${genome_name}/${genome_name}.${read_type}.${mode}.gfa"
        set val("${genome_name}"), file("${genome_name}/${genome_name}.${read_type}.${mode}.fasta.gz"), val("${mode}"), val("long_reads") into contigsLongReadsOnly
        file "${genome_name}/${genome_name}.${read_type}.${mode}.log"

        afterScript "rm -rf *"

    """
    set -e

    mkdir ${genome_name}

    unicycler \
        -l ${long_reads} \
        -o ${genome_name} \
        --min_fasta_length ${min_fasta_length} \
        --keep 0 \
        --mode ${mode} \
        -t ${threads}

    mv ${genome_name}/assembly.gfa ${genome_name}/${genome_name}.${read_type}.${mode}.gfa
    mv ${genome_name}/assembly.fasta ${genome_name}/${genome_name}.${read_type}.${mode}.fasta
    mv ${genome_name}/unicycler.log ${genome_name}/${genome_name}.${read_type}.${mode}.log
    gzip ${genome_name}/${genome_name}.${read_type}.${mode}.fasta
    """
    }

    if (params.short_reads){
        // Hybrid assembly takes all input data types
        Channel.from(file(params.sample_sheet)
            .splitCsv(header: true, sep: ","))
            .map { sample ->
                tuple(sample["name"], file(sample["long_reads"]), file(sample["short_R1"]), file(sample["short_R2"]))}
            .set{ hybrid_ch }

        // Run the assembly
        process unicyclerHybrid {
            container "quay.io/biocontainers/unicycler:0.4.7--py37hdbcaa40_1"
            cpus 16
            memory "128 GB"
            errorStrategy "retry"
            publishDir "${params.output_folder}/${genome_name}/${read_type}/${mode}/"

            input:
            val threads from 16
            val read_type from "hybrid"
            val min_fasta_length from params.min_fasta_length
            set val(genome_name), file(long_reads), file(short_R1), file(short_R2) from hybrid_ch
            each mode from unicycler_modes

            output:
            file "${genome_name}/${genome_name}.${read_type}.${mode}.gfa"
            set val("${genome_name}"), file("${genome_name}/${genome_name}.${read_type}.${mode}.fasta.gz"), val("${mode}"), val("hybrid") into contigsHybrid
            file "${genome_name}/${genome_name}.${read_type}.${mode}.log"

            afterScript "rm -rf *"

        """
        set -e

        mkdir ${genome_name}

        unicycler \
            -1 ${short_R1} \
            -2 ${short_R2} \
            -l ${long_reads} \
            -o ${genome_name} \
            --min_fasta_length ${min_fasta_length} \
            --keep 0 \
            --mode ${mode} \
            -t ${threads}

        mv ${genome_name}/assembly.gfa ${genome_name}/${genome_name}.${read_type}.${mode}.gfa
        mv ${genome_name}/assembly.fasta ${genome_name}/${genome_name}.${read_type}.${mode}.fasta
        mv ${genome_name}/unicycler.log ${genome_name}/${genome_name}.${read_type}.${mode}.log
        gzip ${genome_name}/${genome_name}.${read_type}.${mode}.fasta
        """
        }
    }
}

// Combine all channels into one
if (params.short_reads && params.long_reads){
    contigsShortReadsOnly
    .mix(contigsLongReadsOnly,contigsHybrid)
    .into { contigsForSummary; contigsForProkka }
} else {
    if (params.short_reads){
        contigsShortReadsOnly
        .into { contigsForSummary; contigsForProkka }
    } else {
        contigsLongReadsOnly
        .into { contigsForSummary; contigsForProkka }
    }
}

// Summarize the assembly
process summarizeAssemblies {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    cpus 1
    memory "4 GB"
    errorStrategy "retry"
    publishDir "${params.output_folder}/${genome_name}/${read_type}/${mode}/"

    input:
    set val(genome_name), file(contigs_fasta_gz), val(mode), val(read_type) from contigsForSummary

    output:
    file "${genome_name}.${read_type}.${mode}.json" into assembly_summaries_ch

    afterScript "rm -rf *"

"""
#!/usr/bin/env python3
import gzip
import json
import pandas as pd

# Count the size, depth, and circularity of each contig
contig_info = []
contig_lengths = dict()
length_buffer = 0
contig_name = None
for l in gzip.open("${contigs_fasta_gz}", "rt"):
    if l.startswith(">"):
        if contig_name is not None:
            contig_lengths[contig_name] = length_buffer
            length_buffer = 0
        if " " in l:
            contig_name, contig_dict = l.lstrip(">").rstrip("\\n").split(" ", 1)
            contig_dict = dict([
                (i.split("=")[0], i.split("=")[1])
                for i in contig_dict.split(" ")
            ])
            contig_dict["name"] = contig_name
            contig_dict["circular"] = contig_dict.get("circular", "false")
            contig_info.append(contig_dict)
        else:
            contig_name = l.lstrip(">").rstrip("\\n")
    else:
        length_buffer += len(l.rstrip("\\n"))
# Add the final contig
contig_lengths[contig_name] = length_buffer

# Make into a DataFrame
if len(contig_info) > 0:
    contig_info = pd.DataFrame(contig_info)
    contig_info["length"] = contig_info["length"].apply(int)
    contig_info["depth"] = contig_info["depth"].apply(lambda s: float(s.rstrip("x")))
    contig_info["circular"] = contig_info["circular"].fillna("false") == "true"
else:
    contig_info = pd.DataFrame(dict([("length", contig_lengths)])).reset_index()
    contig_info["depth"] = 1
    contig_info["circular"] = False
contig_info.sort_values(by="length", ascending=False, inplace=True)

# Calculate N50
running_total = 0
n50 = None
for nbases in contig_info["length"].values:
    running_total += nbases
    if running_total >= contig_info["length"].sum() / 2.:
        n50 = int(nbases)
        break
assert n50 is not None

# Summarize these contigs
output = {
    "mode": "${mode}",
    "read_type": "${read_type}",
    "genome_name": "${genome_name}",
    "num_contigs": int(contig_info.shape[0]),
    "num_circular_contigs": int(contig_info["circular"].sum()),
    "circular_contig_lengths": ", ".join(contig_info.query("circular")["length"].apply(str).tolist()),
    "linear_contig_lengths": ", ".join(contig_info.query("circular == False")["length"].apply(str).tolist()),
    "num_bases": int(contig_info["length"].sum()),
    "longest_contig": int(contig_info["length"].max()),
    "num_over_1Mb": int((contig_info["length"] >= 1000000).sum()),
    "num_100kb_to_1Mb": int(((contig_info["length"] < 1000000) & (contig_info["length"] >= 100000)).sum()),
    "num_10kb_to_100kb": int(((contig_info["length"] < 100000) & (contig_info["length"] >= 10000)).sum()),
    "num_1kb_to_10kb":    int(((contig_info["length"] < 10000) & (contig_info["length"] >= 1000)).sum()),
    "num_under_1kb":        int((contig_info["length"] < 1000).sum()),
    "N50": n50,
}

json.dump(output, open("${genome_name}.${read_type}.${mode}.json", "wt"), indent=4)

"""
}

if (params.prokka) {
    process prokkaAnnotate {
        container 'quay.io/fhcrc-microbiome/prokka@sha256:bfa2e2c7cfb5d010cc488f3f8ed8519762f0af61d0d7b32b83ac6bac1eed22a4'
        cpus 16
        memory "120 GB"
        errorStrategy "retry"
        publishDir "${params.output_folder}/${genome_name}/${read_type}/${mode}/"

        input:
        set val(genome_name), file(contigs_fasta_gz), val(mode), val(read_type) from contigsForProkka
        val threads from 16
        
        output:
        file "prokka/${genome_name}.${read_type}.${mode}.faa.gz"
        file "prokka/${genome_name}.${read_type}.${mode}.gff.gz"
        file "prokka/${genome_name}.${read_type}.${mode}.tsv.gz"
        
        """
    #!/bin/bash
    set -e;

    # Decompress the assembly
    gunzip -c ${contigs_fasta_gz} > ${contigs_fasta_gz.simpleName}.fasta

    prokka \
        --outdir prokka/ \
        --prefix ${genome_name}.${read_type}.${mode} \
        --cpus ${threads} \
        --compliant \
        ${contigs_fasta_gz.simpleName}.fasta


    gzip prokka/${genome_name}.${read_type}.${mode}.faa &&
    gzip prokka/${genome_name}.${read_type}.${mode}.gff &&
    gzip prokka/${genome_name}.${read_type}.${mode}.tsv
        """
    }
}

// Summarize the assembly
process summarizeAll {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    cpus 1
    memory "4 GB"
    errorStrategy "retry"
    publishDir "${params.output_folder}/"

    input:
    file all_summary_jsons from assembly_summaries_ch.collect()

    output:
    file "assembly_summary_table.csv"

    afterScript "rm -rf *"

"""
#!/usr/bin/env python3
import gzip
import json
import pandas as pd
import os

all_summary_jsons = "${all_summary_jsons}".split(" ")

# Make sure that the right columns exist
col_names = ["genome_name", "read_type", "mode", "num_contigs", "num_bases", "longest_contig"]

df = []
for fp in all_summary_jsons:
    assert os.path.exists(fp), "File not found in working directory (%s) -- try running again" % (fp)

    print("Reading in %s" % (fp))
    dat = json.load(open(fp, "rt"))
    assert isinstance(dat, dict), "%s is not the correct format" % (fp)

    for k in col_names:
        assert k in dat, "%s not found in summary (%s)" % (k, fp)

    df.append(dat)

print("Making a single summary table")
df = pd.DataFrame(df).sort_values(by=["genome_name", "read_type", "mode"])

df = df.reindex(columns=[
    "genome_name",
    "read_type",
    "mode",
    "num_contigs",
    "num_circular_contigs",
    "num_bases",
    "N50",
    "circular_contig_lengths",
    "linear_contig_lengths",
    "longest_contig",
    "num_over_1Mb",
    "num_100kb_to_1Mb",
    "num_10kb_to_100kb",
    "num_1kb_to_10kb",
    "num_under_1kb",
])

print("Writing out the summary table")
df.to_csv("assembly_summary_table.csv", index=None)

"""
}