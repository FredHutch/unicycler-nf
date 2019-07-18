# Unicycler - Nextflow

_De novo_ assembly using UniCycler, within a Nextflow workflow

```
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
```
