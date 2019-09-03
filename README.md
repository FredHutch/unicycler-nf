# Unicycler - Nextflow

_De novo_ assembly using UniCycler, within a Nextflow workflow

```
Usage:

nextflow run fredhutch/unicycler-nf <ARGUMENTS>

Arguments:
  --sample_sheet       CSV file listing samples to analyze
  --output_folder      Folder to place outputs

Options:
  --short_reads        Sample sheet contains short read data (`short_R1` and `short_R2`)
  --long_reads         Sample sheet contains long read data (`long_reads`)
  --min_fasta_length   Minimum contig length (default: 100)
  --help               Display this message

Sample Sheet:
  The sample_sheet is a CSV with a header indicating which samples correspond to which files.
  The file must contain the column `name`, and `long_reads`, `short_R1`, `short_R2` as appropriate.
```
