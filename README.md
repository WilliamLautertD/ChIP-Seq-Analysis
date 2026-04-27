# ChIP-Seq-Analysis
Generic ChIP-seq Nextflow workflow

## Generic ChIP-seq Nextflow workflow

For ChIP-seq, a separate generic Nextflow workflow is available in `chipseq_main.nf` with configuration in `chipseq_nextflow.config`.

### Why this is generic

- Input FASTQ files are fully controlled by `config/chipseq_samples.tsv` (no fixed filename pattern is required).
- Sample names can be any value in the `sample` column.
- The `control_sample` column pairs each ChIP sample with the exact Input sample ID for downstream `bamCompare` output.
- Paths and run parameters are set in `chipseq_nextflow.config`, including reference genome, mapping filters, bigWig options, and `bamCompare` settings (ratio mode with duplicate ignoring, no scale-factor normalization).

### Files

- `chipseq_main.nf` - ChIP-seq QC, trimming, mapping/filtering, optional bigWig generation, optional ChIP/Input `bamCompare`, and MultiQC.
- `chipseq_nextflow.config` - ChIP-seq pipeline parameters and runtime profiles.
- `config/chipseq_samples.tsv` - sample sheet template; edit file paths, sample IDs, and `control_sample` pairing.

### Quick start

```bash
nextflow run chipseq_main.nf -c chipseq_nextflow.config -profile conda -resume
```

For SLURM:

```bash
nextflow run chipseq_main.nf -c chipseq_nextflow.config -profile slurm -resume
```
