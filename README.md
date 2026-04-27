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

## Dashboard (for non-terminal users)

A Streamlit dashboard is available at `dashboard/app.py` to help users who are not comfortable with HPC terminals.

### What it supports

- Edit sample metadata and file paths directly in a table (sample name, FASTQ pairs, control pairing).
- Enter HPC login target (username, host, remote working directory) so users can run the pipeline with generated SSH command syntax.
- Set core run parameters (reference paths, MAPQ/filter settings, output directory, profile).
- Save `config/chipseq_samples.tsv` from the UI.
- Generate a dashboard override config at `config/chipseq_dashboard_override.config`.
- Launch the Nextflow command from the dashboard (optional).
- Optionally run through SSH (`user@host`) so execution happens on HPC rather than local machine.
- View a separate **QC & Processed Files** tab with counts/status for FastQC, fastp, BAM, bigWig, bamCompare outputs, and MultiQC presence.

### Run the dashboard

```bash
pip install streamlit pandas
streamlit run dashboard/app.py
```

The dashboard is designed so that users review QC and file status there, and download `.bw` files later using their normal HPC connection method (SCP/SFTP/Globus).
