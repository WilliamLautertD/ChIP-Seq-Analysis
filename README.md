# ChIP-Seq-Analysis
Generic ChIP-seq Nextflow workflow

## Generic ChIP-seq Nextflow workflow
A generic Nextflow workflow is available in `chipseq_main.nf` with configuration in `chipseq_nextflow.config`. It covers **mapping, signal-track generation, and deepTools QC**. Read-level QC (FastQC / fastp / MultiQC) is handled by a separate QC workflow `UPLOAD`.

### Scope
- **Mapping** with a choice of two aligners: `bwa-mem2` or `bowtie2` (selected by `params.mapper`). Reads are mapped directly from the raw FASTQs (no trimming step).
- **Post-alignment processing** (identical for both aligners): mate-fixing, coordinate sorting, duplicate marking and removal, SAM-flag filtering, optional MAPQ filtering, indexing. Produces `flagstat` and `markdup` stats.
- **Signal tracks** via deepTools: per-sample `bamCoverage` (raw and RPKM) and ChIP-vs-Input `bamCompare` (ratio).
- **QC** via deepTools: `plotPCA`, `plotCorrelation`, `plotFingerprint`, and `multiBigwigSummary` count matrices.

### Why this is generic
- Input FASTQ files are fully controlled by `config/chipseq_samples.tsv` (no fixed filename pattern is required).
- Sample names can be any value in the `sample` column.
- The `control_sample` column pairs each ChIP sample with the exact Input sample ID for `bamCompare`. A startup check fails fast if any `control_sample` does not name a real sample.
- Paths and run parameters are set in `chipseq_nextflow.config`: reference index, mapping/filter settings, bigWig options, `bamCompare` settings (ratio mode, no scale-factor normalization: `--scaleFactorsMethod None`), and `multiBigwigSummary` bin sizes.

### Mapping and reference index
Set `params.mapper` to `'bwa-mem2'` or `'bowtie2'` and provide the matching index basename. The two index types are **not** interchangeable:

- **bwa-mem2** — build with `bwa-mem2 index genome.fa`; set `bwamem2_index_prefix` to the FASTA basename **including** the extension (e.g. `/path/to/hg19.fa`). Expects `.0123 .amb .ann .bwt.2bit.64 .pac`.
- **bowtie2** — build with `bowtie2-build ref.fa hg19`; set `bowtie2_index_prefix` to the basename **without** `.fa` (e.g. `/path/to/hg19`). Expects `.1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2`.

The index is validated at launch (skipped automatically under `-stub-run`).

### Sample sheet
`config/chipseq_samples.tsv` is tab-separated with these columns:

| column | meaning |
|---|---|
| `sample` | unique sample ID (any string) |
| `fastq_r1`, `fastq_r2` | paths to the paired FASTQs |
| `mark` | histone mark / target, or `INPUT` for control samples |
| `genotype` | condition label (e.g. `EV`, `MARKER_OE`) |
| `replicate` | replicate number |
| `read_group_library`, `read_group_unit` | read-group `LB` / `PU` values |
| `control_sample` | for ChIP rows, the `sample` ID of the matched Input; leave blank for `INPUT` rows |

### Files
- `chipseq_main.nf` — mapping (bwa-mem2 or bowtie2), duplicate marking/removal and filtering, per-sample bigWigs (coverage + RPKM), ChIP/Input `bamCompare` (ratio), and deepTools QC (PCA, correlation, fingerprint, `multiBigwigSummary` matrices).
- `chipseq_nextflow.config` — pipeline parameters and runtime profiles (`conda`, `slurm`).
- `config/chipseq_samples.tsv` — sample sheet template; edit file paths, sample IDs, `mark`/`genotype`, and `control_sample` pairing.

### Quick start
The config file is **not** auto-loaded (it is not named `nextflow.config`), so pass it with `-c` on every run — or symlink it once (`ln -s chipseq_nextflow.config nextflow.config`).

Bowtie2 example:
```bash
nextflow run chipseq_main.nf -c chipseq_nextflow.config -profile conda -resume \
  --mapper bowtie2 \
  --bowtie2_index_prefix /path/to/hg19
```

bwa-mem2 example:
```bash
nextflow run chipseq_main.nf -c chipseq_nextflow.config -profile conda -resume \
  --mapper bwa-mem2 \
  --bwamem2_index_prefix /path/to/hg19.fa
```

For SLURM, swap the profile:
```bash
nextflow run chipseq_main.nf -c chipseq_nextflow.config -profile slurm -resume \
  --mapper bowtie2 --bowtie2_index_prefix /path/to/hg19
```

Wiring test (no index or aligner needed; runs process stubs only):
```bash
nextflow run chipseq_main.nf -c chipseq_nextflow.config -profile conda -stub-run
```

### Output structure
```
results/chipseq_pipeline/
  mapping/
    bam/                              *.filtered.bam (+ .bai)
    stats/                            *.flagstat.tsv, *.markdup.stats.txt, (bowtie2 *.log)
  BigWigs/
    bamCoverage/                      *.coverage.bw        (per sample, normalizeUsing None)
    bamCoverage_RPKM/                 *.rpkm.bw            (per sample, RPKM)
    bamCompare_Bin_200_Smooth_600/    *.ratio.bw          (ChIP / Input ratio)
  QC/
    PCA/            pca.png, pca_data.tsv
    Correlation/    correlation_spearman_heatmap.png, correlation_spearman_matrix.tsv
    Fingerprint/    fingerprint.png, fingerprint_metrics.tsv
    Matrixes/       multiBigwigSummary_QC.npz                        (10 kb -> PCA/correlation)
                    multiBigwigSummary_bin1000.npz + .rawcounts.tsv       (raw coverage)
                    multiBigwigSummary_bin1000.RPKM.npz + .rawcounts.tsv  (RPKM, for later use)
```

The `bamCompare` folder name is derived from `bamcompare_binsize` / `bamcompare_smooth_length`, so it renames itself if you change those values.

### Notes
- **MAPQ**: `bam_min_mapq` defaults to `0` (no MAPQ filter) so that repeat-region reads are retained. Set it to `30` in the config for the standard uniquely-mapped filter. It is applied once during mapping; deepTools does no further MAPQ filtering.
- **Genome**: defaults target **hg19** (`effective_genome_size = 2864785220`). `effective_genome_size` is only used when `bigwig_normalization = 'RPGC'`. Set `blacklist_bed` to an hg19 blacklist to exclude artifact regions from bigWigs, bamCompare, and fingerprint.
- **Conda env**: `envs/qc_mapping_cnv.yaml` must provide the chosen aligner (`bwa-mem2` and/or `bowtie2`), `samtools`, and `deeptools`.

## Dashboard (for non-terminal users) - not available 
A Streamlit dashboard is available at `dashboard/app.py` to help users who are not comfortable with HPC terminals.

### What it supports
- Edit sample metadata and file paths directly in a table (sample name, FASTQ pairs, `mark`/`genotype`, control pairing).
- Enter HPC login target (username, host, remote working directory) so users can run the pipeline with generated SSH command syntax.
- Set core run parameters (reference index prefix, mapper choice, MAPQ/filter settings, output directory, profile).
- Save `config/chipseq_samples.tsv` from the UI.
- Generate a dashboard override config at `config/chipseq_dashboard_override.config`.
- Launch the Nextflow command from the dashboard (optional).
- Optionally run through SSH (`user@host`) so execution happens on HPC rather than the local machine.
- View a separate **QC & Processed Files** tab with counts/status for BAMs, bigWigs (coverage + RPKM), `bamCompare` ratio tracks, and deepTools QC outputs (PCA, correlation, fingerprint, `multiBigwigSummary` matrices).

> Note: this pipeline does not produce FastQC/fastp/MultiQC outputs; read-level QC lives in the separate QC workflow. If `dashboard/app.py` still references those, update it to point at the deepTools QC outputs above.

### Run the dashboard
```bash
pip install streamlit pandas
streamlit run dashboard/app.py
```
The dashboard is designed so that users review QC and file status there, and download `.bw` files later using their normal HPC connection method (SCP/SFTP/Globus).
