#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// ===========================================================================
//  ChIP-seq pipeline
//    STAGE 1 : mapping (bwa-mem2 | bowtie2)   <-- this build
//    STAGE 3 : signal (bamCoverage / bamCompare) + ChIPseqSpikeInFree  [next]
//
//  QC lives in a SEPARATE pipeline (by request). No fastqc / fastp / deepTools
//  here: the aligners consume the RAW fastqs directly.
//
//  Both aligners emit an IDENTICAL BAM contract:
//    coordinate-sorted, mate-fixed, duplicate-marked-and-removed,
//    primary-only, indexed.  MAPQ is NOT baked in (applied per-mark in stage 3).
// ===========================================================================

process BWA_MEM2_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/mapping/bam",  mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/mapping/stats", mode: 'copy', pattern: '*.{tsv,txt}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: bam
    path("${meta.id}.flagstat.tsv"),      emit: flagstat
    path("${meta.id}.markdup.stats.txt"), emit: markdup

    script:
    def idx = params.bwamem2_index_prefix?.toString()?.trim()
    def rg  = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.library}\\tPU:${meta.unit}\\tPL:ILLUMINA"
    """
    ${params.bwamem2_bin} mem -t ${task.cpus} \\
        -R "${rg}" ${params.bwamem2_extra} \\
        ${idx} ${r1} ${r2} \\
      | samtools collate -@ ${task.cpus} -O -u - \\
      | samtools fixmate  -@ ${task.cpus} -m -u - - \\
      | samtools sort     -@ ${task.cpus} -u - \\
      | samtools markdup  -@ ${task.cpus} -f ${meta.id}.markdup.stats.txt - - \\
      | samtools view     -@ ${task.cpus} -b -q ${params.bam_min_mapq} -F ${params.exclude_flags} - \\
      | samtools sort     -@ ${task.cpus} -o ${meta.id}.filtered.bam -

    samtools index    -@ ${task.cpus} ${meta.id}.filtered.bam
    samtools flagstat -@ ${task.cpus} -O tsv ${meta.id}.filtered.bam > ${meta.id}.flagstat.tsv
    """

    stub:
    """
    touch ${meta.id}.filtered.bam ${meta.id}.filtered.bam.bai
    touch ${meta.id}.flagstat.tsv ${meta.id}.markdup.stats.txt
    """
}

process BOWTIE2_ALIGN {
    tag "${meta.id}"
    publishDir "${params.outdir}/mapping/bam",  mode: 'copy', pattern: '*.bam*'
    publishDir "${params.outdir}/mapping/stats", mode: 'copy', pattern: '*.{tsv,txt,log}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), emit: bam
    path("${meta.id}.flagstat.tsv"),      emit: flagstat
    path("${meta.id}.markdup.stats.txt"), emit: markdup
    path("${meta.id}.bowtie2.log"),       emit: log

    script:
    def idx = params.bowtie2_index_prefix?.toString()?.trim()
    """
    bowtie2 -x ${idx} -1 ${r1} -2 ${r2} \\
        --threads ${task.cpus} \\
        --rg-id ${meta.id} \\
        --rg SM:${meta.id} --rg LB:${meta.library} --rg PU:${meta.unit} --rg PL:ILLUMINA \\
        ${params.bowtie2_extra} \\
        2> ${meta.id}.bowtie2.log \\
      | samtools collate -@ ${task.cpus} -O -u - \\
      | samtools fixmate  -@ ${task.cpus} -m -u - - \\
      | samtools sort     -@ ${task.cpus} -u - \\
      | samtools markdup  -@ ${task.cpus} -f ${meta.id}.markdup.stats.txt - - \\
      | samtools view     -@ ${task.cpus} -b -q ${params.bam_min_mapq} -F ${params.exclude_flags} - \\
      | samtools sort     -@ ${task.cpus} -o ${meta.id}.filtered.bam -

    samtools index    -@ ${task.cpus} ${meta.id}.filtered.bam
    samtools flagstat -@ ${task.cpus} -O tsv ${meta.id}.filtered.bam > ${meta.id}.flagstat.tsv
    """

    stub:
    """
    touch ${meta.id}.filtered.bam ${meta.id}.filtered.bam.bai
    touch ${meta.id}.flagstat.tsv ${meta.id}.markdup.stats.txt ${meta.id}.bowtie2.log
    """
}

// ===========================================================================
//  STAGE 3 : signal (bamCoverage / bamCompare) + deepTools QC
//  BAMs are already dedup'd + flag-filtered, so deepTools does no extra
//  read filtering here (no --minMappingQuality downstream).
// ===========================================================================

// Per-file coverage bigWig
process BAM_COVERAGE {
    tag "${meta.id}"
    publishDir "${params.outdir}/BigWigs/bamCoverage", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.coverage.bw"), emit: bigwig

    script:
    def norm = params.bigwig_normalization?.toString()?.trim() ?: 'None'
    def egs  = (norm == 'RPGC') ? "--effectiveGenomeSize ${params.effective_genome_size}" : ''
    def bl   = params.blacklist_bed?.toString()?.trim() ? "--blackListFileName ${params.blacklist_bed}" : ''
    """
    bamCoverage \\
      --bam ${bam} \\
      --outFileName ${meta.id}.coverage.bw \\
      --binSize ${params.bigwig_binsize} \\
      --normalizeUsing ${norm} \\
      ${egs} ${bl} \\
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}.coverage.bw
    """
}

// ChIP vs its matched input: ratio track, no between-sample scale factor
process BAM_COMPARE {
    tag "${meta.id}"
    publishDir "${params.outdir}/BigWigs/bamCompare_Bin_${params.bamcompare_binsize}_Smooth_${params.bamcompare_smooth_length}", mode: 'copy'

    input:
    tuple val(controlId), val(meta), path(chipBam), path(chipBai), path(inputBam), path(inputBai)

    output:
    tuple val(meta), path("${meta.id}.ratio.bw"), emit: bigwig

    script:
    def bl = params.blacklist_bed?.toString()?.trim() ? "--blackListFileName ${params.blacklist_bed}" : ''
    """
    bamCompare \\
      -b1 ${chipBam} \\
      -b2 ${inputBam} \\
      --operation ratio \\
      --scaleFactorsMethod None \\
      --binSize ${params.bamcompare_binsize} \\
      --smoothLength ${params.bamcompare_smooth_length} \\
      ${bl} \\
      --outFileName ${meta.id}.ratio.bw \\
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}.ratio.bw
    """
}

// multiBigwigSummary #1 -- feeds PCA + correlation
process MULTIBIGWIG_SUMMARY_QC {
    publishDir "${params.outdir}/QC/Matrixes", mode: 'copy'

    input:
    tuple val(labels), path(bigwigs)

    output:
    path("multiBigwigSummary_QC.npz"), emit: npz

    script:
    def lab = labels.join(' ')
    """
    multiBigwigSummary bins \\
      --bwfiles ${bigwigs} \\
      --labels ${lab} \\
      --binSize ${params.mbws_qc_binsize} \\
      --outFileName multiBigwigSummary_QC.npz \\
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch multiBigwigSummary_QC.npz
    """
}

// multiBigwigSummary #2 -- 1000 bp bins + raw counts, for your later use
process MULTIBIGWIG_SUMMARY_COUNTS {
    publishDir "${params.outdir}/QC/Matrixes", mode: 'copy'

    input:
    tuple val(labels), path(bigwigs)

    output:
    path("multiBigwigSummary_bin${params.mbws_counts_binsize}.npz"),      emit: npz
    path("multiBigwigSummary_bin${params.mbws_counts_binsize}.rawcounts.tsv"), emit: counts

    script:
    def lab = labels.join(' ')
    """
    multiBigwigSummary bins \\
      --bwfiles ${bigwigs} \\
      --labels ${lab} \\
      --binSize ${params.mbws_counts_binsize} \\
      --outFileName multiBigwigSummary_bin${params.mbws_counts_binsize}.npz \\
      --outRawCounts multiBigwigSummary_bin${params.mbws_counts_binsize}.rawcounts.tsv \\
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch multiBigwigSummary_bin${params.mbws_counts_binsize}.npz
    touch multiBigwigSummary_bin${params.mbws_counts_binsize}.rawcounts.tsv
    """
}

// --- RPKM variant of the second matrix -------------------------------------
// multiBigwigSummary has no RPKM option; RPKM comes from the bigWigs it reads.
// So: build RPKM coverage tracks, then summarize them at the counts bin size.
process BAM_COVERAGE_RPKM {
    tag "${meta.id}"
    publishDir "${params.outdir}/BigWigs/bamCoverage_RPKM", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.rpkm.bw"), emit: bigwig

    script:
    def bl = params.blacklist_bed?.toString()?.trim() ? "--blackListFileName ${params.blacklist_bed}" : ''
    """
    bamCoverage \\
      --bam ${bam} \\
      --outFileName ${meta.id}.rpkm.bw \\
      --binSize ${params.bigwig_binsize} \\
      --normalizeUsing RPKM \\
      ${bl} \\
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch ${meta.id}.rpkm.bw
    """
}

process MULTIBIGWIG_SUMMARY_RPKM {
    publishDir "${params.outdir}/QC/Matrixes", mode: 'copy'

    input:
    tuple val(labels), path(bigwigs)

    output:
    path("multiBigwigSummary_bin${params.mbws_counts_binsize}.RPKM.npz"),          emit: npz
    path("multiBigwigSummary_bin${params.mbws_counts_binsize}.RPKM.rawcounts.tsv"), emit: counts

    script:
    def lab = labels.join(' ')
    """
    multiBigwigSummary bins \\
      --bwfiles ${bigwigs} \\
      --labels ${lab} \\
      --binSize ${params.mbws_counts_binsize} \\
      --outFileName multiBigwigSummary_bin${params.mbws_counts_binsize}.RPKM.npz \\
      --outRawCounts multiBigwigSummary_bin${params.mbws_counts_binsize}.RPKM.rawcounts.tsv \\
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch multiBigwigSummary_bin${params.mbws_counts_binsize}.RPKM.npz
    touch multiBigwigSummary_bin${params.mbws_counts_binsize}.RPKM.rawcounts.tsv
    """
}

process PLOT_PCA {
    publishDir "${params.outdir}/QC/PCA", mode: 'copy'

    input:
    path npz

    output:
    path("pca.png"),      emit: plot
    path("pca_data.tsv"), emit: data

    script:
    """
    plotPCA \\
      --corData ${npz} \\
      --plotFile pca.png \\
      --outFileNameData pca_data.tsv \\
      --plotTitle "PCA of binned coverage"
    """

    stub:
    """
    touch pca.png pca_data.tsv
    """
}

process PLOT_CORRELATION {
    publishDir "${params.outdir}/QC/Correlation", mode: 'copy'

    input:
    path npz

    output:
    path("correlation_spearman_heatmap.png"), emit: heatmap
    path("correlation_spearman_matrix.tsv"),  emit: matrix

    script:
    """
    plotCorrelation \\
      --corData ${npz} \\
      --corMethod spearman \\
      --whatToPlot heatmap \\
      --skipZeros --plotNumbers \\
      --colorMap RdYlBu \\
      --plotTitle "Spearman correlation" \\
      --plotFile correlation_spearman_heatmap.png \\
      --outFileCorMatrix correlation_spearman_matrix.tsv
    """

    stub:
    """
    touch correlation_spearman_heatmap.png correlation_spearman_matrix.tsv
    """
}

process PLOT_FINGERPRINT {
    publishDir "${params.outdir}/QC/Fingerprint", mode: 'copy'

    input:
    tuple val(labels), path(bams), path(bais)

    output:
    path("fingerprint.png"),         emit: plot
    path("fingerprint_metrics.tsv"), emit: metrics

    script:
    def lab = labels.join(' ')
    def bl  = params.blacklist_bed?.toString()?.trim() ? "--blackListFileName ${params.blacklist_bed}" : ''
    """
    plotFingerprint \\
      --bamfiles ${bams} \\
      --labels ${lab} \\
      --skipZeros \\
      --plotFile fingerprint.png \\
      --outQualityMetrics fingerprint_metrics.tsv \\
      ${bl} \\
      --numberOfProcessors ${task.cpus}
    """

    stub:
    """
    touch fingerprint.png fingerprint_metrics.tsv
    """
}

// ===========================================================================
//  Validation
// ===========================================================================

def normMapper() {
    return params.mapper?.toString()?.toLowerCase()?.trim() ?: 'bwa-mem2'
}

def isInput(meta) {
    return meta.mark.toString().equalsIgnoreCase('INPUT')
}

// Collect a (meta, bigwig) channel into ONE ordered (labels, bigwigs) tuple.
// Sorting by id keeps labels aligned with files.
def collectBigwigs(ch) {
    return ch
        .map { meta, bw -> tuple(meta.id, bw) }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows -> tuple(rows.collect { it[0] }, rows.collect { it[1] }) }
}

// Collect a (meta, bam, bai) channel into ONE ordered (labels, bams, bais) tuple.
def collectBams(ch) {
    return ch
        .map { meta, bam, bai -> tuple(meta.id, bam, bai) }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows -> tuple(rows.collect { it[0] }, rows.collect { it[1] }, rows.collect { it[2] }) }
}

def validateReference() {
    def mapper = normMapper()

    if (mapper == 'bwa-mem2') {
        def ref = params.bwamem2_index_prefix?.toString()?.trim()
        if (!ref)
            error "For bwa-mem2, set params.bwamem2_index_prefix to the FASTA basename you indexed (e.g. /path/to/hg38.fa)"
        def exts    = ['.0123', '.amb', '.ann', '.bwt.2bit.64', '.pac']
        def missing = exts.findAll { ext -> !file("${ref}${ext}").exists() }
        if (missing)
            error "bwa-mem2 index missing for '${ref}': ${missing.join(', ')}. Build once: bwa-mem2 index ${ref}"
    }
    else if (mapper == 'bowtie2') {
        def ref = params.bowtie2_index_prefix?.toString()?.trim()
        if (!ref)
            error "For bowtie2, set params.bowtie2_index_prefix (e.g. /path/to/hg38)"
        def exts    = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']
        def missing = exts.findAll { ext -> !file("${ref}${ext}").exists() }
        if (missing)
            error "bowtie2 index missing for '${ref}': ${missing.join(', ')}. Build once: bowtie2-build ref.fa ${ref}"
    }
    else {
        error "Invalid mapper '${mapper}'. Choose 'bwa-mem2' or 'bowtie2'."
    }
}

// ===========================================================================
//  Workflow
// ===========================================================================

workflow {
    // -stub-run is a wiring test; don't require real index files on disk for it.
    if( !workflow.stubRun )
        validateReference()

    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sampleId = row.sample?.trim()
            if (!sampleId)
                error "Every row in ${params.samples} needs a non-empty 'sample'"
            if (!row.fastq_r1 || !row.fastq_r2)
                error "Sample ${sampleId} is missing fastq_r1 or fastq_r2"

            def mark = row.mark?.trim()
            if (!mark)
                error "Sample ${sampleId} is missing 'mark' (e.g. H3K36me3, H3K9me3, or INPUT)"

            def genotype      = row.genotype?.trim() ?: 'unspecified'
            def controlSample = row.control_sample?.trim() ?: ''

            if (mark.equalsIgnoreCase('INPUT')) {
                controlSample = ''
            } else if (!controlSample) {
                error "Sample ${sampleId} is missing control_sample (the INPUT sample ID to pair against)"
            }

            def meta = [
                id            : sampleId,
                mark          : mark,
                genotype      : genotype,
                condition     : "${genotype}_${mark}",
                replicate     : row.replicate ?: '1',
                library       : row.read_group_library ?: 'lib1',
                unit          : row.read_group_unit ?: 'unit1',
                control_sample: controlSample
            ]
            tuple(meta, file(row.fastq_r1), file(row.fastq_r2))
        }

    // Cross-check every control_sample names a real sample (matters for stage 3
    // bamCompare pairing; a typo here would silently yield zero pairs later).
    samples_ch
        .map { meta, r1, r2 -> meta }
        .toList()
        .subscribe { metas ->
            def ids = metas.collect { it.id } as Set
            metas.findAll { it.control_sample }.each { m ->
                if (!ids.contains(m.control_sample))
                    error "Sample '${m.id}' lists control_sample='${m.control_sample}', " +
                          "which is not a 'sample' in ${params.samples}. Known: ${ids.sort().join(', ')}"
            }
        }

    // ---- mapping: RAW fastqs straight into the chosen aligner ----
    def mapper = normMapper()
    if (mapper == 'bwa-mem2') {
        BWA_MEM2_ALIGN(samples_ch)
        mapped_bam = BWA_MEM2_ALIGN.out.bam
    }
    else {
        BOWTIE2_ALIGN(samples_ch)
        mapped_bam = BOWTIE2_ALIGN.out.bam
    }

    // mapped_bam shape: tuple(meta, filtered.bam, filtered.bam.bai)

    // ---- BigWigs ----
    // per-file coverage track
    BAM_COVERAGE(mapped_bam)

    // ChIP vs matched input -> ratio track (pair on control_sample == input id)
    inputs_keyed_ch = mapped_bam
        .filter { meta, bam, bai -> isInput(meta) }
        .map    { meta, bam, bai -> tuple(meta.id, bam, bai) }

    chips_keyed_ch = mapped_bam
        .filter { meta, bam, bai -> !isInput(meta) }
        .map    { meta, bam, bai -> tuple(meta.control_sample, meta, bam, bai) }

    chip_input_pairs_ch = chips_keyed_ch.combine(inputs_keyed_ch, by: 0)

    BAM_COMPARE(chip_input_pairs_ch)

    // ---- QC ----
    // multiBigwigSummary #1 (QC bins) off coverage bigWigs -> PCA + correlation
    cov_collected_ch = collectBigwigs(BAM_COVERAGE.out.bigwig)

    MULTIBIGWIG_SUMMARY_QC(cov_collected_ch)
    PLOT_PCA(MULTIBIGWIG_SUMMARY_QC.out.npz)
    PLOT_CORRELATION(MULTIBIGWIG_SUMMARY_QC.out.npz)

    // multiBigwigSummary #2 (1 kb + raw counts), in TWO flavors for later use:
    //   raw   -> coverage (normalizeUsing None)
    //   RPKM  -> coverage (normalizeUsing RPKM)
    MULTIBIGWIG_SUMMARY_COUNTS(cov_collected_ch)

    BAM_COVERAGE_RPKM(mapped_bam)
    MULTIBIGWIG_SUMMARY_RPKM( collectBigwigs(BAM_COVERAGE_RPKM.out.bigwig) )

    // fingerprint runs on the BAMs (all samples, ChIP + input)
    PLOT_FINGERPRINT( collectBams(mapped_bam) )
}

