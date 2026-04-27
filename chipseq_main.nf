#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC_RAW {
    tag "${meta.id}"
    publishDir "${params.outdir}/qc/raw/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.raw_fastqc.done")
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${reads}
    touch ${meta.id}.raw_fastqc.done
    """

    stub:
    """
    touch ${meta.id}.raw_fastqc.done
    touch ${meta.id}_R1_fastqc.html ${meta.id}_R1_fastqc.zip
    touch ${meta.id}_R2_fastqc.html ${meta.id}_R2_fastqc.zip
    """
}

process FASTP_TRIM {
    tag "${meta.id}"
    publishDir "${params.outdir}/fastq_trimmed", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.trimmed.R1.fastq.gz"), path("${meta.id}.trimmed.R2.fastq.gz"), emit: trimmed
    path("${meta.id}.fastp.html"), emit: html
    path("${meta.id}.fastp.json"), emit: json

    script:
    """
    fastp \
      --in1 ${r1} \
      --in2 ${r2} \
      --out1 ${meta.id}.trimmed.R1.fastq.gz \
      --out2 ${meta.id}.trimmed.R2.fastq.gz \
      --html ${meta.id}.fastp.html \
      --json ${meta.id}.fastp.json \
      --thread ${task.cpus} \
      ${params.fastp_extra}
    """

    stub:
    """
    touch ${meta.id}.trimmed.R1.fastq.gz
    touch ${meta.id}.trimmed.R2.fastq.gz
    touch ${meta.id}.fastp.html
    touch ${meta.id}.fastp.json
    """
}

process FASTQC_TRIMMED {
    tag "${meta.id}"
    publishDir "${params.outdir}/qc/trimmed/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.trimmed_fastqc.done")
    path("*_fastqc.html"), emit: html
    path("*_fastqc.zip"), emit: zip

    script:
    """
    fastqc -t ${task.cpus} ${r1} ${r2}
    touch ${meta.id}.trimmed_fastqc.done
    """

    stub:
    """
    touch ${meta.id}.trimmed_fastqc.done
    touch ${meta.id}.trimmed.R1_fastqc.html ${meta.id}.trimmed.R1_fastqc.zip
    touch ${meta.id}.trimmed.R2_fastqc.html ${meta.id}.trimmed.R2_fastqc.zip
    """
}

process BWA_MEM_FILTER {
    tag "${meta.id}"
    publishDir "${params.outdir}/mapping/bam", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), path("${meta.id}.filtered.bam.bai"), path("${meta.id}.flagstat.tsv"), emit: bam

    script:
    def ref = params.bwa_index_prefix ?: params.reference_fasta
    """
    bwa mem -t ${task.cpus} ${ref} ${r1} ${r2} \
      | samtools collate -@ ${task.cpus} -O -u - \
      | samtools fixmate -@ ${task.cpus} -m -u - - \
      | samtools sort -@ ${task.cpus} -u - \
      | samtools markdup -@ ${task.cpus} - - \
      | samtools view -@ ${task.cpus} -b -q ${params.min_mapq} -F ${params.exclude_flags} - \
      | samtools sort -@ ${task.cpus} -o ${meta.id}.filtered.bam -

    samtools index -@ ${task.cpus} ${meta.id}.filtered.bam
    samtools flagstat -@ ${task.cpus} -O tsv ${meta.id}.filtered.bam > ${meta.id}.flagstat.tsv
    """

    stub:
    """
    touch ${meta.id}.filtered.bam
    touch ${meta.id}.filtered.bam.bai
    touch ${meta.id}.flagstat.tsv
    """
}

process BAM_COVERAGE {
    tag "${meta.id}"
    publishDir "${params.outdir}/coverage/bigwig", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(flagstat)

    output:
    tuple val(meta), path("${meta.id}.bw"), emit: bw

    script:
    """
    bamCoverage \
      --bam ${bam} \
      --outFileName ${meta.id}.bw \
      --outFileFormat bigwig \
      --binSize ${params.bigwig_binsize} \
      --normalizeUsing ${params.bigwig_normalization} \
      --effectiveGenomeSize ${params.effective_genome_size} \
      --numberOfProcessors ${task.cpus} 
    """

    stub:
    """
    touch ${meta.id}.bw
    """
}


process BAM_COMPARE {
    tag "${meta.id}_vs_${controlId}"
    publishDir "${params.outdir}/coverage/bamcompare", mode: 'copy'

    input:
    tuple val(controlId), val(meta), path(chipBam), path(chipBai), path(inputBam), path(inputBai)

    output:
    tuple val(meta), path("${meta.id}.vs.${controlId}.ratio.bw"), emit: bw

    script:
    """
    bamCompare \
      -b1 ${chipBam} \
      -b2 ${inputBam} \
      -o ${meta.id}.vs.${controlId}.ratio.bw \
      --ignoreDuplicates \
      --numberOfProcessors ${task.cpus} \
      --operation ratio \
      --smoothLength ${params.bamcompare_smooth_length} \
      --binSize ${params.bamcompare_binsize}
    """

    stub:
    """
    touch ${meta.id}.vs.${controlId}.ratio.bw
    """
}

process MULTIQC {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path qc_inputs

    output:
    path 'multiqc_report.html', emit: report

    script:
    """
    multiqc . --filename multiqc_report.html
    """

    stub:
    """
    touch multiqc_report.html
    """
}


def validateBwaReference() {
    def ref = (params.bwa_index_prefix ?: params.reference_fasta)?.toString()?.trim()
    if (!ref) {
        error "Set either params.bwa_index_prefix or params.reference_fasta in chipseq_nextflow.config"
    }

    def expectedExts = ['.amb', '.ann', '.bwt', '.pac', '.sa']
    def missing = expectedExts.findAll { ext -> !file("${ref}${ext}").exists() }
    if (missing) {
        error "BWA index files are missing for '${ref}'. Missing: ${missing.join(', ')}. Set params.bwa_index_prefix to the index basename (e.g. /path/to/hg19) or params.reference_fasta to a FASTA with generated BWA index files."
    }
}

workflow {
    validateBwaReference()

    samples_ch = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sampleId = row.sample?.trim()
            if (!sampleId) {
                error "Each row in ${params.samples} must include a non-empty 'sample' value"
            }
            if (!row.fastq_r1 || !row.fastq_r2) {
                error "Sample ${sampleId} is missing fastq_r1 or fastq_r2"
            }

            def condition = row.condition ?: 'unspecified'
            def controlSample = row.control_sample?.trim()
            if (condition.toString().equalsIgnoreCase('Input')) {
                controlSample = ''
            } else if (!controlSample) {
                error "Sample ${sampleId} is missing control_sample. Provide the Input sample ID to pair for bamCompare."
            }

            def meta = [
                id: sampleId,
                condition: condition,
                replicate: row.replicate ?: '1',
                library: row.read_group_library ?: 'lib1',
                unit: row.read_group_unit ?: 'unit1',
                control_sample: controlSample
            ]
            tuple(meta, file(row.fastq_r1), file(row.fastq_r2))
        }

    raw_reads_ch = samples_ch.map { meta, r1, r2 -> tuple(meta, [r1, r2]) }

    FASTQC_RAW(raw_reads_ch)
    FASTP_TRIM(samples_ch)
    FASTQC_TRIMMED(FASTP_TRIM.out.trimmed)
    BWA_MEM_FILTER(FASTP_TRIM.out.trimmed)

    if (params.make_bigwig) {
        BAM_COVERAGE(BWA_MEM_FILTER.out.bam)
    }

    if (params.make_bamcompare) {
        input_bams_ch = BWA_MEM_FILTER.out.bam
            .filter { meta, bam, bai, flagstat -> meta.condition.toString().equalsIgnoreCase('Input') }
            .map { meta, bam, bai, flagstat -> tuple(meta.id, bam, bai) }

        chip_bams_ch = BWA_MEM_FILTER.out.bam
            .filter { meta, bam, bai, flagstat -> !meta.condition.toString().equalsIgnoreCase('Input') }
            .map { meta, bam, bai, flagstat -> tuple(meta.control_sample, meta, bam, bai) }

        chip_input_pairs_ch = chip_bams_ch.combine(input_bams_ch, by: 0)
            .map { controlId, meta, chipBam, chipBai, inputBam, inputBai ->
                tuple(controlId, meta, chipBam, chipBai, inputBam, inputBai)
            }

        BAM_COMPARE(chip_input_pairs_ch)
    }

    multiqc_inputs = FASTQC_RAW.out.html
        .mix(FASTQC_RAW.out.zip)
        .mix(FASTQC_TRIMMED.out.html)
        .mix(FASTQC_TRIMMED.out.zip)
        .mix(FASTP_TRIM.out.html)
        .mix(FASTP_TRIM.out.json)
        .mix(BWA_MEM_FILTER.out.bam.map { meta, bam, bai, flagstat -> flagstat })
        .collect()

    MULTIQC(multiqc_inputs)
}
