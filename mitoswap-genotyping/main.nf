/*  GATK4 Variant Calling Pipeline 
 *  Usage: nextflow run /path/to/main.nf
 *
 *  Author: Mohammed Khalfan < mkhalfan@nyu.edu >
 *  NYU Center for Genetics and System Biology 2020
 */

// Setting some defaults here,
// can be overridden in config or via command line
params.out = "${params.outdir}/out"
params.tmpdir = "${params.outdir}/gatk_temp"
params.snpeff_data = "${params.outdir}/snpeff_data"

// Define modules here
BWA = 'bwa/intel/0.7.17'
PICARD = 'picard/2.17.11'
GATK = 'gatk/4.1.9.0'
R = 'r/intel/4.0.3'
SAMTOOLS = 'samtools/intel/1.11'
SNPEFF = 'snpeff/4.3t'
DEEPTOOLS = 'deeptools/3.5.0'
PYPAIRIX = 'pypairix/0.3.7'
HTSLIB = 'htslib/intel/1.11.0'
JVARKIT = 'jvarkit/base'
QUALIMAP = 'qualimap/2.2.1'
BCFTOOLS = 'bcftools/intel/1.11'
MULTIQC = 'multiqc/1.9'
TRIMMOMATIC = 'trimmomatic/0.36'

// Print some stuff here 
println "reads: $params.reads"
println "ref: $params.ref"
println "output: $params.out"
println "gatk temp dir: $params.tmpdir"
println "snpeff db: $params.snpeff_db"
println "snpeff data: $params.snpeff_data"

// Setup the reference file
ref = file(params.ref)

/* Prepare the fastq read pairs for input.
 * Use the size parameter to not auto-group, and instead 
 * use getBaseName() and remove two regexs to get the ID. 
 * This is custom for NYU CGSB sequence data file naming format 
 * While doing this, count number of input samples
 */
num_samples = 0
Channel
    .fromFilePairs( params.reads, size: -1)
    { file -> file.getBaseName() - ~/${params.fcid}_/ - ~/n0[12]_/ - ~/.fastq/ }
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}"  }
    .tap { read_pairs_ch }
    .subscribe({ num_samples += 1 })

process trim {
    publishDir "${params.out}/trimmed", mode:'copy'

    input:
    set pair_id,
        file(reads) from read_pairs_ch

    output:
    set val(pair_id),
	file("${pair_id}_trimmed_1.fq.gz"),
	file("${pair_id}_trimmed_2.fq.gz") \
	into trimmed_ch

    script:
    """
    module load $TRIMMOMATIC
    java -jar \$TRIMMOMATIC_JAR \
	PE \
	-phred33 \
	-threads ${task.cpus} \
	${reads[0]} \
	${reads[1]} \
	${pair_id}_trimmed_1.fq.gz \
	${pair_id}.unpair_trimmed_1.fq.gz \
	${pair_id}_trimmed_2.fq.gz \
	${pair_id}.unpair_trimmed_2.fq.gz \
	ILLUMINACLIP:${params.adapters}:2:30:10:8:true \
	LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:20
    """
}

process align {
    publishDir "${params.out}/aligned_reads", mode:'copy'
	
    input:
    set pair_id, 
	file(read_1),
	file(read_2) from trimmed_ch
     
    output:
    set val(pair_id), file("${pair_id}_aligned_reads.sam") \
	into aligned_reads_ch
	
    script:
    readGroup = \
	"@RG\\tID:${pair_id}\\tLB:${pair_id}\\tPL:${params.pl}\\tPM:${params.pm}\\tSM:${pair_id}"
    """
    module load $BWA
    bwa mem \
	-K 100000000 \
	-v 3 \
	-t ${task.cpus} \
	-Y \
	-R \"${readGroup}\" \
	$ref \
	$read_1 \
	$read_2 \
	> ${pair_id}_aligned_reads.sam
    """
}

process markDuplicatesSpark {
    publishDir "${params.out}/dedup_sorted", mode:'copy'

    input:
    set val(pair_id), file(aligned_reads) from aligned_reads_ch

    // If we're doing this step, it's the first round
    // so we set val(1) (round = 1).
    output:
    set val(pair_id), \
	val(1), \
	file("${pair_id}_sorted_dedup.bam") \
	into bam_for_variant_calling, \
	sorted_dedup_ch_for_metrics, \
	bam_for_bqsr
    set val(pair_id), \
	file ("${pair_id}_dedup_metrics.txt") \
	into dedup_qc_ch
    set val(pair_id),
        file("${pair_id}_sorted_dedup.bam"),
        file("${pair_id}_sorted_dedup.bam.bai") \
        into full_bam_bw_ch, downsample_bam_ch, qualimap_ch

    script:
    """
    module load $GATK
    gatk MarkDuplicatesSpark \
	-I $aligned_reads \
	-M ${pair_id}_dedup_metrics.txt \
	-O ${pair_id}_sorted_dedup.bam \
	--tmp-dir \${TMPDIR}
    """ 
}

process downsample_bam{
    publishDir "${params.out}/downsampled_bam", mode:'copy'

    input:
    set val(pair_id),
	file(bam),
	file(bam_index) \
	from downsample_bam_ch

    output:
    set file("${pair_id}_downsampled.cram"),
        file("${pair_id}_downsampled.cram.crai") into jbrowse_bam_ch

    when:
    params.do_jbrowse

    script:
    """
    module load $JVARKIT
    module load samtools/intel/1.14
    java -Xmx\${SLURM_MEM_PER_NODE}M -jar \${SORTSAMREFNAME_JAR} \
        --bamcompression 0 \
        --tmpDir \${TMPDIR} \
        --samoutputformat BAM \
        ${bam} | \
    java -Xmx\${SLURM_MEM_PER_NODE}M -jar \${BIOSTAR_JAR} \
        --bamcompression 0 \
        -n 75 \
        --samoutputformat BAM | \
    samtools sort \
        -l 0 \
        --threads \${SLURM_CPUS_PER_TASK} \
        -T \${TMPDIR} \
        --output-fmt BAM | \
    samtools view \
        --threads \${SLURM_CPUS_PER_TASK} \
        --reference $ref \
        -C \
        -o ${pair_id}_downsampled.cram \
        --write-index
    """

}

process qualimap{
    input:
    set val(pair_id),
        file(bam),
	file(bam_index) from qualimap_ch

    output:
    file('*') into multiqc_qualimap_ch

    script:
    """
    module load $QUALIMAP
    qualimap BamQC -bam $bam \
      -outdir ${pair_id} \
      -outformat HTML \
      -nt \${SLURM_CPUS_PER_TASK} \
      --java-mem-size=\${SLURM_MEM_PER_NODE}m 
    """
}

process getMetrics{
    publishDir "${params.out}/metrics", mode:'copy'

    input:
    set val(pair_id), \
	val(round), \
	file(sorted_dedup_reads) \
	from sorted_dedup_ch_for_metrics

    output:
    set val(pair_id), 
	file("${pair_id}_alignment_metrics.txt"), \
	file("${pair_id}_insert_metrics.txt"), \
	file("${pair_id}_insert_size_histogram.pdf"), \
	file("${pair_id}_depth_out.txt") \
	into metrics_qc_ch, metrics_multiqc_ch

    script:
    """
    module load $PICARD
    module load $R
    module load $SAMTOOLS
    java -jar \$PICARD_JAR \
        CollectAlignmentSummaryMetrics \
	R=${params.ref} \
        I=${sorted_dedup_reads} \
	O=${pair_id}_alignment_metrics.txt
    java -jar \$PICARD_JAR \
        CollectInsertSizeMetrics \
        INPUT=${sorted_dedup_reads} \
	OUTPUT=${pair_id}_insert_metrics.txt \
        HISTOGRAM_FILE=${pair_id}_insert_size_histogram.pdf 
    samtools depth -a ${sorted_dedup_reads} > ${pair_id}_depth_out.txt
    """
}

