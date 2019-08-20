#!/usr/bin/env nextflow
/*
========================================================================================
                         SteadyFlow - Steady Sate Transcription PIPELINE
========================================================================================
Steady Sate Analysis Pipeline. Started 2018-06-21.
 #### Homepage / Documentation
 https://github.com/Dowell-Lab/RNA-seq-Flow
 #### Authors
 Margaret Gruca <magr0763@colorado.edu>
========================================================================================
========================================================================================

Pipeline steps:

    1. Pre-processing sra/fastq
        1a. SRA tools -- fasterq-dump sra to generate fastq file
        1b. FastQC (pre-trim) -- perform pre-trim FastQC on fastq files

    2. Trimming
        2a. BBDuk -- trim fastq files for quality and adapters
        2b. FastQC (post-trim) -- perform post-trim FastQC on fastq files (ensure trimming performs as expected)

    3. Mapping w/ HISAT2 -- map to genome reference file

    4. SAMtools -- convert SAM file to BAM, index BAM, flagstat BAM

    5. Quality control
        5a. preseq -- estimate library complexity
        5b. RSeQC -- calculate genomic coverage relative to a reference file, infer experiement (single- v. paired-end), read duplication
        5c. Pileup.sh : BBMap Suite -- genomic coverage by chromosome, GC content, pos/neg reads, intron/exon ratio

    6. Coverage files
        6d. BEDTools : non-normalized & nornmalized bedgraphs
        6b. BEDTools and kentUtils : 5' bigwigs for dREG & normalized bigwigs for genome browser
        
    7. Normalizing bigwigs for Genome Browser use

    8. IGV Tools : bedGraph --> tdf
    
    9. RSeQC Counts : Read counts for a provided RefSeq annotation file

    10. MultiQC : generate QC report for pipeline

*/


def helpMessage() {
    log.info"""
    =========================================
     SteadyFlow v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf -profile slurm --fastqs '/project/*_{R1,R2}*.fastq' --outdir '/project/'

    Required arguments:
         -profile                      Configuration profile to use. <base, slurm>
         --fastqs                      Directory pattern for fastq files: /project/*{R1,R2}*.fastq (Required if --sras not specified)
         --sras                        Directory pattern for SRA files: /project/*.sras (Required if --fastqs not specified)
         --workdir                     Nextflow working directory where all intermediate files are saved.
         --email                       Where to send workflow report email.

    Performance options:
        --threadfqdump                 Runs multi-threading for fastq-dump for sra processing.

    Input File options:
        --singleEnd                    Specifies that the input files are not paired reads (default is paired-end).
        --flip                         Reverse complements each strand. Necessary for some library preps.
        --flipR2                       Reverse complements R2 only.       

    Save options:
        --outdir                       Specifies where to save the output from the nextflow run.
        --savefq                       Compresses and saves raw fastq reads.
        --saveTrim                     Compresses and saves trimmed fastq reads.
        --skipBAM                      Skip saving BAM files. Only CRAM files will be saved with this option.
        --saveAll                      Compresses and saves all fastq reads.
        --savebw                       Saves pos/neg bigwig files for UCSC genome browser.        
        --savebg                       Saves concatenated pos/neg bedGraph file.     
        
    QC Options:
        --skipMultiQC                  Skip running MultiQC.
        --skipRSeQC                    Skip running RSeQC.
        
    Analysis Options:
        --count                        Run RSeQC FPKM count over RefSeq annotated genes.

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

// Configurable variables
params.name = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false
params.bbmap_adapters = "$baseDir/bin/adapters.fa"
params.bedGraphToBigWig = "$baseDir/bin/bedGraphToBigWig"
params.rcc = "$baseDir/bin/rcc.py"
params.workdir = "./nextflowTemp"

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs

if ( params.genome ){
    genome = file(params.genome)
    if( !genome.exists() ) exit 1, "Genome directory not found: ${params.genome}"
}

if ( params.chrom_sizes ){
    chrom_sizes = file(params.chrom_sizes)
    if( !chrom_sizes.exists() ) exit 1, "Genome chrom sizes file not found: ${params.chrom_sizes}"
 }

if ( params.hisat2_indices ){
    hisat2_indices = file("${params.hisat2_indices}")
}

if ( params.genome_refseq ){
    genome_refseq = file("${params.genome_refseq}")
}

if ( params.bbmap_adapters){
    bbmap_adapters = file("${params.bbmap_adapters}")
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}



/*
 * Create a channel for input read files
 */
if (params.fastqs) {
    if (params.singleEnd) {
        fastq_reads_qc = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
        fastq_reads_trim = Channel
                            .fromPath(params.fastqs)
                            .map { file -> tuple(file.simpleName, file) }
    } else {
        Channel
            .fromFilePairs( params.fastqs, size: params.singleEnd ? 1 : 2 )
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
            .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
    }
}

else {
    Channel
        .empty()
        .into { fastq_reads_qc; fastq_reads_trim; fastq_reads_gzip }
}

if (params.sras) {
    println("pattern for SRAs provided")
    read_files_sra = Channel
                        .fromPath(params.sras)
                        .map { file -> tuple(file.baseName, file) }
}

else {
    read_files_sra = Channel.empty()
}


// Header log info
log.info """=======================================================
NascentFlow v${params.version}"
======================================================="""
def summary = [:]
summary['Pipeline Name']    = 'NascentFlow'
summary['Help Message']     = params.help
summary['Pipeline Version'] = params.version
summary['Run Name']         = custom_runName ?: workflow.runName
if(params.reads) summary['Reads']     = params.reads
if(params.fastqs) summary['Fastqs']   = params.fastqs
if(params.sras) summary['SRAs']       = params.sras
summary['Genome Ref']       = params.genome
summary['Thread fqdump']    = params.threadfqdump ? 'YES' : 'NO'
summary['Data Type']        = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Save All fastq']   = params.saveAllfq ? 'YES' : 'NO'
summary['Save BAM']         = params.skipBAM ? 'NO' : 'YES'
summary['Save BigWig']      = params.savebw ? 'YES' : 'NO'
summary['Save bedGraph']    = params.savebg ? 'YES' : 'NO'
summary['Save fastq']       = params.savefq ? 'YES' : 'NO'
summary['Save Trimmed']     = params.saveTrim ? 'YES' : 'NO'
summary['Reverse Comp']     = params.flip ? 'YES' : 'NO'
summary['Reverse Comp R2']  = params.flipR2 ? 'YES' : 'NO'
summary['Run RSeQC']        = params.skipRSeQC ? 'NO' : 'YES'
summary['Run Count']        = params.count ? 'NO' : 'YES'
summary['Run MultiQC']      = params.skipMultiQC ? 'NO' : 'YES'
summary['Max Memory']       = params.max_memory
summary['Max CPUs']         = params.max_cpus
summary['Max Time']         = params.max_time
summary['Output dir']       = params.outdir
summary['Working dir']      = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']     = "$HOME"
summary['Current user']     = "$USER"
summary['Current path']     = "$PWD"
summary['Output dir']       = params.outdir
summary['Script dir']       = workflow.projectDir
summary['Config Profile']   = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "======================================================="

// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    validExitStatus 0,1,127
    publishDir "${params.outdir}/software_versions/", mode: 'copy', pattern: '*.txt'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file '*.txt' into software_versions_text

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    bbversion.sh --version > v_bbduk.txt
    hisat2 --version > v_hisat2.txt
    samtools --version > v_samtools.txt
    fastq-dump --version > v_fastq-dump.txt
    preseq > v_preseq.txt
    bedtools --version > v_bedtools.txt
    igvtools version > v_igv-tools.txt

    for X in `ls *.txt`; do
        cat \$X >> all_versions.txt;
    done
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/*
 * Step 1a -- get fastq files from downloaded sras
 */

process sra_dump {
    tag "$prefix"
    if (params.threadfqdump) {
        cpus 8 }
    else {
        cpus 1
    }
    if (params.savefq || params.saveAllfq) {
        publishDir "${params.outdir}/fastq", mode: 'copy'
    }
    
    input:
    set val(prefix), file(reads) from read_files_sra

    output:
    set val(prefix), file("*.fastq.gz") into fastq_reads_qc_sra, fastq_reads_trim_sra, fastq_reads_gzip_sra
   

    script:
    prefix = reads.baseName
    if (!params.threadfqdump) {
        """
        echo ${prefix}

        fastq-dump ${reads} --gzip
        """
    } else if (!params.singleEnd) {
         """
        export PATH=~/.local/bin:$PATH

        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --split-3 \
            --sra-id ${reads}
        """
    } else if (!params.threadfqdump && !params.singleEnd) {
        """
        echo ${prefix}

        fastq-dump --split-3 ${reads} --gzip
        """
    } else {
        """
        export PATH=~/.local/bin:$PATH

        parallel-fastq-dump \
            --threads 8 \
            --gzip \
            --sra-id ${reads}
        """
    }
}

/*
 * STEP 1b - FastQC
 */

process fastQC {
    validExitStatus 0,1
    tag "$prefix"
    memory '8 GB'
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(prefix), file(reads) from fastq_reads_qc.mix(fastq_reads_qc_sra)

    output:
    file "*.{zip,html,txt}" into fastqc_results

    script:
    """
    fastqc $reads
    """
}


/*
 * STEP 2a - Trimming
 */

process bbduk {
    validExitStatus 0,1
    tag "$name"
    cpus 16
    time '2h'
    memory '20 GB'
    publishDir "${params.outdir}/qc/trimstats", mode: 'copy', pattern: "*.txt"
    if (params.saveTrim || params.saveAllfq) {
        publishDir "${params.outdir}/fastq_trimmed", mode: 'copy', pattern: "*.fastq.gz"
    }      

    input:
    set val(name), file(reads) from fastq_reads_trim.mix(fastq_reads_trim_sra)

    output:
    set val(name), file ("*.trim.fastq.gz") into trimmed_reads_fastqc, trimmed_reads_hisat2
    file "*.txt" into trim_stats

    script:
    if (!params.singleEnd && params.flip) {
        """
        echo ${name}

        reformat.sh -Xmx20g \
                t=16 \
                in=${name}_R1.fastq.gz \
                in2=${name}_R2.fastq.gz \
                out=${name}_R1.flip.fastq.gz \
                out2=${name}_R2.flip.fastq.gz \
                rcomp=t

        bbduk.sh -Xmx20g \
                t=16 \
                in=${name}_R1.flip.fastq.gz \
                in2=${name}_R2.flip.fastq.gz \
                out=${name}_R1.flip.trim.fastq.gz \
                out2=${name}_R2.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${name}.trimstats.txt \
                refstats=${name}.refstats.txt
        """
    } else if (!params.singleEnd && params.flipR2) {
                """
        echo ${name}

        reformat.sh -Xmx20g \
                t=16 \
                in=${name}_R1.fastq.gz \
                in2=${name}_R2.fastq.gz \
                out=${name}_R1.flip.fastq.gz \
                out2=${name}_R2.flip.fastq.gz \
                rcompmate=t

        bbduk.sh -Xmx20g \
                t=16 \
                in=${name}_R1.flip.fastq.gz \
                in2=${name}_R2.flip.fastq.gz \
                out=${name}_R1.flip.trim.fastq.gz \
                out2=${name}_R2.flip.trim.fastq.gz \
                ref=${bbmap_adapters} \
                ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                nullifybrokenquality=t \
                maq=10 minlen=25 \
                tpe tbo \
                literal=AAAAAAAAAAAAAAAAAAAAAAA \
                stats=${name}.trimstats.txt \
                refstats=${name}.refstats.txt
        """
    } else if (params.flip) {
        """
        echo ${name}


        reformat.sh -Xmx20g \
                t=16 \
                in=${name}.fastq.gz \
                out=${name}.flip.fastq.gz \
                rcomp=t

        
        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}.flip.fastq.gz \
                  out=${name}.flip.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  nullifybrokenquality=t \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    }
        else if (!params.singleEnd) {
        """
        echo ${name}      

        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}_R1.fastq.gz \
                  in2=${name}_R2.fastq.gz \
                  out=${name}_R1.trim.fastq.gz \
                  out2=${name}_R2.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  nullifybrokenquality=t \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    } else {
        """
        echo ${name}
        
        bbduk.sh -Xmx20g \
                  t=16 \
                  in=${name}.fastq.gz \
                  out=${name}.trim.fastq.gz \
                  ref=${bbmap_adapters} \
                  ktrim=r qtrim=10 k=23 mink=11 hdist=1 \
                  nullifybrokenquality=t \
                  maq=10 minlen=25 \
                  tpe tbo \
                  literal=AAAAAAAAAAAAAAAAAAAAAAA \
                  stats=${name}.trimstats.txt \
                  refstats=${name}.refstats.txt
        """
    }
}


/*
 * STEP 2b - Trimmed FastQC
 */

process fastqc_trimmed {
    validExitStatus 0,1
    tag "$name"
    memory '4 GB'
    publishDir "${params.outdir}/qc/fastqc/", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(trimmed_reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html,txt}" into trimmed_fastqc_results

    script:
    """
    echo ${name}

    fastqc $trimmed_reads
    """
}


/*
 * STEP 3 - Map reads to reference genome
 */

process hisat2 {
    tag "$name"
    validExitStatus 0
    cpus 32
    memory '100 GB'
    time '2h'
    publishDir "${params.outdir}/qc/hisat2_mapstats", mode: 'copy', pattern: "*.txt"

    input:
    file(indices) from hisat2_indices
    val(indices_path) from hisat2_indices
    set val(name), file(trimmed_reads) from trimmed_reads_hisat2

    output:
    set val(name), file("*.sam") into hisat2_sam
    file("*.txt") into hisat2_mapstats    

    script:
    if (!params.singleEnd) {
        """
        echo ${name}
    
        hisat2  -p 32 \
                --very-sensitive \
                -x ${indices_path} \
                --pen-noncansplice 14 \
                --mp 1,0 \
                --sp 3,1 \
                -1 ${name}_R1.trim.fastq.gz \
                -2 ${name}_R2.trim.fastq.gz \
                --new-summary \
                > ${name}.sam \
                2> ${name}.hisat2_mapstats.txt                
        """
    } else {
        """
        echo ${name}
    
        hisat2  -p 32 \
                --very-sensitive \
                --pen-noncansplice 14 \
                --mp 1,0 \
                --sp 3,1 \
                -x ${indices_path}\
                -U ${trimmed_reads} \
                --new-summary \
                > ${name}.sam \
                2> ${name}.hisat2_mapstats.txt                
        """
    }
}

/*
 * STEP 4 - Convert to BAM format and sort
 */

process samtools {
    tag "$name"
    memory '100 GB'
    cpus 16
    publishDir "${params.outdir}" , mode: 'copy',
    saveAs: {filename ->
             if ((filename.indexOf("sorted.bam") > 0) & !params.skipBAM)                                                                                                                             "mapped/bams/$filename"
        else if ((filename.indexOf("sorted.bam.bai") > 0) & !params.skipBAM)                                                                                                                         "mapped/bams/$filename"
        else if (filename.indexOf("flagstat") > 0)                    "qc/mapstats/$filename"
        else if (filename.indexOf("millionsmapped") > 0)              "qc/mapstats/$filename"
        else if (filename.indexOf("sorted.cram") > 0)                 "mapped/crams/$filename"
        else if (filename.indexOf("sorted.cram.crai") > 0)            "mapped/crams/$filename"
    }

    input:
    set val(name), file(mapped_sam) from hisat2_sam

    output:
    set val(name), file("${name}.sorted.bam") into sorted_bam_ch
    set val(name), file("${name}.sorted.bam.bai") into sorted_bam_indices_ch
    set val(name), file("${name}.flagstat") into bam_flagstat
    set val(name), file("${name}.millionsmapped") into bam_milmapped_bedgraph
    set val(name), file("${name}.sorted.cram") into cram_out
    set val(name), file("${name}.sorted.cram.crai") into cram_index_out

    script:
    if (!params.singleEnd) {
    """

    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x40 ${name}.sorted.bam | cut -f1 | sort | uniq | wc -l > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
    } else {
    """

    samtools view -@ 16 -bS -o ${name}.bam ${mapped_sam}
    samtools sort -@ 16 ${name}.bam > ${name}.sorted.bam
    samtools flagstat ${name}.sorted.bam > ${name}.flagstat
    samtools view -@ 16 -F 0x904 -c ${name}.sorted.bam > ${name}.millionsmapped
    samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
    samtools view -@ 16 -C -T ${genome} -o ${name}.cram ${name}.sorted.bam
    samtools sort -@ 16 -O cram ${name}.cram > ${name}.sorted.cram
    samtools index -c ${name}.sorted.cram ${name}.sorted.cram.crai
    """
    }
}

sorted_bam_ch
   .into {sorted_bams_for_bedtools_bedgraph; sorted_bams_for_preseq; sorted_bams_for_rseqc; sorted_bams_for_dreg_prep; sorted_bams_for_pileup; sorted_bams_for_rseqc_count}

sorted_bam_indices_ch
    .into {sorted_bam_indices_for_bedtools_bedgraph; sorted_bam_indices_for_bedtools_normalized_bedgraph; sorted_bam_indicies_for_pileup; sorted_bam_indices_for_preseq; sorted_bam_indices_for_rseqc; sorted_bam_indices_for_rseqc_count}

/*
 *STEP 5a - Plot the estimated complexity of a sample, and estimate future yields
 *         for complexity if the sample is sequenced at higher read depths.
 */

process preseq {
    tag "$name"
    memory '20 GB'
    time '8h'
    errorStrategy 'ignore'
    publishDir "${params.outdir}/qc/preseq/", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_preseq
    file(bam_indices) from sorted_bam_indices_for_preseq

    output:
    file("*.txt") into preseq_results

    script:
    """
    preseq c_curve -B -o ${name}.c_curve.txt \
           ${bam_file}

    preseq lc_extrap -B -o ${name}.lc_extrap.txt \
           ${bam_file}
    """
 }


/*
 *STEP 5b - Analyze read distributions using RSeQC
 */

process rseqc_qc {
    tag "$name"
    time '8h'
    validExitStatus 0,143
    memory '40 GB'
    publishDir "${params.outdir}/qc/rseqc" , mode: 'copy',
        saveAs: {filename ->
             if (filename.indexOf("infer_experiment.txt") > 0)                       "infer_experiment/$filename"
        else if (filename.indexOf("read_distribution.txt") > 0)                      "read_distribution/$filename"
        else if (filename.indexOf("read_duplication.DupRate_plot.pdf") > 0)          "read_duplication/$filename"
        else if (filename.indexOf("read_duplication.DupRate_plot.r") > 0)            "read_duplication/rscripts/$filename"
        else if (filename.indexOf("read_duplication.pos.DupRate.xls") > 0)           "read_duplication/dup_pos/$filename"
        else if (filename.indexOf("read_duplication.seq.DupRate.xls") > 0)           "read_duplication/dup_seq/$filename"
        else if (filename.indexOf("junction_annotation.junction.xls") > 0)           "junction_annotation/anno/$filename"
        else if (filename.indexOf("junction_annotation.junction_plot.r") > 0)        "junction_annotation/rscripts/$filename"
        else if (filename.indexOf("junction_annotation.splice_events.pdf") > 0)      "junction_annotation/$filename"
        else if (filename.indexOf("junction_annotation.splice_junction.pdf") > 0)    "junction_annotation/$filename"
        else if (filename.indexOf("bam_stat.txt") > 0)                               "bam_stat/$filename"
        else if (filename.indexOf("junction_saturation.junctionSaturation_plot.r") > 0)                                                                                                                            "junction_saturation/rscripts/$filename"
        else if (filename.indexOf("junction_saturation.junctionSaturation_plot.pdf") > 0)                                                                                                                          "junction_saturation/$filename"
        else filename
        }
    
    when:
    !params.skipRSeQC

    input:
    set val(name), file(bam_file) from sorted_bams_for_rseqc
    file(bam_indices) from sorted_bam_indices_for_rseqc

    output:
    file "*.{txt,pdf,r,xls}" into rseqc_results

    script:
    """
    
    export PATH=~/.local/bin:$PATH
    
    read_distribution.py -i ${bam_file} \
                         -r ${genome_refseq} \
                         > ${name}.read_distribution.txt

    read_duplication.py -i ${bam_file} \
                        -o ${name}.read_duplication

    infer_experiment.py -i ${bam_file} \
                        -r ${genome_refseq} \
                        > ${name}.infer_experiment.txt
                        
    bam_stat.py -i ${bam_file} \
                        > ${name}.bam_stat.txt
                        
    junction_annotation.py -i ${bam_file} \
                           -o ${name}.junction_annotation \
                           -r ${genome_refseq}
                           
    junction_saturation.py -i ${bam_file} \
                           -o ${name}.junction_saturation \
                           -r ${genome_refseq}                       
    """
 }

process rseqc_count {
    tag "$name"
    time '8h'
    validExitStatus 0,143
    memory '40 GB'
    publishDir "${params.outdir}/rseqc_counts" , mode: 'copy', pattern: "*FPKM.xls"
    
    when:
    params.count

    input:
    set val(name), file(bam_file) from sorted_bams_for_rseqc_count
    file(bam_indices) from sorted_bam_indices_for_rseqc_count

    output:
    file "*FPKM.xls" into rseqc_counts

    script:
    """
    export PATH=~/.local/bin:$PATH
    
    FPKM_count.py -i ${bam_file} \
                           -o ${name}.counts \
                           -r ${genome_refseq} \
                           -e
    """
 }



/*
 *STEP 5c - Analyze coverage using pileup.sh
 */

process pileup {
    tag "$name"
    memory '50 GB'
    publishDir "${params.outdir}/qc/pileup", mode: 'copy', pattern: "*.txt"

    input:
    set val(name), file(bam_file) from sorted_bams_for_pileup
    file(bam_indices) from sorted_bam_indicies_for_pileup

    output:
    file("*.txt") into pileup_results

    script:
    """

    pileup.sh -Xmx20g \
              in=${bam_file} \
              out=${name}.coverage.stats.txt \
              hist=${name}.coverage.hist.txt
    """
 }

/*
 *STEP 6a - Create non-normalzied bedGraphs for analysis using FStitch/Tfit
 */

process bedgraphs {
    validExitStatus 0,143
    tag "$name"
    memory '80 GB'
    time '4h'
    if (params.savebg) {
            publishDir "${params.outdir}/mapped/bedgraphs", mode: 'copy', pattern: "${name}.bedGraph"
    } 

    input:
    set val(name), file(bam_file) from sorted_bams_for_bedtools_bedgraph
    set val(name), file(bam_indices) from sorted_bam_indices_for_bedtools_bedgraph
    set val(name), file(millions_mapped) from bam_milmapped_bedgraph

    output:
    set val(name), file("*pos.bedGraph") into pos_non_normalized_bedgraphs
    set val(name), file("*neg.bedGraph") into neg_non_normalized_bedgraphs
    set val(name), file("${name}.bedGraph") into non_normalized_bedgraphs
    set val(name), file("${name}.rcc.bedGraph") into bedgraph_tdf
    set val(name), file("${name}.pos.rcc.bedGraph") into bedgraph_bigwig_pos
    set val(name), file("${name}.neg.rcc.bedGraph") into bedgraph_bigwig_neg

    script:
    """

    genomeCoverageBed \
                     -bg \
                     -strand + \
                     -g hg38 \
                     -ibam ${bam_file} \
                     -split \
                     > ${name}.pos.bedGraph

    genomeCoverageBed \
                     -bg \
                     -strand - \
                     -g hg38 \
                     -ibam ${bam_file} \
                     -split \
                     > ${name}.tmp.neg.bedGraph

     awk 'BEGIN{FS=OFS="\t"} {\$4=-\$4}1' ${name}.tmp.neg.bedGraph \
        > ${name}.neg.bedGraph
        rm ${name}.tmp.neg.bedGraph

    cat ${name}.pos.bedGraph \
        ${name}.neg.bedGraph \
        > ${name}.unsorted.bedGraph

    sortBed \
             -i ${name}.unsorted.bedGraph \
             > ${name}.bedGraph

    rm ${name}.unsorted.bedGraph

    python ${params.rcc} \
        ${name}.bedGraph \
        ${millions_mapped} \
        ${name}.rcc.bedGraph \

    python ${params.rcc} \
        ${name}.pos.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.pos.rcc.bedGraph

    sortBed -i ${name}.unsorted.pos.rcc.bedGraph > ${name}.pos.rcc.bedGraph
    rm ${name}.unsorted.pos.rcc.bedGraph

    python ${params.rcc} \
        ${name}.neg.bedGraph \
        ${millions_mapped} \
        ${name}.unsorted.neg.rcc.bedGraph

    sortBed -i ${name}.unsorted.neg.rcc.bedGraph > ${name}.neg.rcc.bedGraph
    rm ${name}.unsorted.neg.rcc.bedGraph

    """
 }

/*
 *STEP 7 - Normalize bigWigs by millions of reads mapped for visualization
 */

process normalized_bigwigs {
    validExitStatus 0
    tag "$name"
    memory '30 GB'
    publishDir "${params.outdir}/mapped/rcc_bigwig", mode: 'copy'
    
    when:
    params.savebw

    input:
    set val(name), file(neg_bedgraph) from bedgraph_bigwig_neg
    set val(name), file(pos_bedgraph) from bedgraph_bigwig_pos
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.rcc.bw") into normalized_bigwig

    script:
    """
    ${params.bedGraphToBigWig} ${pos_bedgraph} ${chrom_sizes} ${name}.pos.rcc.bw
    ${params.bedGraphToBigWig} ${neg_bedgraph} ${chrom_sizes} ${name}.neg.rcc.bw

    """
}

/*
 *STEP 8 - IGV Tools : generate tdfs for optimal visualization in Integrative Genomics Viewer (IGV)
 */

process igvtools {
    tag "$name"
    memory '200 GB'
    time '1h'
    // This often blows up due to a ceiling in memory usage, so we can ignore
    // and re-run later as it's non-essential.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/mapped/tdfs", mode: 'copy', pattern: "*.tdf"

    input:
    set val(name), file(normalized_bg) from bedgraph_tdf
    file(chrom_sizes) from chrom_sizes

    output:
    set val(name), file("*.tdf") into tiled_data_ch

    script:
    """
    igvtools toTDF ${normalized_bg} ${name}.rcc.tdf ${chrom_sizes}
    """
 }



/*
 * STEP 9 - MultiQC
 */
process multiQC {
    validExitStatus 0,1,143
    errorStrategy 'ignore'
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "multiqc_report.html"
    publishDir "${params.outdir}/multiqc/", mode: 'copy', pattern: "*_data"

    when:
    !params.skipMultiQC

    input:
    file multiqc_config
    file (fastqc:'qc/fastqc/*') from fastqc_results.collect()
    file ('qc/fastqc/*') from trimmed_fastqc_results.collect()
    file ('qc/trimstats/*') from trim_stats.collect()
    file ('qc/mapstats/*') from bam_flagstat.collect()
    file ('qc/rseqc/*') from rseqc_results.collect()
    file ('qc/preseq/*') from preseq_results.collect()
    file ('software_versions/*') from software_versions_yaml
    file ('qc/hisat2_mapstats/*') from hisat2_mapstats.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data" into multiqc_report_files

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''

    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config
    """
}

/*
 * STEP 10 - Output Description HTML
 */
//
//process output_documentation {
//    tag "$prefix"
//    publishDir "${params.outdir}/pipeline_info/", mode: 'copy'
//
//    input:
//    file output_docs
//
//    output:
//    file "results_description.html"
//
//    script:
//    """
//    markdown_to_html.r $output_docs results_description.html
//    """
//}
//


/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[SteadyFlow] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[SteadyFlow] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[SteadyFlow] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[SteadyFlow] Pipeline Complete"

}
