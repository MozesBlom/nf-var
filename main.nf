#!/usr/bin/env nextflow

log.info """\

        ================================================
        NF-variant - Variant calling for whole-genome resequence data
        https://github.com/MozesBlom/nf-variant    
        Author: Mozes P.K. Blom
        ================================================
        |indivfile    : ${params.indivs_file}
        |chromos      : ${params.chromos_file}
        |reference    : ${params.ref_file}
        |bams_fn      : ${params.bams_fn}
        |bamsdir      : ${params.inputdir}
        |outputdir    : ${params.outputdir}
        |
        |Sex chromos  : ${params.sex_chromos}
        |
        |Bam file needs indexing?       : ${params.index_bam_needed}
        |Obtain coverage stats?         : ${params.calc_coverage}
        |Call variants by individual?   : ${params.indiv_var_call}
        |Call variants by population?   : ${params.pops_var_call}
        |Calcuate missing data?         : ${params.calc_missing_data}
        |Call consensus sequences?      : ${params.call_consensus}
        |
        |Variant calling & filtering
        |---------------------------
        |Mask heterozygous positions?   : ${params.mask_hets}
        |Mask low/high coverage sites?  : ${params.mask_cov}
        |Min depth per position?        : ${params.mask_min_cov}
        |Max depth per position?        : ${params.mask_max_cov}
        |
        ================================================
        """
        .stripIndent()

/*
====================================================
~ ~ ~ > *  Input channels and file check  * < ~ ~ ~ 
====================================================
*/

/*
 * Create two lists:
 * - All individuals to include (expects that the bam files are called: indiv.bam)
 * - All chromosomes to include (must match exactly with reference genome fa ids)
 */

indivs = file(params.indivs_file).readLines()
indivs_ch = Channel.fromList(indivs)
indivs_bam_ch = indivs_ch.map{ it ->
                                    def indiv = it
                                    def bam_fn = file("${params.inputdir}/${it}${params.bams_fn}")
                                    
                                    [indiv, bam_fn]

                                    }
chromos = file(params.chromos_file).readLines()
chromos_ch = Channel.fromList(chromos)

/*
 * Set data channels:
 * - Reference genome against which the individuals have been mapped
 */
ref_ch = Channel.fromPath(params.ref_file)
                .ifEmpty { error "No reference sequence found from: ${params.ref_file}" }


/*
===============================
~ ~ ~ > *  Processes  * < ~ ~ ~ 
===============================
*/

process index_bam  {

/*
 * If need be index each of the bam files before proceeding.
 */

    tag "Index bam file"

    input:
    tuple val(indiv), path(indiv_bam)

    output:
    tuple val(indiv), path(indiv_bam), path("${indiv_bam}.bai")


    script:
    """

    samtools index -@ ${task.cpus} ${indiv_bam}

    """
}


process cov_estimate {

/*
 * Estimate the mean coverage for each individual and each of the focal chromosomes/scaffolds.
 * NOTE: I tried to lump all chromos per individual into one job but it didn't work.
 * so for now it will result in a larger number of short jobs (which is not optimal for a cluster)
 */

    tag "Estimate coverage"

    input:
    tuple val(indiv), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(indiv), path("${chromo}.tsv")


    script:
    def out_fn = "${chromo}.tsv"
    """

    samtools coverage --region ${chromo} --output ${out_fn} ${indiv_bam}

    """
}


process cov_summary_INDIV {

/*
 * Accurately summarise the coverage across all target chromos per individual
 */

    tag "Coverage summary per individual"
    publishDir "${params.outputdir}/00.coverage/", mode:'copy'

    input:
    tuple val(indiv), path(chromo_cov_tsv_list)

    output:
    path("${indiv}_coverage.tsv")


    script:
    """

    00_cov_stats_INDIV.py \
    -c ${chromo_cov_tsv_list} \
    -i ${indiv}

    """
}


process cov_summary_ALL {

/*
 * Accurately summarise the coverage across all target chromos per individual
 */

    tag "Coverage summary for all indivs"
    publishDir "${params.outputdir}/00.coverage/", mode:'copy'

    input:
    path(indivs_cov_tsv_list)

    output:
    file('*')


    script:
    """

    01_cov_stats_SUMMARY.py \
    -c ${indivs_cov_tsv_list} \
    -s ${params.sex_chromos} \
    -p ${params.plots_per_row}

    """
}


process call_variants_CHROMO {

/*
 * Call variants for each individual by chromosome and an initial round of filtering takes place.
 *
 * We use vcfallelicprimitives to deconstruct mnps to snps.
 * My original script has the -g flag after vcf ap. Double check if this is still supported, since it's not listed in the -h list
 */

    tag "Coverage summary for all indivs"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'
    label 'Endurance'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), path(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt.vcf.gz")


    script:
    """

    freebayes -f ${reference} --region ${chromo} -m 10 -p 2 ${indiv_bam}| \
    vcffilter -f "QUAL < 20" -f "( AB > 0 ) & ( AB < 0.2 )" --invert --or | \
    vcfallelicprimitives -k -g | \
    bgzip -c > ${chromo}_vars_filt.vcf.gz

    """
}


process remove_indels {

/*
 * Optional: Remove indel variation from vcf file
 *
 */

    tag "Remove indels"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_vars_filt_indels.vcf.gz")


    script:
    """

    bgzip -d -c ${var_vcf} | \
    vcffilter -f 'TYPE = ins' -f 'TYPE = del' -f 'TYPE = complex' --invert --or | \
    bgzip -c > ${chromo}_vars_filt_indels.vcf.gz
    
    """
}


process mask_hets {

/*
 * Optional: Create mask file for heterozygous sites
 *
 */

    tag "Generate mask for het sites"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

    input:
    tuple val(individual), val(chromo), path(var_vcf)

    output:
    tuple val(individual), val(chromo), path("${chromo}_hets.tsv")


    script:
    """

    bgzip -d -c ${var_vcf} | \
    vcffilter -f '( AF < 1 ) & ( AB < 0.8 )' | \
    cut -f 1,2 > ${chromo}_hets.tsv

    """
}


process mask_cov {

/*
 * Optional: Create mask file for low coverage sites
 *
 */

    tag "Generate mask for low or excess coverage sites"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

    input:
    tuple val(individual), path(indiv_bam), path(indiv_bam_bai), val(chromo), file(reference)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov.tsv")


    script:
    """

    samtools depth -aa -Q 10 -r ${chromo} ${indiv_bam} | \
    awk '(\$3 < ${params.mask_min_cov} || \$3 > ${params.mask_max_cov}) {print \$1,\$2}' > ${chromo}_cov.tsv

    """
}


process mask_merge {

/*
 * Optional: Merge two mask files
 *
 */

    tag "Merge low cov and het mask files"
    publishDir "${params.outputdir}/01.variants/${individual}", mode:'copy'

    input:
    tuple val(individual), val(chromo), path(het_bed), path(cov_bed)

    output:
    tuple val(individual), val(chromo), path("${chromo}_cov_hets.tsv")


    script:
    """

    sed -i '/#CHROM/d' ${cov_bed}
    sed -i '/##/d' ${het_bed}
    sed -i 's/ /\t/g' ${cov_bed}
    sed -i 's/ /\t/g' ${het_bed}

    # Merge two mask files and sort by position. Remove duplicate entries
    cat ${cov_bed} ${het_bed} | \
    sort -Vk1 -Vk2 | \
    uniq > ${chromo}_cov_hets.tsv

    # Add a header for each column to make it a proper bed file. NOTE, for some reason bcftools consensus has a problem with BED naming, hence tsv.
    sed -i '1i #CHROM\tPOS' ${chromo}_cov_hets.tsv

    """
}


process call_consensus {

/*
 * Optional: Call a consensus sequence for each individual by chromosome WITHOUT MASKING (not recommended)
 */

    tag "Call consensus without masking"
    publishDir "${params.outputdir}/02.consensus/${individual}", mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")


    script:
    """
    tabix -p vcf ${vcf_fn}

    samtools faidx ${params.ref_file} ${chromo} | 
    bcftools consensus ${vcf_fn} -o ${individual}_${chromo}_cons.fa

    sed -i 's/${chromo}/${individual}/g' ${individual}_${chromo}_cons.fa

    """
}


process call_consensus_MASK {

/*
 * Optional: Call a consensus sequence for each individual by chromosome WIT MASKING (recommended)
 */

    tag "Call consensus with masking"
    publishDir "${params.outputdir}/02.consensus/${individual}", mode:'copy'

    input:
    tuple val(individual), val(chromo), path(vcf_fn), path(mask_fn)

    output:
    tuple val(individual), val(chromo), path("${individual}_${chromo}_cons.fa")


    script:
    """
    tabix -p vcf ${vcf_fn}

    samtools faidx ${params.ref_file} ${chromo} | 
    bcftools consensus ${vcf_fn} -m ${mask_fn} -o ${individual}_${chromo}_cons.fa

    sed -i 's/${chromo}/${individual}/g' ${individual}_${chromo}_cons.fa

    """
}


process calc_missing_data_INDIV {

/*
 * Optional: Estimate the amount of missing data for each individual and each chromosome
 * NOTE: Within the framework of the current pipeline this can be optimised by using the mask files.
 * However, I already had existing scripts borrowed from nf-phylo to make this work.
 */

    tag "Calculate missing data per individual and chromosome"
    publishDir "${params.outputdir}/02.consensus/${individual}", mode:'copy'

    input:
    tuple val(individual), path(cons_fn_list)

    output:
    file("${individual}_missing_data.tsv")


    script:
    """

    02_cons_stats_INDIV.py \
    -c ${cons_fn_list} \
    -i ${individual}

    """
}


process calc_missing_data_SUMMARY {

/*
 * Optional: Summarise the amount of missing data across all individuals
 * NOTE: Within the framework of the current pipeline this can be optimised by using the mask files.
 * However, I already had existing scripts borrowed from nf-phylo to make this work.
 */

    tag "Calculate missing data across all individual"
    publishDir "${params.outputdir}/02.consensus/", mode:'copy'

    input:
    path(cons_fn_list)

    output:
    file("all_indivs_missing_data.tsv")
    file("all_indivs_missing_data.pdf")


    script:
    """
    unset DISPLAY

    03_cons_stats_SUMMARY.py \
    -c ${cons_fn_list}

    """
}

/*
============================================
~ ~ ~ > *  Pipeline specification  * < ~ ~ ~ 
============================================
*/

workflow {
    // If need be index the bam file and create a standard channel for downstream workflow:
    // - indiv, bam, bai, chromo, reference (the .bai was added because otherwise it wouldn't be copied to the work folder)
    if (params.index_bam_needed) {
        indivs_bam_indexed_ch = index_bam(indivs_bam_ch)
        indivs_bam_chromo_ref_ch = indivs_bam_indexed_ch.combine(chromos_ch).combine(ref_ch)
    } else {
        indivs_bam_indexed_ch = indivs_bam_ch.map{ it ->
                                    def indiv = it[0]
                                    def bam_fn = it[1]
                                    def bai_fn = file("${bam_fn}.bai")
                                    
                                    [indiv, bam_fn, bai_fn]

                                    }
        indivs_bam_chromo_ref_ch = indivs_bam_indexed_ch.combine(chromos_ch).combine(ref_ch)
    }

    // Estimate and visualise coverage stats
    if (params.calc_coverage) {
        cov_estimate(indivs_bam_chromo_ref_ch)
        cov_estimate.out.groupTuple() | cov_summary_INDIV
        cov_summary_INDIV.out.collect() | cov_summary_ALL
    } else {
    }

    // Call variants and filter (with or without indels)
    if (params.indiv_var_call == true && params.filt_indels == true) {
        vars_filt_fn_ch = call_variants_CHROMO(indivs_bam_chromo_ref_ch) | remove_indels
           // This 'else if' rather than 'else' statement is included in case you only want to estimate the coverage
    } else if (params.indiv_var_call == true && params.filt_indels == false) {
        vars_filt_fn_ch = call_variants_CHROMO(indivs_bam_chromo_ref_ch)
    }

    // Create mask files (optional, but highly recommended). This can only be run in combination with variant calling
    if (params.indiv_var_call == true && params.mask_hets == true && params.mask_cov == true) {
        mask_het_fn_ch = vars_filt_fn_ch | mask_hets
        mask_cov_fn_ch = mask_cov(indivs_bam_chromo_ref_ch)
        mask_combined_ch = mask_het_fn_ch.combine(mask_cov_fn_ch, by: [0,1])
        mask_fn_ch = mask_merge(mask_combined_ch)
    } else if (params.indiv_var_call == true && params.mask_hets == true && params.mask_cov == false) {
        mask_fn_ch = vars_filt_fn_ch | mask_hets
    } else if (params.indiv_var_call == true && params.mask_hets == false && params.mask_cov == true) {
        mask_fn_ch = mask_cov(indivs_bam_chromo_ref_ch)
    }

    // Call consensus (optional, this can only be run in combination with indiv variant calling)
    if (params.indiv_var_call == true && params.call_consensus == true) {
        if (params.mask_hets == true || params.mask_cov == true) {
            vars_mask_ch = vars_filt_fn_ch.combine(mask_fn_ch, by: [0,1])
            consensus_fn = call_consensus_MASK(vars_mask_ch)
        } else {
            consensus_fn = call_consensus_MASK(vars_filt_fn_ch)
        }
        // Calculate missing data?
        if (params.calc_missing_data == true) {
            consensus_fn.groupTuple(by: 0)
                        .map{ it ->

                                def indiv = it[0]
                                def cons_fn = it[2].flatten()

                                [indiv, cons_fn]
                            } | calc_missing_data_INDIV
            calc_missing_data_INDIV.out.collect() | calc_missing_data_SUMMARY
        }
    }
}


