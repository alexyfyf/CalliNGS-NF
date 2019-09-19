/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */
 
  
/* 
 * 'CalliNGS-NF' - A Nextflow pipeline for variant calling with NGS data
 * 
 * This pipeline that reproduces steps from the GATK best practics of SNP 
 * calling with RNAseq data procedure:
 * https://software.broadinstitute.org/gatk/guide/article?id=3891
 * 
 * Anna Vlasova 
 * Emilio Palumbo 
 * Paolo Di Tommaso
 * Evan Floden 
 */
 
nextflow.preview.dsl=2

/*
 * Define the default parameters
 */ 

params.genome     = "$baseDir/data/genome.fa"
params.variants   = "$baseDir/data/known_variants.vcf.gz"
params.blacklist  = "$baseDir/data/blacklist.bed" 
params.reads      = "$baseDir/data/reads/rep1_{1,2}.fq.gz"
params.results    = "results"
params.gatk       = '/usr/local/bin/GenomeAnalysisTK.jar'
params.gatk_launch = "java -jar $params.gatk" 

log.info """\
C A L L I N G S  -  N F    v 2.0 
================================
genome   : $params.genome
reads    : $params.reads
variants : $params.variants
blacklist: $params.blacklist
results  : $params.results
gatk     : $params.gatk
"""

include './callings' params(params)

workflow {
  reads_ch = Channel.fromFilePairs(params.reads)

  p1A_prepare_genome_samtools(params.genome)
  p1B_prepare_genome_picard(params.genome)
  p1C_prepare_star_genome_index(params.genome)
  p1D_prepare_vcf_file(
    params.variants, 
    params.blacklist
  )
  p2_rnaseq_mapping_star(
    params.genome, 
    p1C_prepare_star_genome_index.out, 
    reads_ch
  )
  p3_rnaseq_gatk_splitNcigar(
    params.genome, 
    p1A_prepare_genome_samtools.out, 
    p1B_prepare_genome_picard.out, 
    p2_rnaseq_mapping_star.out
  )
  p4_rnaseq_gatk_recalibrate(
    params.genome,
    p1A_prepare_genome_samtools.out,
    p1B_prepare_genome_picard.out,
    p3_rnaseq_gatk_splitNcigar.out,
    p1D_prepare_vcf_file.out
  )
  p5_rnaseq_call_variants(
    params.genome,
    p1A_prepare_genome_samtools.out,
    p1B_prepare_genome_picard.out,
    p4_rnaseq_gatk_recalibrate.out.groupTuple()
  )
  p6A_post_process_vcf(
    p5_rnaseq_call_variants.out,
    p1D_prepare_vcf_file.out
  )
  p6B_prepare_vcf_for_ase(
    p6A_post_process_vcf.out
  ) 
  p6C_ASE_knownSNPs(
    params.genome,
    p1A_prepare_genome_samtools.out,
    p1B_prepare_genome_picard.out,
    groupBamVcf(
      p4_rnaseq_gatk_recalibrate.out, 
      p6A_post_process_vcf.out
    )
  )
}