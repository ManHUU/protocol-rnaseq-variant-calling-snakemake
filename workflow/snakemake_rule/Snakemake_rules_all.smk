# GATK RNA-seq Variant Calling Pipeline - Snakemake Workflow
# Optimized with resource parameters aligned to resources.yaml

import os
from pathlib import Path
import glob

# Load configuration
configfile: "config.yaml"


# Define directory
OUTPUT_DIR = config["output_dir"]
LOG_DIR = config["log_dir"]
TEMP_DIR = config["temp_dir"]
ENV_DIR = config["envs"]
FASTQ_DIR = config["fastq_dir"]

# Create directories
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)
os.makedirs(TEMP_DIR, exist_ok=True)


CHROMOSOMES = config["chromosomes"]


# Define reference paths for clarity
REF_FASTA = config["reference"]["fasta"]
REF_DICT = REF_FASTA.replace(".fa", ".dict").replace(".fasta", ".dict")
REF_FAI = REF_FASTA + ".fai"

# Input File Scanning
FOUND_FASTQ = glob.glob(os.path.join(FASTQ_DIR, "*.fastq"))
SAMPLES = []

for fastq_path in FOUND_FASTQ:
    filename = os.path.basename(fastq_path)
    sample_name = filename.split('.')[0]
    SAMPLES.append(sample_name)

SAMPLES.sort()



#-----
rule all:
    input:
#        expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz", sample=SAMPLES),
#        f"{OUTPUT_DIR}/filtered.joined.vcf.gz"
        REF_FAI,
        REF_DICT,
        expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES)



# HIGHLIGHT: Rule to generate .fai and .dict if missing
rule prepare_reference:
    input:
        fasta = REF_FASTA
    output:
        fai = REF_FAI,
        dict = REF_DICT
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/prepare_reference.log"
    shell:
        """
        exec &> {log}
        # Generate FAI if missing
        if [ ! -f {output.fai} ]; then
            samtools faidx {input.fasta}
        fi

        # Generate Dict if missing
        if [ ! -f {output.dict} ]; then
            gatk CreateSequenceDictionary -R {input.fasta}
        fi
        """



# Rule 1: Align FASTQ to Reference and Sort
rule fastq_to_bam:
    input:
        # Assuming single-end based on your file list; 
        # if paired-end, you'd use r1 and r2 here.
        fastq=lambda wildcards: os.path.join(FASTQ_DIR, f"{wildcards.sample}.fastq"),
        star_index = config["reference"]["star_index_dir"]
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.sorted.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.sorted.bam.bai"
    conda:
        os.path.join(ENV_DIR, "star_aligner.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.fastq_to_bam.log"
    params:
        ref=config["reference"]["fasta"],
        sample_id = "{sample}",
        out_prefix = f"{OUTPUT_DIR}/{{sample}}."
    threads:
        config["resources"]["fastq_to_bam"]["threads"]
    resources:
        mem_mb=config["resources"]["fastq_to_bam"]["memory"] * 1024
    shell:
        """
        exec &> {log}
        echo "DATE: $(date)"
        echo "RULE: fastq_to_bam"
        echo "SAMPLE: {params.sample_id}"
        echo "---------------------------------------"

        # STAR alignment
        # Generates coordinate-sorted BAM and adds Read Groups required for GATK
        STAR --runThreadN {threads} \
             --genomeDir {input.star_index} \
             --readFilesIn {input.fastq} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattrRGline ID:{params.sample_id} SM:{params.sample_id} LB:lib1 PL:ILLUMINA PU:unit1 \
             --outFileNamePrefix {params.out_prefix}

        # Rename STAR's default output name to match Snakemake's output
        mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam}

        # Index the BAM
        samtools index -@ {threads} {output.bam}
        echo "---------------------------------------"
        echo "Alignment and indexing complete."
        """       



# Rule 2: Mark Duplicates
rule mark_duplicates:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.sorted.bam"
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.markdup.bam.bai",
        metrics=f"{OUTPUT_DIR}/{{sample}}.markdup_metrics.txt"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.mark_duplicates.log"
    threads: config["resources"]["mark_duplicates"]["threads"]
    resources:
        mem_mb=config["resources"]["mark_duplicates"]["memory"] * 1024
    shell:
        """
        exec &> {log}
    
        gatk --java-options "-Xmx{resources.mem_mb}m" MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --VALIDATION_STRINGENCY SILENT > {log} 2>&1
        samtools index {output.bam}
        """



#Rule 3: Filter non-standard contigs
rule filter_standard_contigs:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.bam"
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.filter_contigs.log"
    threads: config["resources"]["filter_contigs"]["threads"]
    resources:
        mem_mb=config["resources"]["filter_contigs"]["memory"] * 1024

    shell:
        """
        exec &> {log}
        samtools view -h {input.bam} \
        | awk 'BEGIN{{OFS="\t"}} /^@/ || $3 ~ /^chr([1-9]$|1[0-9]$|2[0-2]$|X|Y|M)$/ {{print}}' \
        | samtools view -b -o {output.bam} -
        """




#Rule 4:Index "filtered bam files"
rule index_filtered_bam:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam"
    output:
        bai=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam.bai"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.index_filtered_bam.log"
    threads: config["resources"]["index_bam"]["threads"]
    resources:
        mem_mb=config["resources"]["index_bam"]["memory"] * 1024
    shell:
        """
        exec &> {log}
        samtools index {input.bam} {output.bai}
        """   



#Rule 5: SplitNCigar reads
rule split_n_cigar_reads:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam",
        bai=f"{OUTPUT_DIR}/{{sample}}.markdup.filtered.bam.bai",
        fasta=REF_FASTA,
        fai=REF_FAI,
        dict=REF_DICT
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.split.bam"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.split_n_cigar_reads.log"
    threads: config["resources"]["split_n_cigar_reads"]["threads"]
    resources:
        mem_mb=config["resources"]["split_n_cigar_reads"]["memory"] * 1024
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" SplitNCigarReads \
            -R {config[reference][fasta]} \
            -I {input.bam} \
            -O {output.bam}
        """



#Rule 6: Calculate BQSR (Base Quality Score Recalibration) table
rule base_recalibrator:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.split.bam"
    output:
        table=f"{OUTPUT_DIR}/{{sample}}.recal_data.table"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.base_recalibrator.log"
    threads: config["resources"]["bqsr"]["threads"]
    resources:
        mem_mb=config["resources"]["bqsr"]["memory"] * 1024
    params:
        reference=config["reference"]["fasta"],
        dbsnp=config["reference"]["dbsnp"],
        mills=config["reference"]["mills_indels"]
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" BaseRecalibrator \
            -R {params.reference} \
            -I {input.bam} \
            --known-sites {params.dbsnp} \
            --known-sites {params.mills} \
            -O {output.table}
        """



# Rule 7: Apply  BQSR 
rule apply_BQSR:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.split.bam",
	table=f"{OUTPUT_DIR}/{{sample}}.recal_data.table"
    output:
        bam=f"{OUTPUT_DIR}/{{sample}}.BQSR.bam"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.gather_bqsr_reports.log"
    threads: config["resources"]["bqsr"]["threads"]
    resources:
        mem_mb=config["resources"]["bqsr"]["memory"] * 1024
    params:
        reference=config["reference"]["fasta"],
        dbsnp=config["reference"]["dbsnp"],
        mills=config["reference"]["mills_indels"]

    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" ApplyBQSR \
            -R {params.reference} \
            -I {input.bam} \
	    --bqsr-recal-file {input.table} \
            -O {output.bam}
        """



#Rule 8: Variant Calling 
rule haplotype_caller:
    input:
        bam=f"{OUTPUT_DIR}/{{sample}}.BQSR.bam",
        fasta=REF_FASTA,
        fai=REF_FAI,
        dict=REF_DICT
    output:
        vcf=f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.variant_calling.log"
    threads: config["resources"]["haplotype_caller"]["threads"]
    resources:
        mem_mb=config["resources"]["haplotype_caller"]["memory"] * 1024
    params:
        reference=config["reference"]["fasta"],
        dbsnp=config["reference"]["dbsnp"],
        mills=config["reference"]["mills_indels"],
        intervals=" ".join([f"-L {c}" for c in config["chromosomes"]])
    shell:
        """
        exec &> {log}
        gatk --java-options "-Xmx{resources.mem_mb}m" HaplotypeCaller \
            -ERC GVCF \
            -R {params.reference} \
            -I {input.bam} \
            {params.intervals} \
	    --dont-use-soft-clipped-bases true \
	    --standard-min-confidence-threshold-for-calling 20.0 \
            -O {output.vcf}
        """



# Rule 9: Index gVCF files
rule index_gvcf:
    input:
        vcf=f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz"
    output:
        tbi=f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz.tbi"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/{{sample}}.index_vcf.log"
    threads: config["resources"]["index_bam"]["threads"]
    resources:
        mem_mb=config["resources"]["index_bam"]["memory"] * 1024
    shell:
        """
        set -e

        exec &> {log}
        echo "Index vcf.gz"
        gatk IndexFeatureFile \
            -I {input.vcf}
        """



#Rule 10: generate sample map for next step (create genomicsDB)
rule generate_sample_map:
    input: vcf_gz=expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz",sample=SAMPLES)
    output:
        map_file=f"{OUTPUT_DIR}/sample_map_file",
        done=f"{OUTPUT_DIR}/BARIA_sample_gz_map.done"
    run:
        import os, glob
        vcf_files = glob.glob(f"{OUTPUT_DIR}/*.g.vcf.gz")
        SAMPLES = [os.path.basename(file).replace('.g.vcf.gz', '') for file in vcf_files]

        with open(output.map_file, 'w') as f:
            for sample in SAMPLES:
               f.write(f"{sample}\t{OUTPUT_DIR}/{sample}.g.vcf.gz\n")


        # Touch the done file to signal completion
        with open(output.done, 'w') as f:
            f.write('')



#Rule 11: generate genomicsDB
rule GenomicsDB:
    input:
        gvcfs=expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz.tbi", sample=SAMPLES),
	gvcfs_tbi=expand(f"{OUTPUT_DIR}/{{sample}}.g.vcf.gz", sample=SAMPLES),
        sample_map=f"{OUTPUT_DIR}/sample_map_file"
    output:
        touch(f"{OUTPUT_DIR}/GenomicsDB.done")
    params:
        workspace=f"{OUTPUT_DIR}/genomicsdb",
        chromosomes=" ".join([f"-L {chrom}" for chrom in CHROMOSOMES])
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/genomicsDB.log"
    threads: config["resources"]["GenomicsDB"]["threads"]
    resources:
        mem_mb=config["resources"]["GenomicsDB"]["memory"] * 1024
    shell:
        """
        exec &> {log}
        # Remove the workspace if it exists, otherwise GATK will fail
        rm -rf {params.workspace}
        gatk --java-options "-Xmx{resources.mem_mb}m" GenomicsDBImport \
            --genomicsdb-workspace-path {params.workspace} \
            --sample-name-map {input.sample_map} \
	    --batch-size 3 \
	    {params.chromosomes} \
	    --reader-threads {threads} \
            && touch {output}
        """




#Rule 12: perform joint genotyping 
rule join_genotyping:
    input:
        f"{OUTPUT_DIR}/GenomicsDB.done"
    output:
        f"{OUTPUT_DIR}/Joint_all.vcf.gz"
    log:
	f"{LOG_DIR}/genotypeGVCFs.log"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    params:
        reference=config["reference"]["fasta"],
	database=f"{OUTPUT_DIR}/genomicsdb"
    threads: config["resources"]["join_caller"]["threads"]
    resources:
        mem_mb=config["resources"]["join_caller"]["memory"] * 1024
    shell:
        """
        exec &> {log}
        echo "joint genotyping"
        gatk --java-options "-Xmx{resources.mem_mb}m" GenotypeGVCFs \
           -R {params.reference} \
           -V gendb://{params.database} \
           -O {output} > {log}
        """



#Rule 13: Annotate variant filter 
## This step won't filter out any reads, the number of reads of output will be the same as input
rule annotate_variant_filter:
    input:
        f"{OUTPUT_DIR}/Joint_all.vcf.gz"
    output:
        f"{OUTPUT_DIR}/Filtered.vcf.gz"
    conda:
        os.path.join(ENV_DIR, "gatk.yaml")
    log:
        f"{LOG_DIR}/variant_filtration.log"
    params:
        reference=config["reference"]["fasta"],
        filter_expression=config["variant_filters"]["expression"],
        filter_name=config["variant_filters"]["name"]
    resources:
        mem_mb= 64 * 1024
    shell:
        """
        exec &> {log}
        echo "variant_filtration"
        /home/mhu/miniconda3/envs/gatk_env/bin/gatk --java-options "-Xmx{resources.mem_mb}m" VariantFiltration \
            -R {params.reference} \
            -V {input} \
            -O {output} \
            --filter-name "{params.filter_name}" \
            --filter-expression "{params.filter_expression}"
       """
