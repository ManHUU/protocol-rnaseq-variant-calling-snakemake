# Protocol for a Snakemake-Based Workflow for RNA-Seq Variant Calling Using GATK on Slurm-based HPC Clusters

This repository provides a **fully documented and reproducible computational workflow** for robust variant detection from RNA sequencing (RNA-seq) data.  
The pipeline integrates **STAR**, **SAMtools**, **GATK**, and **BCFtools** using the **Snakemake** workflow management system, enabling scalable and reproducible execution on **Slurm-based High-Performance Computing (HPC) clusters**.

This repository accompanies the manuscript:

**Protocol for a Snakemake-Based Workflow for RNA-Seq Variant Calling Using GATK on Slurm-based HPC Clusters**

---

## Background and Motivation

The increasing volume and affordability of RNA sequencing data provide new opportunities to detect genetic variants beyond traditional differential gene expression analysis. However, variant calling from RNA-seq data remains challenging due to transcript splicing, uneven coverage, and technical biases.

Recent advances in splice-aware aligners, variant calling algorithms, and curated reference resources (e.g., dbSNP, ClinVar, gnomAD) have significantly improved the reliability and interpretability of RNA-derived variant detection.  
This protocol presents a **practical, modular, and reproducible** workflow that adheres to widely accepted best practices while remaining adaptable to diverse HPC environments.

---

## Workflow Overview

| Step | Stage | Description | Tools |
|-----:|------|-------------|-------|
| 1 | Splice-aware alignment | Align RNA-seq reads to the reference genome | STAR |
| 2 | Post-alignment processing | Sort, index BAM files and mark duplicates | SAMtools, GATK |
| 3 | Variant calling | RNA-seq–aware variant detection | GATK |
| 4 | Variant filtering | Quality control and hard filtering | GATK, BCFtools |
| 5 | Visualization & reporting (optional) | Variant statistics and summary plots | BCFtools / custom scripts |

---

## Software Requirements

| Category | Software |
|--------|----------|
| Workflow engine | Snakemake |
| Alignment | STAR |
| BAM processing | SAMtools |
| Variant calling | GATK |
| Variant filtering | BCFtools |
| Environment management | Conda / Mamba|
| Cluster scheduler | Slurm |

All software dependencies are defined using reproducible environments located in `envs/`.

---

## Repository Structure

```text
.
├── README.md
├── LICENSE
├── workflow/
│   ├── Snakefile
│   ├── snakemake_rule/
│   │   └── Snakemake_rules_all.smk
│   ├── envs/
│   │   ├── star_aligner.yaml
│   │   └── gatk.yaml
│   ├── configs/
│   │   ├── config.yaml
│   │   └── resources.yaml
│   ├── download_refs.slurm
│   ├── build_star_index.slurm
│   └── run_smk.sh
├── docs/
│   ├── protocol.md
│   ├── input_output.md
│   └── troubleshooting.md
└── tests/
    └── minimal_test_config.yaml
```

## A quick start
### Step 1: Preparation

#### 1.1 Required Reference Files

| File type | Description |
|----------|-------------|
| Reference genome FASTA | Genome sequence (`.fa`) |
| Genome annotation | Gene annotation file (`.gtf`) |
| Known variants (optional) | dbSNP, ClinVar, gnomAD VCF files |

---

#### 1.2 Download Reference Files

A Slurm helper script is provided to download the reference genome and annotation files:

```bash
sbatch download_refs.slurm
```

#### 1.3 Build STAR Genome Index
```bash
sbatch build_star_index.slurm
```
#### 1.4 Configure Workflow Parameters
Edit the configuration files located in `configs/`:

| File | Purpose |
|------|---------|
| `config.yaml` | Input/output paths, reference files, pipeline options |
| `resources.yaml` | CPU, memory, and runtime requirements |

### Step 2: Execute the Snakemake Workflow
#### 2.1 Dry Run (Recommended)
Validate the workflow without executing any jobs:
```bash
bash run_smk.sh
```

#### 2.2 Run on Slurm-based HPC Cluster
```bash
bash run_smk.sh
```

## Citation

If you use this workflow, please cite the accompanying manuscript:
Protocol for a Snakemake-Based Workflow for RNA-Seq Variant Calling Using GATK on Slurm-based HPC Clusters





