#!/usr/bin/bash

##----------------- Define default Jobs-----------------
# Ensure jobs is an integer; default to 1 if first argument is empty
jobs=${1:-1}

# Re-verify that jobs is a number (HPC users often pass non-integers by mistake)
if ! [[ "$jobs" =~ ^[0-9]+$ ]]; then
    echo "Warning: Job count '$jobs' is not an integer. Defaulting to 1."
    jobs=1
fi



##----------------- Define variables-----------------
# Replace below  with your actual account and path
ACCOUNT="ausei14399"
PARTITION="genoa"
# Use a completely distinct path name
CONDA_PREFIX_DIR="/gpfs/work2/0/prjs1616/gatk_paper_demo/snakemake_conda_factory"
MAIN_CONFIG="/gpfs/work2/0/prjs1616/gatk_paper_demo/workflow/configs/config.yaml"
CLUSTER_CONFIG="/gpfs/work2/0/prjs1616/gatk_paper_demo/workflow/configs/resources.yaml"
SNAKE_DIR="/projects/prjs1616/gatk_paper_demo/workflow/snakemake_rule/Snakemake_rules_all.smk"



# Clean slate
# Only run below mannually if you actually want to reset your environments.
#rm -rf "$CONDA_PREFIX_DIR"
mkdir -p "$CONDA_PREFIX_DIR"

echo "STEP 1: Unlocking"
snakemake --snakefile "$SNAKE_DIR" --configfile "$MAIN_CONFIG" --unlock --cores 1

echo "STEP 2: Building environments with CONDA frontend (if missing)"
# We use 'conda' instead of 'mamba' here to bypass the strict libmamba check
srun --account=$ACCOUNT \
     --partition=$PARTITION \
     --time=00:45:00 \
     --mem=16G \
     --ntasks=1 \
     snakemake --snakefile "$SNAKE_DIR" \
     --configfile "$MAIN_CONFIG" \
     --use-conda \
     --conda-frontend conda \
     --conda-prefix "$CONDA_PREFIX_DIR" \
     --conda-create-envs-only \
     --cores 1


# Actual Execution
echo "STEP 3: Executing"
snakemake --verbose \
    --snakefile "$SNAKE_DIR" \
    --configfile "$MAIN_CONFIG" \
    --cluster-config "$CLUSTER_CONFIG" \
    --cluster "sbatch --account={cluster.account} \
                      --partition={cluster.partition} \
                      --time={cluster.time} \
                      --cpus-per-task={cluster.cpus_per_task} \
                      --mem={cluster.mem} \
                      --output={cluster.output} \
                      --error={cluster.error} \
                      --job-name={cluster.job_name}" \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix "$CONDA_PREFIX_DIR" \
    --jobs "$jobs" \
    --cores "$jobs" \
    --latency-wait 120 \
    --keep-going \
    --rerun-incomplete



