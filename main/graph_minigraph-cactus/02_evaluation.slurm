#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=10-00:00:00
#SBATCH --mem=400GB
#SBATCH --mail-user=cnnamdi@charlotte.edu
#SBATCH --mail-type=ALL

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID : $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "======================================================"
echo ""

module load parallel/20240422
pjl="parallel_jobs_log"

#=================================================================================
# simulating structural variation
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate survivor

cd evaluation/sandbox/

# SURVIVOR simSV parameter_file

SURVIVOR simSV riouncc.fa parameter_file 0 0 simulated

# SURVIVOR eval simulated.vcf simulated.bed 10 eval_res

COMMENT

# slurm-4875594.out

#=================================================================================
# Function for simulating illumina short paired-end reads from a genome
#=================================================================================

simulate_short_reads() {
	echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $1
        echo $sample
        echo "--------------------------------------"

	art_illumina -ss HS25 -p -sam -i $1 \
		-l 150 -f 45 -m 500 -s 50 \
		-o $2/simulated_reads/$sample
	}

export -f simulate_short_reads

#=================================================================================
# simulate illumina short paired-end reads
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate art

find evaluation/sandbox/ -type f \( -name "riouncc.fa" -o -name "simulated.fasta" \) \
	| parallel -j 2 --joblog $pjl/sandbox_simulate_short_reads.txt simulate_short_reads {} "evaluation/sandbox"

COMMENT

# slurm-4875784.out
# slurm-4876026.out

#=================================================================================
# pangenome graph construction using minigraph-cactus
#=================================================================================

: << 'COMMENT'

cactus_dir="/projects/cooper_research2/chinaza/packages/cactus/cactus-bin-v2.9.3/"
source "${cactus_dir}venv-cactus-v2.9.3/bin/activate"
#cactus --version

#toil clean jobstore

cactus-pangenome jobstore \
	evaluation/sandbox/pangenome/seq_file.txt \
        --outDir evaluation/sandbox/pangenome \
        --outName sorghum-pg \
        --reference riouncc \
        --vcf --giraffe --gfa --gbz
# slurm-_________.out


#halStats evaluation/sandbox/pangenome/sorghum-pg.full.hal
# slurm-_________.out

COMMENT

#=================================================================================
# activate conda environment
#=================================================================================

module load anaconda3/2023.09
conda activate vg_env2

source ./01_graph_minigraph-cactus.slurm

#=================================================================================
# recreate minimizer index (that includes zipcodes) and a new zipcodes file
#=================================================================================

: << 'COMMENT'

base="evaluation/sandbox"
map_reads "$base/simulated_reads/riouncc1.fq" "$base/pangenome" "$base/mappings"

COMMENT

#=================================================================================
# map other reads
#=================================================================================

: << 'COMMENT'



COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'
COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'
COMMENT

#=================================================================================
#
#=================================================================================

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
