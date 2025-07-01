#! /bin/bash

#SBATCH --partition=Draco
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10-00:00:00
#SBATCH --mem=350GB
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
# function for simulating illumina short paired-end reads from a genome
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
# slurm-4875865.out

halStats evaluation/sandbox/pangenome/sorghum-pg.full.hal
# slurm-4998990.out

COMMENT

#=================================================================================
# activate conda environment
# activate ./01_graph_minigraph-cactus.slurm script
# create base directory
#=================================================================================

module load anaconda3/2023.09
conda activate vg_env2

source ./01_graph_minigraph-cactus.slurm

base="evaluation/sandbox"

#=================================================================================
# recreate minimizer index (that includes zipcodes) and a new zipcodes file
# map some reads
#=================================================================================

: << 'COMMENT'

map_reads "$base/simulated_reads/riouncc1.fq" "$base/pangenome" "$base/mappings"

COMMENT

# slurm-4999007.out

#=================================================================================
# map other reads
#=================================================================================

: << 'COMMENT'

map_reads2 "$base/simulated_reads/simulated1.fq" "$base/pangenome" "$base/mappings"

COMMENT

# slurm-4999866.out

#=================================================================================
# stats of alignments (gam file)
#=================================================================================

: << 'COMMENT'

find $base/mappings/*gam \
	| parallel -j 2 --joblog $pjl/vg_stats_eval_sandbox.txt vg_stats {}

COMMENT

# slurm-5012567.out

#=================================================================================
# SV calling and genotyping
#=================================================================================

# compute the snarls
#vg snarls $base/pangenome/sorghum-pg.d2.gbz > $base/pangenome/sorghum-pg.d2.snarls
# slurm-5012576.out


: << 'COMMENT'

find $base/mappings/*gam \
	| parallel -j 2 --joblog $pjl/vg_call_eval_sandbox.txt \
	vg_call {} "$base/vg_call" "$base/pangenome"

COMMENT

#slurm-5012581.out

#=================================================================================
# filter VCF file, exclude rows with lowdepth and lowad and rows with SVs less than 50bp
#=================================================================================

: << 'COMMENT'

#rm $base/vg_call/*filtered.vcf

find $base/vg_call/*vcf \
       | parallel -j 2 --joblog $pjl/filter_vcf_eval_sandbox.txt filter_vcf {} "$base/vg_call/"

#slurm-5026737.out

COMMENT

#=================================================================================
# map simulated reads to original pangenome that was constructed without the simulated genome
# genotype mapped reads, filter it
#=================================================================================

# map_reads2 "$base/simulated_reads/simulated1.fq" "sorghum-pg" "$base/mappings/mapped_to_original_pangenome"
# ==> slurm-5026730.out

# vg_call "$base/mappings/mapped_to_original_pangenome/simulated.mapped.gam" \	
# 	"$base/vg_call/mapped_to_original_pangenome" "sorghum-pg"
# ==> slurm-_______.out

# rm $base/vg_call/mapped_to_original_pangenome/*.filtered.vcf
# filter_vcf $base/vg_call/mapped_to_original_pangenome/*.vcf "$base/vg_call/mapped_to_original_pangenome"
# ==> slurm-_______.out

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
