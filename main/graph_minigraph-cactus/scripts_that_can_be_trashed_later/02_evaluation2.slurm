#! /bin/bash

#SBATCH --partition=Orion
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

base="/projects/cooper_research2/chinaza/"
cd $base/structural_variant_calling/graph_minigraph-cactus/evaluation2/sandbox/simulated_svs_n_reads/

# copy insertion sequence from quickstart, might comeback and generate something else that's more relevant for sorghum
#cp $base/packages/varsim/tests/quickstart_test_2/insert_seq.txt $base/structural_variant_calling/graph_minigraph-cactus/evaluation2/sandbox/simulated_svs_n_reads/

$base/packages/varsim/varsim.py --reference riouncc.fa --seed 9560 --id simulated \
	--simulator art --read_length 150 --total_coverage 45 --vc_in_vcf sv.vcf \
	--sv_num_ins 100 --sv_num_del 100 --sv_num_dup 0 --sv_num_inv 0 --vc_prop_het 0.1 \
	--simulator_executable $base/packages/varsim/opt/ART/art_bin_VanillaIceCream/art_illumina \
	--sv_prop_het 0.1 --sv_min_length_lim 50 --sv_insert_seq insert_seq.txt \
	--sv_dgv GRCh37_hg19_supportingvariants_2013-07-23.txt

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
