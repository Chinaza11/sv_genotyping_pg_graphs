#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10-00:00:00
#SBATCH --mem=300GB
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
#=================================================================================
#=================================================================================

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load samtools/1.19

for i in vg_call/*vcf; do
	echo ""
	echo "-----------------------------------"
	sample=$(basename $i | awk -F'.' '{print $1}')
	echo $i
	echo $sample
	echo "-----------------------------------"

	bgzip -c $i > vg_call/$sample.vcf.gz
	tabix -p vcf vg_call/$sample.vcf.gz

	bcftools filter \
		-e 'GT="0" | GT="0/0" | GT~"\." | ABS(STRLEN(REF) - STRLEN(ALT)) < 50' $i.gz \
		| sed 's/^S288C.//g' | bgzip -c > vg_call/$sample.filtered.vcf.gz

	tabix -p vcf vg_call/$sample.filtered.vcf.gz
done

COMMENT

#=================================================================================
# make linear graph and index it
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

vg autoindex --workflow map \
	--prefix evaluation/graphs/empty_sample_graph/linear_graph \
        --ref-fasta input/S288C.genome.fa \
        -T tmp

#slurm-4472970.out

COMMENT

#=================================================================================
# construct sample graphs and index them
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

make_sample_graph() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"
	echo "Running vg construct"
	echo ""
	vg construct -r input/S288C.genome.fa -v $1 \
		-S -a -f -p \
		> evaluation/graphs/sample_graph/$sample.reconstructed.vg

	echo ""
        echo "--------------------------------------"
        echo "Running vg mod"
        echo ""
	vg mod -v $1 evaluation/graphs/sample_graph/$sample.reconstructed.vg \
		> evaluation/graphs/sample_graph/$sample.vg

	echo ""
        echo "--------------------------------------"
        echo "Running vg convert"
        echo ""
	vg convert -f evaluation/graphs/sample_graph/$sample.vg \
		> evaluation/graphs/sample_graph/$sample.gfa

        echo ""
        echo "--------------------------------------"
        echo "Running vg autoindex"
        echo ""
	vg autoindex --workflow map \
		--prefix evaluation/graphs/sample_graph/$sample \
		--gfa evaluation/graphs/sample_graph/$sample.gfa \
        	-T tmp
}

export -f make_sample_graph

find vg_call/*filtered.vcf.gz \
	| parallel -j 12 --joblog $pjl/make_sample_graph.txt make_sample_graph {}

#slurm-4473361.out

COMMENT

#=================================================================================
#
#=================================================================================


#=================================================================================
#
#=================================================================================


#=================================================================================
#
#=================================================================================


#=================================================================================
#
#=================================================================================

#=================================================================================
#
#=================================================================================

#=================================================================================
#
#=================================================================================

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
