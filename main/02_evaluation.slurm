#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10-00:00:00
#SBATCH --mem=500GB
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

module load anaconda3/2023.09
conda activate vg_env

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# create and index linear reference genome
#=================================================================================

: << 'COMMENT'

base_dir="/projects/cooper_research2/chinaza/structural_variant_calling/graph_vcf"

vg autoindex --workflow giraffe \
        --prefix evaluation/graphs/linear_genome/linear_genome \
        --ref-fasta $base_dir/genomes/reference_genome/riouncc.fa \
        -T tmp

#slurm-4561637.out

COMMENT

#=================================================================================
# map reads to graph version of the linear reference genome
#=================================================================================

: << 'COMMENT'

map_reads() {
        echo ""
        echo "--------------------------------------"
        read_name=$(basename $1 | awk -F'_' '{print $1}')
        read_2=$(echo $1 | sed 's/_1/_2/')

        echo $read_name
        echo $1
        echo $read_2
        echo "--------------------------------------"

        vg giraffe -Z $2.giraffe.gbz \
                -m $2.shortread.withzip.min \
                -z $2.shortread.zipcodes \
                -d $2.dist \
                --progress \
                -t 30 \
                -f $1 -f $read_2 > $3/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"

map_reads $reads_path/Grif16309_1.fastq \
	"evaluation/graphs/linear_genome/linear_genome" \
	"evaluation/mappings/linear_genome"
 
#master.sh was used to submit multiple slurm jobs
#slurm-4561655.out,4561656,4561657,4561658,4561659,4561660,4561661,4561662,4561663,4561664,4561665,4561666,4561667,4561669

COMMENT

#=================================================================================
# stats of alignments to the linear reference genome
#=================================================================================

: << 'COMMENT'

vg_stats() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

        vg stats -a $1
}

export -f vg_stats

find evaluation/mappings/linear_genome/*gam | parallel -j 7 --joblog $pjl/vg_stats_linear_genome.txt vg_stats {}
#slurm-4561860.out

COMMENT

#=================================================================================
#
#=================================================================================

#: << 'COMMENT'

# compute the snarls
#vg snarls evaluation/graphs/linear_genome/linear_genome.giraffe.gbz \
#	> evaluation/graphs/linear_genome/linear_genome.snarls
# slurm-4561865.out

vg_call() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

        # compute the read support
        vg pack -x $2.giraffe.gbz \
                -g $1 \
                -o $3/$sample.pack \
                -Q 5 \
                -t 10

        vg call $2.giraffe.gbz \
		-r $2.snarls \
                -k $3/$sample.pack \
                -s $sample \
                -t 10 \
                > $3/$sample.vcf
}

export -f vg_call

find evaluation/mappings/linear_genome/*gam \
	| parallel -j 7 --joblog $pjl/vg_call_linear_genome.txt \
		vg_call {} \
			evaluation/graphs/linear_genome/linear_genome \
			evaluation/calls/linear_genome
#slurm-________.out

#COMMENT

#=================================================================================
#
#=================================================================================

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
