#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10-00:00:00
#SBATCH --mem=200GB
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
# phase vcf, zip it and index it
#=================================================================================

: << 'COMMENT'

path="/projects/cooper_research2/chinaza/structural_variant_calling/graph_vcf/vg_autoindex"
cp $path/input/sv_2.vcf input/

#make phased version of file
sed 's#1/1#1|1#g' input/sv_2.vcf > input/sv_2_phased.vcf
rm input/sv_2.vcf

module load samtools/1.19
bgzip -c input/sv_2_phased.vcf > input/sv_2_phased.vcf.gz
tabix -p vcf input/sv_2_phased.vcf.gz

#slurm-4340403.out

COMMENT

#=================================================================================
# create and index variation graph
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

base_dir="/projects/cooper_research2/chinaza/structural_variant_calling/graph_vcf"

vg autoindex --workflow giraffe \
	--prefix graph_n_index/giraffe \
	--ref-fasta $base_dir/genomes/reference_genome/riouncc.fa \
	--vcf input/sv_2_phased.vcf.gz \
	-T tmp

#slurm-4340554.out

COMMENT

#=================================================================================
# map reads to graph
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

map_reads() {
        echo ""
        echo "--------------------------------------"
        read_name=$(basename $1 | awk -F'_' '{print $1}')
        read_2=$(echo $1 | sed 's/_1/_2/')

        echo $read_name
        echo $1
        echo $read_2
        echo "--------------------------------------"

	vg giraffe -Z graph_n_index/giraffe.giraffe.gbz \
		-m graph_n_index/giraffe.shortread.withzip.min \
		-z graph_n_index/giraffe.shortread.zipcodes \
		-d graph_n_index/giraffe.dist \
		--progress \
		-t 30 \
		-f $1 -f $read_2 >mappings/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"

#find $reads_path/*_1.fastq | parallel -j 5 --joblog $pjl/map_reads.txt map_reads {}
#slurm-_____.out

map_reads $reads_path/Grif16309_1.fastq 
#master.sh was used submit multiple slurm jobs
#slurm-4340607.out,4449637,4449638,4449639,4452105,4452106,4452107,4452108,4452767,4452768,4452769,4453229,4453230,4453231

COMMENT

#=================================================================================
# stats of alignments (gam file)
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

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

find mappings/*gam | parallel -j 7 --joblog $pjl/vg_stats.txt vg_stats {}
#slurm-4561631.out

COMMENT

#=================================================================================
# SV calling and genotyping
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

# compute the snarls
#vg snarls graph_n_index/giraffe.giraffe.gbz > graph_n_index/giraffe.snarls
# slurm-4456489.out


vg_call() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

        # compute the read support
        vg pack -x graph_n_index/giraffe.giraffe.gbz \
		-g $1 \
		-o vg_call/$sample.pack \
		-Q 5 \
		-t 10

        vg call graph_n_index/giraffe.giraffe.gbz \
		-r graph_n_index/giraffe.snarls \
		-k vg_call/$sample.pack \
		-s $sample \
		-t 10 \
                > vg_call/$sample.vcf
}

export -f vg_call

find mappings/*gam | parallel -j 6 --joblog $pjl/vg_call.txt vg_call {}
#slurm-4456494.out

COMMENT

#=================================================================================
#
#=================================================================================


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
