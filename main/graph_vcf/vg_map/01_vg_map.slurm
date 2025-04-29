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

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# awk code below was used to fix error reported in slurm-4336820.out
# every string after SVTYPE=DEL; or SVTYPE=INS; in the info column was deleted
#=================================================================================

: << 'COMMENT'

sv_file_path="/projects/cooper_research2/chinaza/structural_variant_calling/graph_vcf/sv_output/final"
awk 'BEGIN {OFS="\t"} !/^#/ {split($8, info, ";"); $8 = info[1]} {print}' $sv_file_path/sv.vcf > input/sv_2.vcf

module load samtools/1.19
bgzip -c input/sv_2.vcf > input/sv_2.vcf.gz
tabix -p vcf input/sv_2.vcf.gz

#slurm-4336828.out

COMMENT

#=================================================================================
# create and index variation graph
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

base_dir="/projects/cooper_research2/chinaza/structural_variant_calling/graph_vcf"

vg autoindex --workflow map \
	--ref-fasta $base_dir/genomes/reference_genome/riouncc.fa \
	--vcf input/sv_2.vcf.gz \
	-T tmp

#slurm-4336831.out

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
	
	vg map -t 60 -x graph_n_index/index.xg \
		-g graph_n_index/index.gcsa \
		-f $1 -f $read_2 \
		> mappings/$read_name.mapped.gam	
}

export -f map_reads
reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"

#find $reads_path/*_1.fastq | head -n 5 | parallel -j 5 --joblog $pjl/map_reads_1.txt map_reads {}
#slurm-4336967.out

#find $reads_path/*_1.fastq | sed -n '6,10p' | parallel -j 5 --joblog $pjl/map_reads_2.txt map_reads {}
#slurm-4336969.out

#find $reads_path/*_1.fastq | tail -n 4 | parallel -j 4 --joblog $pjl/map_reads_3.txt map_reads {}
#slurm-4336971.out

#was still slow, so I ran it individually
#map_reads $reads_path/pi535995_1.fastq
#slurm-4340318.out,4340317,4340316,4340315,4340313,4340312,4340311,4340310,4340309,4340308,4340307,4340306,4340305,4340302

COMMENT

#=================================================================================
# stats of alignment (gam file)
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

find mappings/*gam | parallel -j 14 --joblog $pjl/vg_stats.txt vg_stats {}
#slurm-4593277.out

COMMENT

#=================================================================================
# SV calling and genotyping
#=================================================================================

#: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

# compute the snarls
#vg snarls graph_n_index/index.xg > graph_n_index/index.snarls
# slurm-4593305.out

vg_call() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

        # compute the read support
	# ignore reads with MAPQ < 5 and positions with base quality < 5
        vg pack -x graph_n_index/index.xg \
                -g $1 \
                -o vg_call/$sample.pack \
                -Q 5 \
                -t 10

        vg call graph_n_index/index.xg \
                -r graph_n_index/index.snarls \
                -k vg_call/$sample.pack \
                -s $sample \
                -t 10 \
                > vg_call/$sample.vcf
}

export -f vg_call

find mappings/*gam | parallel -j 14 --joblog $pjl/vg_call.txt vg_call {}
#slurm-4593327.out

#COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'


COMMENT

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
