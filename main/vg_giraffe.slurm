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
#
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
#
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
#
#=================================================================================

#: << 'COMMENT'

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

#find $reads_path/*_1.fastq | head -n 5 | parallel -j 5 --joblog $pjl/map_reads_1.txt map_reads {}
#slurm-_____.out

map_reads $reads_path/Grif16309_1.fastq

#COMMENT

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
