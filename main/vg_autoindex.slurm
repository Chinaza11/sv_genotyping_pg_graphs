#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=10-00:00:00
#SBATCH --mem=60GB
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
#
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
#
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
	
	vg map -t 60 -x index.xg -g index.gcsa -f $1 -f $read_2 \
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
#
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

count_mapped_reads() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

	vg view -a $1 | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>0' | wc -l > mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=10' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=20' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=30' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=40' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=50' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' | awk '$2>=60' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' | awk '$2>=1' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' | awk '$2>=0.9' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' | awk '$2>=0.5' | wc -l >> mappings/stats/$sample.tsv
	vg view -a $1 | wc -l >> mappings/stats/$sample.tsv
}

export -f count_mapped_reads

find mappings/*.gam | parallel -j 12 --joblog $pjl/count_mapped_reads.txt count_mapped_reads {}

#slurm-4303947.out

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

for file in mappings/stats/*.tsv
do
	echo ""
        echo "--------------------------------------"
        sample=$(basename $file | awk -F'.' '{print $1}')

        echo $sample
        echo $file
        echo "--------------------------------------"

	cat $file | tr '\n' '\t' >> mappings/stats/all.tsv
	echo -ne "construct\t" >> mappings/stats/all.tsv
	grep $sample mappings/stats/samples.txt | awk '{print $3}' >> mappings/stats/all.tsv

done

#slurm-_______.out

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

sort_gam() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

	vg gamsort -t 10 \
		-i mappings_sorted/$sample.mapped.sorted.gam.gai $1 \
		> mappings_sorted/$sample.mapped.sorted.gam
}

export -f sort_gam

find mappings/*.gam | parallel -j 4 --joblog $pjl/sort_gam.txt sort_gam {}

#slurm-4303978.out

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
