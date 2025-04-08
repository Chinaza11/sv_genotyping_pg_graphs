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
# awk code below was used to fix error reported in slurm-4195389.out
# every string after SVTYPE=DEL; or SVTYPE=INS; in the info column was deleted
#=================================================================================

: << 'COMMENT'

#awk 'BEGIN {OFS="\t"} !/^#/ {split($8, info, ";"); $8 = info[1]} {print}' input/sv.vcf > input/sv_2.vcf

#module load samtools/1.19
#bgzip -c input/sv_2.vcf > input/sv_2.vcf.gz
#tabix -p vcf input/sv_2.vcf.gz

#slurm-4195543.out

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

vg autoindex --workflow map \
	--ref-fasta input/S288C.genome.fa \
	--vcf input/sv_2.vcf.gz \
	-T tmp

#slurm-4195548.out

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load sra-tools/3.1.0

cd sequencing_reads
cat runids.txt | parallel -j 12 --tag --joblog job_log.txt fasterq-dump {}

#slurm-4200233.out

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
	
	vg map -t 10 -x index.xg -g index.gcsa -f $1 -f $read_2 \
		> mappings/$read_name.mapped.gam	
}

export -f map_reads

find sequencing_reads/*_1.fastq | parallel -j 12 --joblog $pjl/map_reads.txt map_reads {}

#slurm-4200261.out

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

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

# compute the snarls
#vg snarls index.xg > index.snarls
# slurm-4339452.out

vg_call() {
	echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

	# compute the read support
	vg pack -x index.xg -g $1 -o vg_call/$sample.pack -Q 5 -t 10

	vg call index.xg -r index.snarls -k vg_call/$sample.pack -s $sample -t 10 \
		> vg_call/$sample.vcf
}

export -f vg_call

find mappings/*gam | parallel -j 12 --joblog $pjl/vg_call.txt vg_call {}
#slurm-4339481.out

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

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
