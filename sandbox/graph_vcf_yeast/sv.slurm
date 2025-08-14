#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=20
#SBATCH --time=3-00:00:00
#SBATCH --mem=100GB
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
pjl="parallel_jobs_log

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# minimap2
# https://github.com/lh3/minimap2
# https://lh3.github.io/minimap2/minimap2.html
#=================================================================================

: << 'COMMENT'

module load minimap2/2.28

run_minimap() {
	species_name=$(basename $1 | awk -F'.' '{print $1}')
	minimap2 -cx asm20 --cs -t 15 genomes/reference_genome/S288C.genome.fa $1 > sv_output/minimap2/$species_name.asm20.paf
}

export -f run_minimap

find genomes/*.fa | parallel -j 4 --joblog $pjl/run_minimap.txt run_minimap {}

#=================================================================================
# paftools
#=================================================================================

module load anaconda3/2023.09
conda activate minimap2

#---------------------------------------------------------------------------------
# alignment stats

for file in sv_output/minimap2/*.paf; do
    echo "Processing: $file"
    paftools.js stat "$file"
    echo "--------------------------------------"
done

# slurm-4158369.out

#---------------------------------------------------------------------------------
# variant calling

paftools_vc() {
	echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')
	echo $species_name
	sort -k6,6 -k8,8n $1 | paftools.js call -L 10000 -f genomes/reference_genome/S288C.genome.fa - > sv_output/paftools/$species_name.asm20.vcf
	echo "--------------------------------------"
}

export -f paftools_vc

find sv_output/minimap2/*.paf | parallel -j 4 --joblog $pjl/paftools_vc.txt paftools_vc {}

# slurm-4158384.out

#=================================================================================
# covert paftools vcf to bed
#=================================================================================

for file in sv_output/paftools/*.vcf
do
	echo ""
	echo "--------------------------------------"
	species_name=$(basename $file | awk -F'.' '{print $1}')
	echo $species_name
	echo "--------------------------------------"
	
	# deletions
	grep -ve "^#" $file | awk '
		OFS="\t" { 
			if (length($4)-length($5) > 19 && length($5) < 5) { 
				print $1, $2, $2+(length($4)-1), $3, $6, "Deletion", $4, $8, "'"$species_name"'" 
			} 
		}' > sv_output/paftools.$species_name.del.bed

	# insertions
	grep -ve "^#" $file | awk '
		OFS="\t" { 
			if (length($5)-length($4) > 19 && length($4) < 5) { 
				print $1, $2, $2+(length($5)-1), $3, $6, "Insertion", $5, $8, "'"$species_name"'" 
			} 
		}' > sv_output/paftools.$species_name.ins.bed
done

COMMENT

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# nucmer
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate mummer3

run_nucmer() {
	echo ""
        echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')

        echo $species_name
	echo "--------------------------------------"
	
	# alignment
	nucmer --maxmatch --minmatch 15 --maxgap 500 --mincluster 200 --diagdiff 20 genomes/reference_genome/S288C.genome.fa $1 -prefix sv_output/nucmer/$species_name
	
	# alignment stats
	dnadiff -prefix sv_output/nucmer/alignment_stats/$species_name -d sv_output/nucmer/$species_name.delta
}

export -f run_nucmer

find genomes/*.fa | parallel -j 4 --joblog $pjl/run_nucmer.txt run_nucmer {}

# slurm-4161201.out

#=================================================================================
# assemblytics
#=================================================================================

module load R/4.3.3

run_assemblytics() {
        echo ""
        echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')

        echo $species_name
        echo "--------------------------------------"

	Assemblytics $1 sv_output/assemblytics/$species_name 5000 50 10000
}

export -f run_assemblytics

find sv_output/nucmer/*.delta | parallel -j 4 --joblog $pjl/run_assemblytics.txt run_assemblytics {}

# slurm-4161631.out

#=================================================================================
# convert assemblytics bed to bed
#=================================================================================

for file in sv_output/assemblytics/*.Assemblytics_structural_variants.bed
do
        echo ""
        echo "--------------------------------------"
        species_name=$(basename $file | awk -F'.' '{print $1}')
        echo $species_name
        echo "--------------------------------------"

	awk 'OFS="\t" {
		if ($7 == "Insertion" && $8 > -5 && $8 < 5 && $9 > 19) {
			print $1, $2, $2+$9, $4, ".", $7, $10, ".", "'"$species_name"'"
		} 
		else if ($7 == "Deletion" && $9 >= 0 && $9 < 5 && $8 > 19) {
			print $1, $2, $2+$8, $4, ".", $7, $1":"$2"-"$3, ".", "'"$species_name"'"
		}
	}' $file | sort -k1,1 -k2,2n -k3,3n -k6,6 -k7,7 -u > sv_output/assemblytics/$species_name.coords.bed

done	

#=================================================================================
# replace assemblytics *coords.bed file with nucleotides
# ====>> comeback later doublecheck the code in replace_regions_with_nucleotides.py <<====
#=================================================================================

module load bedtools2/2.29.0

for file in sv_output/assemblytics/*.coords.bed
do
	echo ""
        echo "--------------------------------------"
        species_name=$(basename $file | awk -F'.' '{print $1}')
        echo $species_name
        echo "--------------------------------------"

	python3 replace_regions_with_nucleotides.py \
		$file \
		genomes/reference_genome/S288C.genome.fa \
		genomes/$species_name.genome.fa \
		$species_name > sv_output/assemblytics/$species_name.final.bed	
done

#=================================================================================
# put assemblytics deletions and insertions in separate files
#=================================================================================

for file in sv_output/assemblytics/*.final.bed
do
	echo ""
        echo "--------------------------------------"
        species_name=$(basename $file | awk -F'.' '{print $1}')
        echo $species_name
        echo "--------------------------------------"

	grep Deletion $file > sv_output/assemblytics.$species_name.del.bed
	grep Insertion $file > sv_output/assemblytics.$species_name.ins.bed

done

COMMENT

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# last
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09

conda activate last

dir="genomes/reference_genome"

#lastdb -R01 $dir/lastdb/S288C $dir/S288C.genome.fa

run_last() {
        echo ""
        echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')

        echo $species_name
        echo "--------------------------------------"
	
	dir="genomes/reference_genome"

	lastal -e25 -v -q3 -j4 "$dir/lastdb/S288C" $1 | last-split -s35 -v | gzip > sv_output/last/$species_name.maf.gz
}

export -f process_file

find genomes/*.fa | parallel -j 4 --joblog $pjl/run_last.txt run_last {}

# slurm-4166903.out

#=================================================================================
# asmvar
#=================================================================================

run_asmvar() {
        echo ""
        echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')

        echo $species_name
        echo "--------------------------------------"

	for chr in $(grep ">" genomes/SK1.genome.fa | sed 's/>//')
	do
		echo "--------------------------------------"
		echo $chr
		echo "--------------------------------------"
		ASV_VariantDetector -s $species_name \
			-r $chr \
			-i sv_output/last/$species_name.maf.gz \
			-t genomes/reference_genome/S288C.genome.fa \
			-q $1 \
			-o sv_output/asmvar/$species_name.$chr \
			> sv_output/asmvar/$species_name.$chr.age \
			2> sv_output/asmvar/$species_name.$chr.log
	done
}

export -f run_asmvar

find genomes/*.fa | parallel -j 4 --joblog $pjl/run_asmvar.txt run_asmvar {}

# slurm-4172225.out

#=================================================================================
# cat asmvar
#=================================================================================

for file in genomes/*.fa
do
	echo ""
	echo "--------------------------------------"
	species_name=$(basename $file | awk -F'.' '{print $1}')

	echo $species_name
	echo "--------------------------------------"
	output="sv_output/asmvar"

	cat $output/$species_name*.vcf | grep -v \"#\" | \
		cat <(grep -e \"^#\" $output/$species_name.chrI.vcf) - > $output/$species_name.allchr.vcf
	
done

#=================================================================================
# convert asmvar vcf to bed
#=================================================================================

for file in sv_output/asmvar/*allchr.vcf
do
	echo ""
        echo "--------------------------------------"
        species_name=$(basename $file | awk -F'.' '{print $1}')
        echo $species_name
        echo "--------------------------------------"
	
	# deletions
	grep -ve "^#" $file | awk '
		OFS="\t" {
			if ($7 == "." || $7 == "AGEFALSE") { 
				if (length($4)-length($5)>19 && length($5) < 5) { 
					print $1, $2, $2+(length($4)-1), $3, $6, "Deletion", $4, $8, "'"$species_name"'" 
				} 
			} 
		}' > sv_output/asmvar.$species_name.del.bed
	
	# insertions
	grep -ve "^#" $file | awk '
		OFS="\t" { 
			if ($7 == "." || $7 == "AGEFALSE") { 
				if (length($5)-length($4)>19 && length($4) < 5) { 
					print $1, $2, $2+(length($5)-1), $3, $6, "Insertion", $5, $8, "'"$species_name"'" 
				} 
			} 
		}' > sv_output/asmvar.$species_name.ins.bed

done

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# merge callsets
# Create high-sensitivity set (union of all three callsets) for each strain
#=================================================================================

module load bedtools2/2.29.0

for file in genomes/*.fa
do
        echo ""
        echo "--------------------------------------"
        species_name=$(basename $file | awk -F'.' '{print $1}')
        echo $species_name
        echo "--------------------------------------"

	for i in ins del;
	do
		echo $i
		echo "--------------------------------------"
		cat sv_output/paftools.$species_name.$i.bed \
			<(bedtools intersect \
				-a sv_output/asmvar.$species_name.$i.bed \
				-b sv_output/paftools.$species_name.$i.bed \
				-v -r -f 0.5) \
			<(bedtools intersect \
				-a sv_output/assemblytics.$species_name.$i.bed \
				-b sv_output/asmvar.$species_name.$i.bed \
					sv_output/paftools.$species_name.$i.bed \
				-v -r -f 0.5) \
		> sv_output/unions/union.$species_name.$i.bed
	done
done

#=================================================================================
# merge strains
#=================================================================================

module load bedtools2/2.29.0

for i in ins del;
do
	echo $i
	echo "--------------------------------------"
	cat sv_output/unions/union.UFRJ50816.$i.bed \
		<(bedtools intersect \
			-a sv_output/unions/union.YPS128.$i.bed \
			-b sv_output/unions/union.UFRJ50816.$i.bed \ 
			-v -r -f 0.9) \
		<(bedtools intersect \
			-a sv_output/unions/union.CBS432.$i.bed \
			-b sv_output/unions/union.UFRJ50816.$i.bed \
				sv_output/unions/union.YPS128.$i.bed \
			-v -r -f 0.9) \
		<(bedtools intersect \
			-a sv_output/unions/union.SK1.$i.bed \
			-b sv_output/unions/union.UFRJ50816.$i.bed \
				sv_output/unions/union.YPS128.$i.bed \
				sv_output/unions/union.CBS432.$i.bed \
			-v -r -f 0.9) \
	| sort -k1,1V -k2,2n -k3,3n -k6,6 -k7,7 -u \
	> sv_output/unions/merged.$i.bed
done

#=================================================================================
# convert bed to vcf
#=================================================================================

cat <(
    awk 'OFS="\t" { 
        print $1, $2, $4, $7, substr($7,1,1), $5, "PASS", "SVTYPE=DEL;" $8, "GT", "1/1" 
    }' sv_output/unions/merged.del.bed
) <(
    awk 'OFS="\t" { 
        print $1, $2, $4, substr($7,1,1), $7, $5, "PASS", "SVTYPE=INS;" $8, "GT", "1/1" 
    }' sv_output/unions/merged.ins.bed
) | sort -k1,1V -k2,2n -u | cat vcf_header_gt.vcf - > sv_output/final/sv.vcf

#=================================================================================
# compress sv.vcf and index it
#=================================================================================

module load samtools/1.19

bgzip -c sv_output/final/sv.vcf > sv_output/final/sv.vcf.gz

tabix -p vcf sv_output/final/sv.vcf.gz 

COMMENT

#=================================================================================
#
#=================================================================================



echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
