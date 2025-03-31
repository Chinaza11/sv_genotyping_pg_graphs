#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=10-00:00:00
#SBATCH --mem=50GB
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
# download and unzip genomes
#=================================================================================

: << 'COMMENT'

cd genomes/downloaded_genomes

wget --show-progress -i nonsoftmasked_genomes_download_path.txt
gunzip *fa.gz

cd ../..

#=================================================================================
# select only the 10 chromosomes (exclude scaffolds)
# add "chr" to chromosome name to avoid error when running syri
#=================================================================================

for file in genomes/downloaded_genomes/*.fa
do
	species_name=$(basename $file | cut -d'_' -f2 | cut -d'.' -f1)
	echo ""
        echo $file
        echo $species_name
        awk '/^>/ {n++} n<=10' $file > genomes/$species_name.fa
	sed -i 's/>10/>chr10/;
		s/>9/>chr9/;
		s/>8/>chr8/;
		s/>7/>chr7/;
		s/>6/>chr6/;
		s/>5/>chr5/;
		s/>4/>chr4/;
		s/>3/>chr3/;
		s/>2/>chr2/;
		s/>1/>chr1/' genomes/$species_name.fa
done

mv genomes/riouncc.fa genomes/reference_genome/

COMMENT

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
	species_name=$(basename $1 | cut -d'.' -f1)

	echo ""
	echo $species_name
	minimap2 -cx asm20 --cs -t 15 genomes/reference_genome/riouncc.fa $1 > sv_output/minimap2/$species_name.asm20.paf
}

export -f run_minimap

find genomes/*.fa | parallel -j 4 --joblog $pjl/run_minimap.txt run_minimap {}

# slurm-4293722.out

COMMENT

#=================================================================================
# paftools
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate minimap2

#---------------------------------------------------------------------------------
# alignment stats

for file in sv_output/minimap2/*.paf; do
    echo "Processing: $file"
    paftools.js stat "$file"
    echo "--------------------------------------"
done

#---------------------------------------------------------------------------------
# variant calling

paftools_vc() {
	echo ""
        echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')
        echo $species_name
        sort -k6,6 -k8,8n $1 | paftools.js call -L 10000 -f genomes/reference_genome/riouncc.fa - > sv_output/paftools/$species_name.asm20.vcf
        echo "--------------------------------------"
}

export -f paftools_vc

find sv_output/minimap2/*.paf | parallel -j 6 --joblog $pjl/paftools_vc.txt paftools_vc {}

# slurm-4295838.out

COMMENT

#=================================================================================
# covert paftools vcf to bed
#=================================================================================

: << 'COMMENT'

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
		}' > sv_output/bed/paftools.$species_name.del.bed

	# insertions
	grep -ve "^#" $file | awk '
		OFS="\t" { 
			if (length($5)-length($4) > 19 && length($4) < 5) { 
				print $1, $2, $2+(length($5)-1), $3, $6, "Insertion", $5, $8, "'"$species_name"'" 
			} 
		}' > sv_output/bed/paftools.$species_name.ins.bed
done

#slurm-4299585.out

COMMENT


#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# nucmer
#=================================================================================

: << 'COMMENT'

module load anaconda3
conda activate syri-1.6.3

export PATH_TO_MUMmer="/projects/cooper_research/Programs/mummer-4.0.0beta2"

mummer_n_syri() {
	echo ""
	echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')

        echo $species_name
        echo "--------------------------------------"

	$PATH_TO_MUMmer/nucmer -c 100 -b 500 -l 50 \
		genomes/reference_genome/riouncc.fa $1 \
		--prefix sv_output/nucmer/$species_name

	$PATH_TO_MUMmer/delta-filter -m -i 90 -l 100 \
		sv_output/nucmer/$species_name.delta \
		> sv_output/nucmer/$species_name.filtered.delta	

	$PATH_TO_MUMmer/show-coords -THrd \
		sv_output/nucmer/$species_name.filtered.delta \
		> sv_output/nucmer/$species_name.filtered.coords

	syri -c sv_output/nucmer/$species_name.filtered.coords \
		-d sv_output/nucmer/$species_name.filtered.delta \
		-r genomes/reference_genome/riouncc.fa -q $1 \
		--dir sv_output/syri/ \
		--prefix $species_name.
}

export -f mummer_n_syri

find genomes/*.fa | parallel -j 6 --joblog $pjl/mummer_n_syri.txt mummer_n_syri {}

# slurm-4295825.out

COMMENT

#=================================================================================
# covert syri vcf to bed
#=================================================================================

: << 'COMMENT'

for file in sv_output/syri/*.vcf
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
                }' > sv_output/bed/syri.$species_name.del.bed

        # insertions
        grep -ve "^#" $file | awk '
                OFS="\t" {
                        if (length($5)-length($4) > 19 && length($4) < 5) {
                                print $1, $2, $2+(length($5)-1), $3, $6, "Insertion", $5, $8, "'"$species_name"'"
                        }
                }' > sv_output/bed/syri.$species_name.ins.bed
done

#slurm-4300973.out

COMMENT

#=================================================================================
#
#=================================================================================

#=================================================================================
#
#=================================================================================

#=================================================================================
#=================================================================================
#=================================================================================


#=================================================================================
# last
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate last

export dir="genomes/reference_genome"

#lastdb -R01 $dir/lastdb/riouncc $dir/riouncc.fa

run_last() {
        echo ""
        echo "--------------------------------------"
        species_name=$(basename $1 | awk -F'.' '{print $1}')

        echo $species_name
        echo "--------------------------------------"
	
	#lastal -e25 -v -q3 -j4 "$dir/lastdb/riouncc" $1 | last-split -s35 -v | gzip > sv_output/last/$species_name.maf.gz
	lastal -e25 -v -q3 -j4 "$dir/lastdb/riouncc" $1 | last-split -s35 -v | gzip > sv_output/last2/$species_name.maf.gz
}

export -f run_last

#find genomes/*.fa | parallel -j 4 --joblog $pjl/run_last.txt run_last {}
# slurm-4296422.out

# not enough time was selected above so rerunning it for the last 4 genomes to save time
find genomes -type f \( -name "pi5*.fa" -o -name "pi6*.fa" \) | parallel -j 4 --joblog $pjl/run_last2.txt run_last {}
# slurm-4299940.out

COMMENT

#=================================================================================
# asmvar
#=================================================================================

: << 'COMMENT'

run_asmvar() {
        echo ""
        echo "--------------------------------------"
        echo $species_name
	echo $1
	echo chr$2
        echo "--------------------------------------"

	ASV_VariantDetector -s $species_name \
		-r chr$2 \
		-i sv_output/last/$species_name.maf.gz \
		-t genomes/reference_genome/riouncc.fa \
		-q $1 \
		-o sv_output/asmvar/$species_name.chr$2 \
		> sv_output/asmvar/$species_name.chr$2.age \
		2> sv_output/asmvar/$species_name.chr$2.log
}

export -f run_asmvar
export species_name="pi655972" #was changed to each of the genomes name

parallel --joblog $pjl/run_asmvar_$species_name.txt run_asmvar genomes/$species_name.fa ::: {1..10}

#slurm-4327482.out, 4327485, 4327486, 4327488, 4327489, 4327490, 4327492, 4327494, 4327506-4327513

COMMENT

#=================================================================================
# cat asmvar
#=================================================================================

: << 'COMMENT'

for file in genomes/*.fa
do
	echo ""
	echo "--------------------------------------"
	species_name=$(basename $file | awk -F'.' '{print $1}')

	echo $species_name
	echo "--------------------------------------"
	output="sv_output/asmvar"

	cat $(ls -v $output/$species_name*.vcf) | grep -v \"#\" | \
		cat <(grep -e \"^#\" $output/$species_name.chr1.vcf) - > $output/$species_name.allchr.vcf
	
done

#slurm-4333155.out

COMMENT

#=================================================================================
# convert asmvar vcf to bed
#=================================================================================

: << 'COMMENT'

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
		}' > sv_output/bed/asmvar.$species_name.del.bed
	
	# insertions
	grep -ve "^#" $file | awk '
		OFS="\t" { 
			if ($7 == "." || $7 == "AGEFALSE") { 
				if (length($5)-length($4)>19 && length($4) < 5) { 
					print $1, $2, $2+(length($5)-1), $3, $6, "Insertion", $5, $8, "'"$species_name"'" 
				} 
			} 
		}' > sv_output/bed/asmvar.$species_name.ins.bed

done

#slurm-4333204.out

COMMENT

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# merge callsets
# Create high-sensitivity set (union of all three callsets) for each strain
#=================================================================================

: << 'COMMENT'

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
		cat sv_output/bed/paftools.$species_name.$i.bed \
			<(bedtools intersect \
				-a sv_output/bed/asmvar.$species_name.$i.bed \
				-b sv_output/bed/paftools.$species_name.$i.bed \
				-v -r -f 0.5) \
			<(bedtools intersect \
				-a sv_output/bed/syri.$species_name.$i.bed \
				-b sv_output/bed/asmvar.$species_name.$i.bed \
					sv_output/bed/paftools.$species_name.$i.bed \
				-v -r -f 0.5) \
		> sv_output/unions/union.$species_name.$i.bed
	done
done

#slurm_4333218.out

COMMENT

#=================================================================================
# merge strains
#=================================================================================

: << 'COMMENT'

module load bedtools2/2.29.0

INPUT_DIR="sv_output/unions"
OUTPUT_DIR="sv_output/unions"

SAMPLES=("grassl" "grif16309" "leoti" "pi229841" "pi297155" "pi329250" "pi329311" "pi506069" "pi510757" "pi532566" "pi655972")

for i in ins del; do
    echo "$i"
    echo "--------------------------------------"
    OUTPUT_FILE="$OUTPUT_DIR/merged.$i.bed"
    cat "$INPUT_DIR/union.chineseamber.$i.bed" > "$OUTPUT_FILE.tmp"

    for ((idx=0; idx<${#SAMPLES[@]}; idx++)); do
        SAMPLE=${SAMPLES[$idx]}
		echo "$SAMPLE"
        bedtools intersect -a "$INPUT_DIR/union.$SAMPLE.$i.bed" -b "$OUTPUT_FILE.tmp" -v -r -f 0.9 >> "$OUTPUT_FILE.tmp"
    done

    sort -k1,1V -k2,2n -k3,3n -k6,6 -k7,7 -u "$OUTPUT_FILE.tmp" > "$OUTPUT_FILE"
    rm "$OUTPUT_FILE.tmp"

done

# slurm-4334433.out

COMMENT

#=================================================================================
# convert bed to vcf
#=================================================================================

: << 'COMMENT'

cat <(
    awk 'OFS="\t" { 
        print $1, $2, $4, $7, substr($7,1,1), $5, "PASS", "SVTYPE=DEL;" $8, "GT", "1/1" 
    }' sv_output/unions/merged.del.bed
) <(
    awk 'OFS="\t" { 
        print $1, $2, $4, substr($7,1,1), $7, $5, "PASS", "SVTYPE=INS;" $8, "GT", "1/1" 
    }' sv_output/unions/merged.ins.bed
) | sort -k1,1V -k2,2n -u | cat vcf_header_gt.vcf - > sv_output/final/sv.vcf

# slurm-4334472.out

COMMENT

#=================================================================================
# compress sv.vcf and index it
#=================================================================================

: << 'COMMENT'

module load samtools/1.19

bgzip -c sv_output/final/sv.vcf > sv_output/final/sv.vcf.gz

tabix -p vcf sv_output/final/sv.vcf.gz 

# slurm-4334473.out

COMMENT

#=================================================================================
#
#=================================================================================


echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
