#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:00
#SBATCH --mem=20GB
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

#=================================================================================
# minimap2
# https://github.com/lh3/minimap2
# https://lh3.github.io/minimap2/minimap2.html
#=================================================================================

: << 'COMMENT'

process_file() {
	species_name=$(basename $1 | awk -F'.' '{print $1}')
	echo $species_name
}

export -f process_file

find genomes/*.fa | parallel -j 4 --joblog job_log.txt process_file {}

#=================================================================================
#
#=================================================================================

module load anaconda3/2023.09
conda activate minimap2

minimap2 -cx asm5 --cs -t 15 \
	genomes/reference_genome/S288C.genome.fa genomes/SK1.genome.fa > sandbox/SK1.asm5.paf

echo "--------------------------------------"
echo "alignment stats"
paftools.js stat sandbox/SK1.asm5.paf

echo "--------------------------------------"
echo "SK1 asm5 -L"
sort -k6,6 -k8,8n sandbox/SK1.asm5.paf | \
	paftools.js call -L 10000 -f genomes/reference_genome/S288C.genome.fa - > sandbox/SK1_1.vcf	

echo "--------------------------------------"
echo "SK1 asm5 -l"
sort -k6,6 -k8,8n sandbox/SK1.asm5.paf | \
        paftools.js call -l 10000 -f genomes/reference_genome/S288C.genome.fa - > sandbox/SK1_2.vcf

echo "--------------------------------------"
echo "SK1 asm5 -L"
sort -k6,6 -k8,8n sandbox/SK1.asm5.paf | \
        paftools.js call -L 10000 -f genomes/reference_genome/S288C.genome.fa - > sandbox/SK1_1.vcf

echo "--------------------------------------"
echo "SK1 asm5 -l"
sort -k6,6 -k8,8n sandbox/SK1.asm5.paf | \
        paftools.js call -l 10000 -f genomes/reference_genome/S288C.genome.fa - > sandbox/SK1_2.vcf

# slurm-4158388.out

#=================================================================================
#
#=================================================================================

module load anaconda3/2023.09

#conda activate asmvar_dependencies

#conda install matplotlib -y
#conda install bioconda::tabix

#=================================================================================
#
#=================================================================================

conda create -n asmvar_dependencies_2 python=3.9 -y

conda activate asmvar_dependencies_2

conda install -c conda-forge scikit-learn matplotlib tabix -y

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

cd sandbox/

PATH_TO_MUMmer="/projects/cooper_research/Programs/mummer-4.0.0beta2"

$PATH_TO_MUMmer/nucmer -c 100 -b 500 -l 50 S288C.genome.fa SK1.genome.fa
$PATH_TO_MUMmer/delta-filter -m -i 90 -l 100 out.delta > out.filtered.delta
$PATH_TO_MUMmer/show-coords -THrd out.filtered.delta > out.filtered.coords

module load anaconda3
conda activate syri-1.6.3

syri -c out.filtered.coords \
	-d out.filtered.delta \
	-r S288C.genome.fa \
	-q SK1.genome.fa \
	--dir syri_output/ \
	--prefix SK1.

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load bedtools2/2.29.0

# Define input and output directories
INPUT_DIR="sv_output/unions"
OUTPUT_DIR="sandbox/sv_output/unions"

# Define the reference sample
REFERENCE="UFRJ50816"

# Define all other samples dynamically
SAMPLES=("YPS128" "CBS432" "SK1")

# Loop through insertion and deletion categories
for i in ins del; do
    echo "$i"
    echo "--------------------------------------"

    # Start with the reference file
    OUTPUT_FILE="$OUTPUT_DIR/merged.$i.bed"
    cat "$INPUT_DIR/union.$REFERENCE.$i.bed" > "$OUTPUT_FILE.tmp"

    # Iterate over each sample and filter intersections
    for ((idx=0; idx<${#SAMPLES[@]}; idx++)); do
        SAMPLE=${SAMPLES[$idx]}
        bedtools intersect -a "$INPUT_DIR/union.$SAMPLE.$i.bed" -b "$OUTPUT_FILE.tmp" -v -r -f 0.9 >> "$OUTPUT_FILE.tmp"
    done

    # Sort and remove duplicates
    sort -k1,1V -k2,2n -k3,3n -k6,6 -k7,7 -u "$OUTPUT_FILE.tmp" > "$OUTPUT_FILE"
    rm "$OUTPUT_FILE.tmp"
done

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load bedtools2/2.29.0

# Define input and output directories
INPUT_DIR="sv_output/unions"
OUTPUT_DIR="sandbox/2/sv_output/unions"

# Define the reference sample
REFERENCE="UFRJ50816"

# Define all other samples dynamically
SAMPLES=("YPS128" "CBS432" "SK1")

# Loop through insertion and deletion categories
for i in ins del; do
    echo "$i"
    echo "--------------------------------------"

    # Construct the merged file using a scalable nested bedtools approach
    OUTPUT_FILE="$OUTPUT_DIR/merged.$i.bed"
    MERGE_CMD=("cat" "$INPUT_DIR/union.$REFERENCE.$i.bed")
    INTERSECT_B_FILES=("$INPUT_DIR/union.$REFERENCE.$i.bed")

    for SAMPLE in "${SAMPLES[@]}"; do
        MERGE_CMD+=( "<(bedtools intersect -a \"$INPUT_DIR/union.$SAMPLE.$i.bed\" -b ${INTERSECT_B_FILES[*]} -v -r -f 0.9)" )
        INTERSECT_B_FILES+=("$INPUT_DIR/union.$SAMPLE.$i.bed")
    done

    eval "${MERGE_CMD[@]}" | sort -k1,1V -k2,2n -k3,3n -k6,6 -k7,7 -u > "$OUTPUT_FILE"
done

COMMENT

#=================================================================================
#
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
        > sandbox/3/sv_output/unions/merged.$i.bed
done

#=================================================================================
#
#=================================================================================

cat <(
    awk 'OFS="\t" { 
        print $1, $2, $4, $7, substr($7,1,1), $5, "PASS", "SVTYPE=DEL;" $8, "GT", "1/1" 
    }' sandbox/3/sv_output/unions/merged.del.bed
) <(
    awk 'OFS="\t" { 
        print $1, $2, $4, substr($7,1,1), $7, $5, "PASS", "SVTYPE=INS;" $8, "GT", "1/1" 
    }' sandbox/3/sv_output/unions/merged.ins.bed
) | sort -k1,1V -k2,2n -u | cat vcf_header_gt.vcf - > sandbox/3/sv_output/final/sv.vcf

#=================================================================================
#
#=================================================================================

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
