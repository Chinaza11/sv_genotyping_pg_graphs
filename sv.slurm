#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=22
#SBATCH --time=3-00:00:00
#SBATCH --mem=1000GB
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

#: << 'COMMENT'

module load anaconda3/2023.09
conda activate minimap2

#---------------------------------------------------------------------------------
# alignment stats

for file in sv_output/minimap2/*.paf; do
    echo "Processing: $file"
    paftools.js stat "$file"
    echo "--------------------------------------"
done

# slurm-4159729.out
# slurm-_______.out

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

# slurm-4159743.out
# slurm-_______.out

#COMMENT

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

# slurm-_____.out

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
