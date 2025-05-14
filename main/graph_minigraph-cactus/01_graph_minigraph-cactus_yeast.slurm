#! /bin/bash

#SBATCH --partition=Draco
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

#=================================================================================
# download and unzip genomes: download interactively in the DTN node
#=================================================================================

: << 'COMMENT'

cd genomes

wget --show-progress -i nonsoftmasked_genomes_download_path.txt
gunzip *fa.gz

COMMENT

#=================================================================================
# select only the 10 chromosomes, exclude scaffolds
# add "chr" to chromosomes name
#=================================================================================

: << 'COMMENT'

for file in genomes/*.fa
do
	species_name=$(basename $file | awk -F'.' '{print $1}' | awk -F'_' '{print $2}')
        echo $file
        echo $species_name

        awk '/^>/ {n++} n<=10' $file > genomes/$species_name.fa

	awk '
		/^>/ {
		sub(/^>1$/, ">chr1")
		sub(/^>2$/, ">chr2")
		sub(/^>3$/, ">chr3")
		sub(/^>4$/, ">chr4")
		sub(/^>5$/, ">chr5")
		sub(/^>6$/, ">chr6")
		sub(/^>7$/, ">chr7")
		sub(/^>8$/, ">chr8")
		sub(/^>9$/, ">chr9")
		sub(/^>10$/, ">chr10")
	}
	{ print }
	' genomes/$species_name.fa > genomes/$species_name.tmp
	
	mv genomes/$species_name.tmp genomes/$species_name.fa
done

rm genomes/*.toplevel.fa

COMMENT

#=================================================================================
# pangenome graph construction using minigraph-cactus
#=================================================================================

: << 'COMMENT'

cactus_dir="/projects/cooper_research2/chinaza/packages/cactus/cactus-bin-v2.9.3/"
source "${cactus_dir}venv-cactus-v2.9.3/bin/activate"
#cactus --version

toil clean jobstore

cactus-pangenome jobstore \
	genomes/seq_file.txt \
	--outDir sorghum-pg \
	--outName sorghum-pg \
	--reference riouncc \
	--vcf --giraffe --gfa --gbz
# slurm-4651095.out


halStats sorghum-pg/sorghum-pg.full.hal
# slurm-4666212.out

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

        vg giraffe -Z sorghum-pg/sorghum-pg.d2.gbz \
                -m sorghum-pg/sorghum-pg.d2.min \
                -d sorghum-pg/sorghum-pg.d2.dist \
                --progress \
		--fragment-mean 258.681 \
		--fragment-stdev 104.517 \
                -t 30 \
                -f $1 -f $read_2 >mappings_2/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"

map_reads $reads_path/Grif16309_1.fastq
#master.sh was used submit multiple slurm jobs
#slurm-4659495.out,4659517,4659518,4659519,4659520,4659521,4659516,4659515,4659514,4659513,4659512,4659511,4659509,4659510

#slurm-4659545.out #specified fragment length and stdev

COMMENT

#=================================================================================
#
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

vg_stats mappings_2/Grif16309.mapped.gam


#find mappings/*gam | parallel -j 7 --joblog $pjl/vg_stats.txt vg_stats {}
# slurm-4659544.out, 

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
