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
# map reads to graph: had to switch to the most recent version of vg (vg 1.65.0)
# code below will recreate minimizer index (that includes zipcodes) and a new zipcodes file
# use these newly created files in mapping other reads
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env2

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
                -t 30 \
                -f $1 -f $read_2 >mappings/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"

map_reads $reads_path/Grif16309_1.fastq
# slurm-4698163.out

COMMENT

#=================================================================================
# map other reads to graph: using newly created zipcodes and minimizer index (with zipcodes) files
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env2

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
		-m sorghum-pg/sorghum-pg.d2.shortread.withzip.min \
                -z sorghum-pg/sorghum-pg.d2.shortread.zipcodes \
                -d sorghum-pg/sorghum-pg.d2.dist \
                --progress \
                -t 30 \
                -f $1 -f $read_2 >mappings/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"

map_reads $reads_path/Grif16309_1.fastq
#master.sh was used to submit multiple slurm jobs
# slurm-4698176.out,4698177,4698178,4698179,4698180,4698181,4698182,4698183,4698184,4698185,4698186,4698187,4698175,4698163

COMMENT

#=================================================================================
# stats of alignments (gam file)
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env2

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
# slurm-___________.out 

COMMENT

#=================================================================================
# SV calling and genotyping
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env2

# compute the snarls
vg snarls sorghum-pg/sorghum-pg.d2.gbz > sorghum-pg/sorghum-pg.d2.snarls
# slurm-___________.out


vg_call() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

        # compute the read support
        vg pack -x sorghum-pg/sorghum-pg.d2.gbz \
                -g $1 \
                -o vg_call/$sample.pack \
                -Q 5 \
                -t 10

        vg call sorghum-pg/sorghum-pg.d2.gbz \
                -r sorghum-pg/sorghum-pg.d2.snarls \
                -k vg_call/$sample.pack \
                -s $sample \
                -t 10 \
                > vg_call/$sample.vcf
}

export -f vg_call

find mappings/*gam | parallel -j 6 --joblog $pjl/vg_call.txt vg_call {}
#slurm-_________.out

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
