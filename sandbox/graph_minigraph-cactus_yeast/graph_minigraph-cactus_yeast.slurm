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

module load parallel/20240422
pjl="parallel_jobs_log"

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
	--outDir yeast-pg \
	--outName yeast-pg \
	--reference S288C \
	--vcf --giraffe --gfa --gbz
# slurm-4650071.out

halStats yeast-pg/yeast-pg.full.hal
echo ""
echo "---------------------------------------------------------------------------"
halStats yeast-pg/yeast-pg.full.hal --sequenceStats S288C

#slurm-4666215.out

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

	vg giraffe -Z yeast-pg/yeast-pg.d2.gbz \
                -m yeast-pg/yeast-pg.d2.min \
                -d yeast-pg/yeast-pg.d2.dist \
                --progress \
                -t 30 \
                -f $1 -f $read_2 >mappings/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research2/chinaza/structural_variant_calling/sandbox/graph_vcf_yeast/vg_autoindex/sequencing_reads"

# map_reads $reads_path/SRR4074255_1.fastq
# slurm-4666197.out, slurm-4698150.out

map_reads $reads_path/SRR4074256_1.fastq
# slurm-4698158.out ==> vg crashed

COMMENT

#=================================================================================
# graph stats
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

vg stats -zNEl yeast-pg/yeast-pg.gbz

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

        vg giraffe -Z yeast-pg/yeast-pg.d2.gbz \
		-m yeast-pg/yeast-pg.d2.shortread.withzip.min \
		-z yeast-pg/yeast-pg.d2.shortread.zipcodes \
                -d yeast-pg/yeast-pg.d2.dist \
                --progress \
                -t 30 \
                -f $1 -f $read_2 >mappings_2/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research2/chinaza/structural_variant_calling/sandbox/graph_vcf_yeast/vg_autoindex/sequencing_reads"

find $reads_path/*_1.fastq | parallel -j 12 --joblog $pjl/map_reads.txt map_reads {}
# slurm-4698245.out
# Reads SRR4074411, SRR4074412, and SRR4074413 threw a warning.

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

find mappings_2/*gam | parallel -j 12 --joblog $pjl/vg_stats.txt vg_stats {}
# slurm-4698330.out

COMMENT

#=================================================================================
# SV calling and genotyping
#=================================================================================

#: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env2

# compute the snarls
#vg snarls yeast-pg/yeast-pg.d2.gbz > yeast-pg/yeast-pg.d2.snarls
# slurm-4698344.out


vg_call() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"

        # compute the read support
        vg pack -x yeast-pg/yeast-pg.d2.gbz \
                -g $1 \
                -o vg_call/$sample.pack \
                -Q 5 \
                -t 10

        vg call yeast-pg/yeast-pg.d2.gbz \
                -r yeast-pg/yeast-pg.d2.snarls \
                -k vg_call/$sample.pack \
                -s $sample \
                -t 10 \
                > vg_call/$sample.vcf
}

export -f vg_call

find mappings_2/*gam | parallel -j 12 --joblog $pjl/vg_call.txt vg_call {}
#slurm-_________.out

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
