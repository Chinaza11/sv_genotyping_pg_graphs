#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
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

#=================================================================================
# subset genomes to only one chromosome
#=================================================================================

: << 'COMMENT'

for file in /projects/cooper_research2/chinaza/structural_variant_calling/sandbox/graph_minigraph-cactus_yeast/genomes/*.fa
do
        species_name=$(basename $file | awk -F'.' '{print $1}')
        echo $file
        echo $species_name

        awk '/^>/ {n++} n==1' $file > genomes/$species_name.fa

done

# slurm-4669402.out

COMMENT

#=================================================================================
# pangenome graph construction using minigraph-cactus
#=================================================================================

:<< 'COMMENT'

cactus_dir="/projects/cooper_research2/chinaza/packages/cactus/cactus-bin-v2.9.3/"
source "${cactus_dir}venv-cactus-v2.9.3/bin/activate"
#cactus --version

#toil clean jobstore

cactus-pangenome jobstore \
	genomes/seq_file.txt \
	--outDir yeast-pg \
	--outName yeast-pg \
	--reference S288C \
	--vcf --giraffe --gfa --gbz
# slurm-4670823.out

COMMENT



:<< 'COMMENT'

halStats yeast-pg/yeast-pg.full.hal
echo ""
echo "---------------------------------------------------------------------------"
halStats yeast-pg/yeast-pg.full.hal --sequenceStats S288C

#slurm-___________.out

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

	vg giraffe -Z yeast-pg/yeast-pg.d2.gbz \
                -m yeast-pg/yeast-pg.d2.min \
                -d yeast-pg/yeast-pg.d2.dist \
                --progress \
                -t 30 \
                -f $1 -f $read_2 >mappings/$read_name.mapped.gam
}

export -f map_reads
reads_path="/projects/cooper_research2/chinaza/structural_variant_calling/sandbox/graph_vcf_yeast/vg_autoindex/sequencing_reads"

map_reads $reads_path/SRR4074255_1.fastq
# slurm-4670828.out

COMMENT

#=================================================================================
#
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

vg stats -zNEl yeast-pg/yeast-pg.gbz

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
