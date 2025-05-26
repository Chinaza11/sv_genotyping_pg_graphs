#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --time=10-00:00:00
#SBATCH --mem=700GB
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
# activate conda environment
#=================================================================================

module load anaconda3/2023.09
conda activate vg_env2

#=================================================================================
# map reads to graph: had to switch to the most recent version of vg (vg 1.65.0)
# code below will recreate minimizer index (that includes zipcodes) and a new zipcodes file
# use these newly created files in mapping other reads
#=================================================================================

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

#reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"
#map_reads $reads_path/Grif16309_1.fastq
# slurm-4698163.out

#=================================================================================
# function for mapping other reads to graph:
#	using the newly created zipcodes and minimizer index (with zipcodes) files
#=================================================================================

map_reads2() {
        echo ""
        echo "--------------------------------------"
        read_name=$(basename $1 | awk -F'_' '{print $1}')

	if [[ "$1" == *_1.fastq || "$1" == *_1.fastq.gz ]]; then
		read_name=$(basename $1 | awk -F'_' '{print $1}')
	        read_2=$(echo $1 | sed 's/_1/_2/')
	elif [[ "$1" == *_R1.fastq || "$1" == *_R1.fastq.gz ]]; then
		read_name=$(basename $1 | awk -F'.' '{print $1}' | sed 's/_R1//')
		read_2=$(echo "$1" | sed 's/_R1/_R2/')
	else
		echo "Unrecognized read 1 format: $1" >&2
		return 1
	fi

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
                -f $1 -f $read_2 >$2/$read_name.mapped.gam
}

export -f map_reads2

#=================================================================================
# map other Wild_Sorghum_WGS reads to graph using map_read2 function
#=================================================================================

#reads_path="/projects/cooper_research1/Wild_Sorghum_WGS"
#map_reads2 $reads_path/Grif16309_1.fastq "mappings_wild_sorghum_wgs"

#master.sh was used to submit multiple slurm jobs
# slurm-4698176.out,4698177,4698178,4698179,4698180,4698181,4698182,4698183,4698184,4698185,4698186,4698187,4698175,4698163

#=================================================================================
# stats of alignments (gam file) function
#=================================================================================

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

#=================================================================================
# stats of alignments (gam file): Wild_Sorghum_WGS
#=================================================================================

#find mappings_wild_sorghum_wgs/*gam \
#	| parallel -j 14 --joblog $pjl/vg_stats_wild_sorghum_wgs.txt vg_stats {}
# slurm-4699086.out

#=================================================================================
# SV calling and genotyping function
#=================================================================================

# compute the snarls
#vg snarls sorghum-pg/sorghum-pg.d2.gbz > sorghum-pg/sorghum-pg.d2.snarls
# slurm-4699090.out

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
                -o $2/$sample.pack \
                -Q 5 \
                -t 20

        vg call sorghum-pg/sorghum-pg.d2.gbz \
                -r sorghum-pg/sorghum-pg.d2.snarls \
                -k $2/$sample.pack \
                -s $sample \
                -t 20 \
                > $2/$sample.vcf
}

export -f vg_call

#=================================================================================
# SV calling and genotyping: Wild_Sorghum_WGS
#=================================================================================

#find mappings_wild_sorghum_wgs/*gam \
#	| parallel -j 14 --joblog $pjl/vg_call_wild_sorghum_wgs.txt vg_call {} "vg_call_wild_sorghum_wgs"

#slurm-4699097.out

#=================================================================================
#=================================================================================
#=================================================================================

#=================================================================================
# map TERRA_RAW reads to graph using map_read2 function
#=================================================================================

#reads_path="/projects/cooper_research1/TERRA_RAW"
#map_reads2 $reads_path/placeholder.fastq.gz "mappings_terra_raw"

#master.sh was used to submit multiple slurm jobs
# slurm_log/mapping_terra_raw/

#=================================================================================
# investigate the warnings in the terra_raw mapping files
#=================================================================================

search_logs() {
    local pattern="$1"
    local label="$2"

    echo "===== Searching for: $label ====="
    match_count=0
    total_files=0

    for file in slurm_log/mapping_terra_raw/*; do
        ((total_files++))
        if grep -qi "$pattern" "$file"; then
            echo "------------------------------------------"
            echo "$file"
            grep -i "$pattern" "$file"
            ((match_count++))
        fi
    done

    echo -e "\n$match_count/$total_files files contain: \"$label\"\n"
}

#search_logs "warning" "warning"
# slurm-4719526.out ==> 196/383 files have warnings

#search_logs "warning\[vg::Watchdog\]" "warning[vg::Watchdog]"
# slurm-4719528.out ==> 40/383 files have warnings

#search_logs "Encountered 100000 ambiguously-paired reads before finding enough" "ambiguous paired reads warning"
# slurm-4719533.out ==> 12/383 files have warnings

#search_logs "Refusing to perform too-large rescue alignment" "too-large rescue alignment warning"
# slurm-4719535.out ==> 165/383 files have warnings

#=================================================================================
# stats of alignments (gam file): TERRA_RAW
# uncomment vg_stats function above
#=================================================================================

#find mappings_terra_raw/*gam | parallel -j 20 --joblog $pjl/vg_stats_terra_raw.txt vg_stats {}
# slurm-4719547.out

#=================================================================================
# SV calling and genotyping: TERRA_RAW
#=================================================================================

#find mappings_terra_raw/*gam | parallel -j 20 --joblog $pjl/vg_call_terra_raw.txt vg_call {} "vg_call_terra_raw"
#slurm-_______.out

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

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"
