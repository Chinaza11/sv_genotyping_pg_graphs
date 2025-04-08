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
#
#=================================================================================

: << 'COMMENT'

module load samtools/1.19

for i in vg_call/*vcf; do
	echo ""
	echo "-----------------------------------"
	sample=$(basename $i | awk -F'.' '{print $1}')
	echo $i
	echo $sample
	echo "-----------------------------------"

	bgzip -c $i > vg_call/$sample.vcf.gz
	tabix -p vcf vg_call/$sample.vcf.gz

	bcftools filter \
		-e 'GT="0" | GT="0/0" | GT~"\." | ABS(STRLEN(REF) - STRLEN(ALT)) < 50' $i.gz \
		| sed 's/^S288C.//g' | bgzip -c > vg_call/$sample.filtered.vcf.gz

	tabix -p vcf vg_call/$sample.filtered.vcf.gz
done

COMMENT

#=================================================================================
# make linear graph and index it
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

vg autoindex --workflow map \
	--prefix evaluation/graphs/empty_sample_graph/linear_graph \
        --ref-fasta input/S288C.genome.fa \
        -T tmp

#slurm-4472970.out

COMMENT

#=================================================================================
# construct sample graphs and index them
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

make_sample_graph() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
        echo "--------------------------------------"
	echo "Running vg construct"
	echo ""
	vg construct -r input/S288C.genome.fa -v $1 \
		-S -a -f -p \
		> evaluation/graphs/sample_graph/$sample.reconstructed.vg

	echo ""
        echo "--------------------------------------"
        echo "Running vg mod"
        echo ""
	vg mod -v $1 evaluation/graphs/sample_graph/$sample.reconstructed.vg \
		> evaluation/graphs/sample_graph/$sample.vg

	echo ""
        echo "--------------------------------------"
        echo "Running vg convert"
        echo ""
	vg convert -f evaluation/graphs/sample_graph/$sample.vg \
		> evaluation/graphs/sample_graph/$sample.gfa

        echo ""
        echo "--------------------------------------"
        echo "Running vg autoindex"
        echo ""
	vg autoindex --workflow map \
		--prefix evaluation/graphs/sample_graph/$sample \
		--gfa evaluation/graphs/sample_graph/$sample.gfa \
        	-T tmp
}

export -f make_sample_graph

find vg_call/*filtered.vcf.gz \
	| parallel -j 12 --joblog $pjl/make_sample_graph.txt make_sample_graph {}

#slurm-4473361.out

COMMENT

#=================================================================================
# map reads to sample graph
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

map_to_sample_graph() {
        echo ""
        echo "--------------------------------------"
        read_name=$(basename $1 | awk -F'_' '{print $1}')
        read_2=$(echo $1 | sed 's/_1/_2/')

        echo $read_name
        echo $1
        echo $read_2
        echo "--------------------------------------"

        vg map -t 10 \
		-d evaluation/graphs/sample_graph/$read_name \
		-f $1 -f $read_2 \
                > evaluation/mappings/sample_graph/$read_name.mapped.gam

}

export -f map_to_sample_graph

find sequencing_reads/*_1.fastq \
        | parallel -j 12 --joblog $pjl/map_to_sample_graph.txt map_to_sample_graph {}

#slurm-4473790.out

COMMENT

#=================================================================================
# map reads to empty sample graph
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

map_to_empty_sample_graph() {
        echo ""
        echo "--------------------------------------"
        read_name=$(basename $1 | awk -F'_' '{print $1}')
        read_2=$(echo $1 | sed 's/_1/_2/')

        echo $read_name
        echo $1
        echo $read_2
        echo "--------------------------------------"

        vg map -t 10 \
                -d evaluation/graphs/empty_sample_graph/linear_graph \
                -f $1 -f $read_2 \
                > evaluation/mappings/empty_sample_graph/$read_name.mapped.gam

}

export -f map_to_empty_sample_graph

find sequencing_reads/*_1.fastq \
        | parallel -j 12 --joblog $pjl/map_to_empty_sample_graph.txt map_to_empty_sample_graph {}

#slurm-4474040.out

COMMENT

#=================================================================================
# collect identities and mapqs from sample graphs and linear graphs
#=================================================================================

: << 'COMMENT'

module load anaconda3/2023.09
conda activate vg_env

collect_identities_n_mapqs() {
        echo ""
        echo "--------------------------------------"
        sample=$(basename $1 | awk -F'.' '{print $1}')

        echo $sample
        echo $1
	echo "sample_graph"
        echo "--------------------------------------"

	echo "collecting identities"
	vg view -aj $1 \
		| jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' \
		> evaluation/mappings/sample_graph/$sample.identites.tsv 

	echo ""
	echo "collecting mapqs"
	vg view -aj $1 \
		| jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' \
		> evaluation/mappings/sample_graph/$sample.mapqs.tsv

	echo ""
	echo ""
	echo "--------------------------------------"
	echo "empty sample graph"
	echo ""
	echo "collecting identities"
	vg view -aj evaluation/mappings/empty_sample_graph/$sample.mapped.gam \
                | jq -rc '[.name, if .identity == null then 0 else .identity end ] | @tsv' \
		> evaluation/mappings/empty_sample_graph/$sample.identites.tsv

	echo ""
        echo "collecting mapqs"
        vg view -aj evaluation/mappings/empty_sample_graph/$sample.mapped.gam \
                | jq -rc '[.name, if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv' \
                > evaluation/mappings/empty_sample_graph/$sample.mapqs.tsv
}

export -f collect_identities_n_mapqs

find evaluation/mappings/sample_graph/*.gam \
        | parallel -j 12 --joblog $pjl/collect_identities_n_mapqs.txt collect_identities_n_mapqs {}

#slurm-4474280.out

COMMENT

#=================================================================================
# sort identities and mapqs files
#=================================================================================

: << 'COMMENT'

for i in evaluation/mappings/sample_graph/*.tsv evaluation/mappings/empty_sample_graph/*.tsv; do
	echo ""
	echo $i
	output=$(echo $i | sed 's/tsv/sorted.tsv/')
	echo $output	

	sort -k 1b,1 $i > $output
done

COMMENT

#=================================================================================
# join identities and mapqs; filter different 
#=================================================================================

: << 'COMMENT'

join_identities_n_mapqs() {
	echo ""
	echo ""
        echo "--------------------------------------"
        sample_name=$(basename $1 | awk -F'.' '{print $1}')
        empty_sample_graph_identities=$(echo $1 | sed 's/sample_graph/empty_sample_graph/')

	echo $sample_name
	echo $1
	echo $empty_sample_graph_identities
        echo "--------------------------------------"
	echo "Joining identities"
	join -j 1 -t $'\t' $1 $empty_sample_graph_identities > evaluation/mappings/$sample_name.identities.tsv

	echo ""
        echo "--------------------------------------"
	sample_graph_mapqs=$(echo $1 | sed 's/identites/mapqs/')
        empty_sample_graph_mapqs=$(echo $sample_graph_mapqs | sed 's/sample_graph/empty_sample_graph/')

        echo $sample_graph_mapqs
        echo $empty_sample_graph_mapqs
        echo "--------------------------------------"
        echo "Joining mapqs"
	join -j 1 -t $'\t' $sample_graph_mapqs $empty_sample_graph_mapqs > evaluation/mappings/$sample_name.mapqs.tsv

	echo ""
	echo "Joining identities and mapqs"
	join -j 1 -t $'\t' evaluation/mappings/$sample_name.identities.tsv evaluation/mappings/$sample_name.mapqs.tsv \
		> evaluation/mappings/$sample_name.ids_n_mapqs.tsv

	#rm evaluation/mappings/$sample_name.identities.tsv
	#rm evaluation/mappings/$sample_name.mapqs.tsv

	# filter different
	awk '!($2 == $3 && $4 == $5)' evaluation/mappings/$sample_name.ids_n_mapqs.tsv \
		> evaluation/mappings/$sample_name.ids_n_mapqs.different.tsv
}

export -f join_identities_n_mapqs

find evaluation/mappings/sample_graph/*identites.sorted.tsv \
	| parallel -j 12 --joblog $pjl/join_identities_n_mapqs.txt join_identities_n_mapqs {}

#slurm-______.out

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
