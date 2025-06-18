#! /bin/bash

#=================================================================================
# run map_reads2 function in 01_graph_minigraph-cactus_yeast.slurm on Wild_Sorghum_WGS files
#=================================================================================

: << "COMMENT"

for sample in pi302118_1 pi302267_1 pi329250_1 pi329251_1 pi329252_1 pi369484_1 pi369487_1 pi524718_1 pi532564_1 pi532565_1 pi532566_1 pi532568_1 pi535995_1;
do
        echo $sample
        sed "s#Grif16309_1#$sample#g" 01_graph_minigraph-cactus.slurm > $sample.01_graph_minigraph-cactus.slurm
        sbatch $sample.01_graph_minigraph-cactus.slurm
        rm $sample.01_graph_minigraph-cactus.slurm
done

COMMENT

#=================================================================================
# run map_reads2 function in 01_graph_minigraph-cactus_yeast.slurm on TERRA_RAW files
#=================================================================================

: << "COMMENT"

for file in /projects/cooper_research1/TERRA_RAW/*_R1.fastq.gz;
do
	sample=$(basename $file | awk -F'.' '{print $1}')
        echo $sample
        sed "s#placeholder#$sample#g" 01_graph_minigraph-cactus.slurm > $sample.01_graph_minigraph-cactus.slurm
        sbatch $sample.01_graph_minigraph-cactus.slurm
        rm $sample.01_graph_minigraph-cactus.slurm
done

COMMENT

#=================================================================================
# run vg_call function in 01_graph_minigraph-cactus_yeast.slurm on TERRA_RAW files
#=================================================================================

#: << "COMMENT"

for file in mappings_terra_raw/*.mapped.gam;
do
	echo $file
	sample=$(basename "$file" | awk -F'.' '{print $1}')
	echo $sample
        sed "s#placeholder#$sample#g" 01_graph_minigraph-cactus.slurm > $sample.01_graph_minigraph-cactus.slurm
        sbatch $sample.01_graph_minigraph-cactus.slurm
        rm $sample.01_graph_minigraph-cactus.slurm
done

#COMMENT

#=================================================================================
#
#=================================================================================

#=================================================================================
#
#=================================================================================

#=================================================================================
#
#=================================================================================
