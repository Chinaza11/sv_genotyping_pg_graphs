#! /bin/bash

#=================================================================================
# run map_reads function in 01_graph_minigraph-cactus_yeast.slurm
#=================================================================================

#: << COMMENT

for sample in pi302118_1 pi302267_1 pi329250_1 pi329251_1 pi329252_1 pi369484_1 pi369487_1 pi524718_1 pi532564_1 pi532565_1 pi532566_1 pi532568_1 pi535995_1;
do
        echo $sample
        sed "s#Grif16309_1#$sample#g" 01_graph_minigraph-cactus_yeast.slurm > $sample.01_graph_minigraph-cactus_yeast.slurm
        sbatch $sample.01_graph_minigraph-cactus_yeast.slurm
        rm $sample.01_graph_minigraph-cactus_yeast.slurm
done

#COMMENT

#=================================================================================
#
#=================================================================================

#=================================================================================
#
#=================================================================================
