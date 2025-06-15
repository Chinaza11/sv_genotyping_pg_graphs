# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 14:11:56 2025

@author: Chinaza
"""

import pandas as pd, os

os.chdir("G:/My Drive/PhD/project/structural_variant_genotyping_with_pangenome_graph/data_analysis/data_and_result/main/graph_minigraph-cactus/terra_raw")

df_genotype_info = pd.read_excel("Songsomboon_table_revised1.xlsx", sheet_name="Supplement table 1", skiprows=2)

df_eigenvec = pd.read_csv("pca_no_filter/plink.eigenvec", sep="\t")
df_ids = df_eigenvec.iloc[:, [0]]

# =============================================================================
# 
# =============================================================================

for index, genotype in df_genotype_info['Genotype'].items():
    if genotype.startswith("PI"):
        matches = df_ids[df_ids["#IID"].str.contains(genotype.split(" ")[1], case=False, na=False)]
        for i in range(len(matches)):
            df_ids.loc[matches.index[i], 'Genotype'] = genotype
    else:
        matches = df_ids[df_ids["#IID"].str.contains(genotype, case=False, na=False)]
        for i in range(len(matches)):
            df_ids.loc[matches.index[i], 'Genotype'] = genotype

df_ids.to_csv('df_ids.csv', index=False)
