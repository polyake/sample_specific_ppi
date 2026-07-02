import pandas as pd
import numpy as np
import h5py
import os

print("It works")

mapping = pd.read_csv("data/drosophila_m_gene_mapping.csv")
mapping_forward = dict(zip(mapping["alias"], mapping["#string_protein_id"]))
mapping_back = dict(zip(mapping["#string_protein_id"], mapping["alias"]))
ppi = pd.read_csv("data/7227.protein.physical.links.full.v12.0.txt", sep=" ")

files = os.listdir("data/GSE121160_RAW/")


for file in files:
    data = pd.read_csv("data/GSE121160_RAW/"+file, sep="\t", header=None)
    data = data.set_index(0).T
    data = data.rename(mapping_forward, axis=1)

    ppi_sub = ppi[(ppi["protein1"].isin(data.columns))&(ppi["protein2"].isin(data.columns))]
    ppi_sub_high = ppi_sub[(ppi_sub["experiments"]>700)|(ppi_sub["database"]>700)][["protein1","protein2","experiments","database","combined_score"]]
    #ppi_sub_high['protein1_name'] = ppi_sub_high['protein1'].map(mapping_back)
    #ppi_sub_high['protein2_name'] = ppi_sub_high['protein2'].map(mapping_back)

    for i in range(len(data)):
        aux = ppi_sub_high.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein1", right_index=True)
        ppi_0 = aux.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein2", right_index=True)
        #ppi_0 = ppi_0.rename({str(i)+"_x":"protein1_freq", str(i)+"_y":"protein2_freq"}, axis=1)
        ppi_0.to_csv("./network_data/drosophila_m/"+file[:-7]+"_"+str(i)+"_network.csv", index=False)
    
    print(file)