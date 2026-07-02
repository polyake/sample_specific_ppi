import pandas as pd
import numpy as np
import h5py

print("It works")

data_all = pd.read_csv("data/1-9DDA_mapped.csv")

string_meta = pd.read_csv("9606.protein.info.v12.0.txt", sep="\t")
mapping_back = dict(zip(data_all["#string_protein_id"], data_all["alias"]))
ppi = pd.read_csv("9606.protein.physical.links.full.v12.0.txt", sep=" ")


samples = ["g-1 ", "g-3", "g-4", "g-5", "g-6", "g-7", "g-8", "g-9"]

for sample in samples:
    data = data_all[["#string_protein_id", sample]]
    data = data.set_index("#string_protein_id")
    data = data.T
    data = data.fillna(0)

    ppi_sub = ppi[(ppi["protein1"].isin(data.columns))&(ppi["protein2"].isin(data.columns))]
    ppi_sub_high = ppi_sub[(ppi_sub["experiments"]>700)|(ppi_sub["database"]>700)][["protein1","protein2","experiments","database","combined_score"]]
    ppi_sub_high['protein1_name'] = ppi_sub_high['protein1'].map(mapping_back)
    ppi_sub_high['protein2_name'] = ppi_sub_high['protein2'].map(mapping_back)

    for i in range(len(data)):
        aux = ppi_sub_high.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein1", right_index=True)
        ppi_0 = aux.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein2", right_index=True)
        #ppi_0 = ppi_0.rename({str(i)+"_x":"protein1_freq", str(i)+"_y":"protein2_freq"}, axis=1)
        ppi_0.to_csv("./network_data/protein_"+sample+"_"+str(i)+"_network.csv", index=False)
    
        print(i)
    print(sample)