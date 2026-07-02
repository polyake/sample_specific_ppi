import pandas as pd
import numpy as np
import h5py

print("It works")


ppi = pd.read_csv("data/7227.protein.physical.links.full.v12.0.txt", sep=" ")
data_all = pd.read_csv("data/drosophila_m_lfq_mapped_to_string.csv")

samples = data_all.columns[1:]

for sample in samples:
    data = data_all[["#string_protein_id", sample]]
    data = data.set_index("#string_protein_id")
    data = data.T
    data = data.fillna(0)

    ppi_sub = ppi[(ppi["protein1"].isin(data.columns))&(ppi["protein2"].isin(data.columns))]
    ppi_sub_high = ppi_sub[(ppi_sub["experiments"]>700)|(ppi_sub["database"]>700)][["protein1","protein2","experiments","database","combined_score"]]

    for i in range(len(data)):
        aux = ppi_sub_high.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein1", right_index=True)
        ppi_0 = aux.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein2", right_index=True)
        #ppi_0 = ppi_0.rename({str(i)+"_x":"protein1_freq", str(i)+"_y":"protein2_freq"}, axis=1)
        ppi_0.to_csv("./network_data/drosophila_m/protein_"+sample[14:]+"_"+str(i)+"_network.csv", index=False)

    print(sample)