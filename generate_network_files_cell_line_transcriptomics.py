import pandas as pd
import numpy as np
import h5py

print("It works")

string_meta = pd.read_csv("9606.protein.info.v12.0.txt", sep="\t")
mapping = pd.read_csv("human_MHCC97H_gene_to_string_mapping.csv")
mapping_forward = dict(zip(mapping["preferred_name"], mapping["#string_protein_id"]))
mapping_back = dict(zip(mapping["#string_protein_id"], mapping["preferred_name"]))
ppi = pd.read_csv("9606.protein.physical.links.full.v12.0.txt", sep=" ")


samples = ["GSM7454068_V350003741_L01_81.txt", "GSM7454069_V350003741_L01_82.txt", "GSM7454070_V350003741_L01_83.txt",
          "GSM7454071_V350003741_L01_84.txt", "GSM7454072_V350003741_L01_85.txt", "GSM7454073_V350003741_L01_86.txt",
          "GSM7454074_V350003741_L01_87.txt", "GSM7454075_V350003741_L01_88.txt", "GSM7454076_V350003741_L01_73.txt"]

for sample in samples:
    data = pd.read_csv("data/"+sample, sep=" ")
    data = data.set_index("Gene")
    data = data.T
    data = data.rename(mapping_forward, axis=1)

    ppi_sub = ppi[(ppi["protein1"].isin(data.columns))&(ppi["protein2"].isin(data.columns))]
    ppi_sub_high = ppi_sub[(ppi_sub["experiments"]>700)|(ppi_sub["database"]>700)][["protein1","protein2","experiments","database","combined_score"]]
    ppi_sub_high['protein1_name'] = ppi_sub_high['protein1'].map(mapping_back)
    ppi_sub_high['protein2_name'] = ppi_sub_high['protein2'].map(mapping_back)

    for i in range(len(data)):
        aux = ppi_sub_high.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein1", right_index=True)
        ppi_0 = aux.merge(pd.DataFrame(data.iloc[i][data.iloc[i]!=0]), left_on="protein2", right_index=True)
        ppi_0 = ppi_0.rename({str(i)+"_x":"protein1_freq", str(i)+"_y":"protein2_freq"}, axis=1)
        ppi_0.to_csv("./network_data/"+sample[:-4]+"_"+str(i)+"_network.csv", index=False)
    
        print(i)
    print(sample)