import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import os
import statsmodels.api as sm
import scipy
import resource

import warnings
warnings.filterwarnings("ignore")

rsrc = resource.RLIMIT_AS
resource.setrlimit(rsrc, (41474836480, 41474836480)) ## all usage limit set to ~40GB

def fit_separate_models(data_t, data_p, filename):
    common_nodes = set(data_t.columns)&set(data_p.columns)-set(["file_name", "replicate", "Unnamed: 0"])
    
    df_final = pd.DataFrame()
    for i,node in enumerate(common_nodes):
        Y = data_t[node]
        X = sm.add_constant(data_t[["file_name"]])
        model = sm.OLS(Y,X)
        results_t = model.fit()
    
        Y = data_p[node]
        X = sm.add_constant(data_p[["file_name"]])
        model = sm.OLS(Y,X)
        results_p = model.fit()
        row = {"node":[node],
               "transcriptomic_coef":[results_t.params["file_name"]], "transcriptomic_pvalue":[results_t.pvalues["file_name"]],
               "proteomic_coef":[results_p.params["file_name"]], "proteomic_pvalue":[results_p.pvalues["file_name"]]}
        df_final = pd.concat([df_final, pd.DataFrame(row)])
        
    df_final.to_csv("local_tests/"+filename+".csv", index=False)
    
    
def interaction_test(data_t, data_p, filename):
    common_nodes = set(data_t.columns)&set(data_p.columns)-set(["file_name", "replicate", "Unnamed: 0"])
    
    df_final = pd.DataFrame()
    for i,node in enumerate(common_nodes):
        aux = data_t[["file_name",node]]
        aux["group"] = [0 for i in range(len(aux))]

        aux2 = data_p[["file_name",node]]
        aux2["group"] = [1 for i in range(len(aux2))]

        aux_final = pd.concat([aux,aux2], ignore_index=True)
        aux_final["interaction"] = aux_final["file_name"]*aux_final["group"]
    
        Y = aux_final[node]
        X = sm.add_constant(aux_final[["file_name","group","interaction"]])
        model = sm.OLS(Y,X)
        results = model.fit()
        row = {"node":[node],
               "interaction_coef":[results.params["interaction"]], "interaction_pvalue":[results.pvalues["interaction"]]}
        df_final = pd.concat([df_final, pd.DataFrame(row)])
        
    df_final.to_csv("local_tests/"+filename+".csv", index=False)
    

def random_label_switch(df, dataset, seed=0):
    g1, g2 = df["group"].unique()
    switch_map = {g1:g2, g2:g1}
    
    np.random.seed(seed)
    
    if dataset=="cell_line":
        pair_switch = df[["file_name"]].drop_duplicates().assign(switch=lambda x:np.random.choice([True,False], size=len(x)))
        df = df.merge(pair_switch, on=["file_name"])
    elif dataset=="embryo":
        pair_switch = df[["file_name", "replicate"]].drop_duplicates().assign(switch=lambda x:np.random.choice([True,False], size=len(x)))
        df = df.merge(pair_switch, on=["file_name", "replicate"])
    
    mask = df["switch"]
    df.loc[mask, "group"] = df.loc[mask, "group"].map(switch_map)
    
    return df

def random_slope_gen(df, node, dataset, seed):
    switched_df = random_label_switch(df, dataset, seed)
    
    Y = switched_df[switched_df["group"]==0][node]
    X = sm.add_constant(switched_df[switched_df["group"]==0][["file_name"]])
    model = sm.OLS(Y,X)
    results = model.fit()
    slope1 = results.params["file_name"]
    
    Y = switched_df[switched_df["group"]==1][node]
    X = sm.add_constant(switched_df[switched_df["group"]==1][["file_name"]])
    model = sm.OLS(Y,X)
    results = model.fit()
    slope2 = results.params["file_name"]
    
    return abs(slope1-slope2)
    
def run_randomized_test(data_t, data_p, filename, dataset):
    common_nodes = set(data_t.columns)&set(data_p.columns)-set(["file_name", "replicate", "Unnamed: 0"])
    
    df_final = pd.DataFrame()
    for i,node in enumerate(common_nodes):
        #################################
        ## getting the random slope differences
        if dataset=="cell_line":
            aux = data_t[["file_name",node]]
            aux["group"] = [0 for i in range(len(aux))]
            aux2 = data_p[["file_name",node]]
            aux2["group"] = [1 for i in range(len(aux2))]
        elif dataset=="embryo":
            aux = data_t[["file_name", "replicate",node]]
            aux["group"] = [0 for i in range(len(aux))]
            aux2 = data_p[["file_name", "replicate",node]]
            aux2["group"] = [1 for i in range(len(aux2))]
            
        merged_df = pd.concat([aux,aux2], ignore_index=True)
        
        random_slope_diffs = []
        for seed in range(1000):
            slope_diff = random_slope_gen(merged_df, node, dataset, seed)
            random_slope_diffs.append(slope_diff)
        
        ################################
        ## getting the "null" -- observed slope difference
        Y = merged_df[merged_df["group"]==0][node]
        X = sm.add_constant(merged_df[merged_df["group"]==0][["file_name"]])
        model = sm.OLS(Y,X)
        results = model.fit()
        slope1 = results.params["file_name"]
    
        Y = merged_df[merged_df["group"]==1][node]
        X = sm.add_constant(merged_df[merged_df["group"]==1][["file_name"]])
        model = sm.OLS(Y,X)
        results = model.fit()
        slope2 = results.params["file_name"]
    
        null = abs(slope1-slope2)
        ###############################
        count = len([x for x in random_slope_diffs if x < null])
        pvalue = count/1000
        
        row = {"node":[node], "count":[count], "p_value":[pvalue]}
        df_final = pd.concat([df_final, pd.DataFrame(row)])
        
    df_final.to_csv("local_tests/"+filename+".csv", index=False)
    
    
################# run
metrics = ["degree_centralities", "clustering_coefs", "node_bcs"]

## cell_line data
file_names = {"GSM7454068_V350003741_L01_81_0_network.csv":1, "GSM7454069_V350003741_L01_82_0_network.csv":2,
              "GSM7454070_V350003741_L01_83_0_network.csv":3,
          "GSM7454071_V350003741_L01_84_0_network.csv":4, "GSM7454072_V350003741_L01_85_0_network.csv":5,
              "GSM7454073_V350003741_L01_86_0_network.csv":6,
          "GSM7454074_V350003741_L01_87_0_network.csv":7, "GSM7454075_V350003741_L01_88_0_network.csv":8,
              "GSM7454076_V350003741_L01_73_0_network.csv":9}

file_names_protein = {"protein_g-1 _network.csv":1, "protein_g-3_network.csv":3,
          "protein_g-4_network.csv":4, "protein_g-5_network.csv":5,
              "protein_g-6_network.csv":6,
          "protein_g-7_network.csv":7, "protein_g-8_network.csv":8,
              "protein_g-9_network.csv":9}

for metric in metrics:
    data = pd.read_csv("metrics_data/data/"+metric+".csv", index_col=False)
    data_protein = pd.read_csv("metrics_data/data/protein_"+metric+".csv", index_col=False)
    data["file_name"] = data["file_name"].map(file_names)
    data_protein["file_name"] = data_protein["file_name"].map(file_names_protein)
    
    data = data.fillna(0)
    data_protein = data_protein.fillna(0)
    
    fit_separate_models(data, data_protein, "cell_line_"+metric+"_tendency.csv")
    print("Cell line tendency done!")
    interaction_test(data, data_protein, "cell_line_"+metric+"_interaction.csv")
    print("Cell line interaction done!")
    run_randomized_test(data, data_protein, "cell_line_"+metric+"_randomized.csv", "cell_line")
    print("Cell line randomized done!")
    
## embryo data
for metric in metrics:
    data = pd.read_csv("metrics_data/data/drosophila_m/"+metric+".csv", index_col=False)
    data_protein = pd.read_csv("metrics_data/data/drosophila_m/protein_"+metric+".csv", index_col=False)
    data["replicate"] = data["file_name"].str.split("_", expand=True)[4].astype(int)
    data_protein["replicate"] = data_protein["file_name"].str[-5].astype(int)
    data["file_name"] = data["file_name"].str.split("_",expand=True)[3].str[:-1].astype(int)
    data_protein["file_name"] = data_protein["file_name"].str[-9:-7].astype(int)
    
    data = data.fillna(0)
    data_protein = data_protein.fillna(0)
    
    fit_separate_models(data, data_protein, "embryo_"+metric+"_tendency.csv")
    print("Embryo tendency done!")
    interaction_test(data, data_protein, "embryo_"+metric+"_interaction.csv")
    print("Embryo interaction done!")
    run_randomized_test(data, data_protein, "embryo_"+metric+"_randomized.csv", "embryo")
    print("Embryo randomized done!")