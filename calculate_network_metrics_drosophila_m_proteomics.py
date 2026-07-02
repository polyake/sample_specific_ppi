import pandas as pd
import numpy as np
import networkx as nx
import math
from collections import defaultdict
import os
import resource
from multiprocessing import Pool

rsrc = resource.RLIMIT_AS
resource.setrlimit(rsrc, (41474836480, 41474836480)) ## all usage limit set to ~40GB

## https://stackoverflow.com/questions/70858169/networkx-entropy-of-subgraphs-generated-from-detected-communities
def degree_distribution(G):
    vk = dict(G.degree())
    vk = list(vk.values()) # we get only the degree values
    
    Pk=defaultdict(float)
    for k in vk:
        Pk[k] = Pk[k] + 1
    Pk_values = np.asarray(list(Pk.values()))
    Pk_final = Pk_values/sum(Pk_values) # the sum of the elements of P(k) must to be equal to one
    
    return Pk_final,vk

def network_entropy(G, Pk, vk):
    ## Shannon -- for the degree distributiom
    H = 0
    for p in Pk:
        if(p > 0):
            H = H - p*math.log(p, 2)
##            
    ## random_walker -- for the degrees
    vk_nonzero = [i for i in vk if i != 0]
    ln_vk = np.log2(vk_nonzero)
    sum_ln_vk = sum(ln_vk)
    const = len(vk) * np.log2(len(vk)-1)
    H_rw = sum_ln_vk/const
    
    A = nx.adjacency_matrix(G)
    A = A.toarray()
    lbd = max(np.linalg.eigvals(A))
    H_ks = math.log(lbd)
    
    return H, H_rw, H_ks

def density_calc(G, degrees, num_nodes):
    avg = np.average(degrees)
    const = num_nodes-1
    return avg/const

def centralization_calc(G, degrees, num_nodes):
    ### ~ max(degree)-avg(degree)
    density = density_calc(G, degrees, num_nodes)
    
    if num_nodes>2:
        norm_const = num_nodes/(num_nodes-2)
        centralization = ((max(degrees)/(num_nodes-1))-density)*norm_const
    else:
        centralization = np.nan
    return density, centralization

def heterogeneity_calc(G, degrees):
    num = math.sqrt(np.var(degrees))
    denom = np.average(degrees)
    return num/denom

def skewness_of_degdist(G, degrees):
    m = np.average(degrees)
    sigma = np.std(degrees)
    skewness = (np.average((np.power(degrees - m, 3))))/(sigma**3)
    return skewness


def metric_calculator(G, file_name):
    ### whole graph
    Pk, degrees = degree_distribution(G)
    num_nodes = len(degrees)
    density, centralization = centralization_calc(G, degrees, num_nodes)
    heterogeneity = heterogeneity_calc(G, degrees)
    entropy_shannon, entropy_rw, entropy_ks = network_entropy(G, Pk, degrees)
    clustering = nx.clustering(G)
    clustering_coefs = list(clustering.values())
    node_bc_dict = nx.betweenness_centrality(G)
    node_bc = list(node_bc_dict.values())
    edge_bc = list(nx.edge_betweenness_centrality(G).values())

    dat = [file_name,nx.number_of_nodes(G), nx.number_of_edges(G), nx.number_connected_components(G), 
           np.average(degrees), max(degrees), density, centralization, heterogeneity, nx.transitivity(G),
           np.average(clustering_coefs), max(clustering_coefs), skewness_of_degdist(G, degrees),
           nx.degree_assortativity_coefficient(G), entropy_shannon, entropy_rw, entropy_ks,
          max(node_bc), np.average(node_bc), max(edge_bc), np.average(edge_bc)]
    column_names = ['file_name', 'num_nodes', 'num_edges', 'num_components', 
                    'avg_deg', 'max_deg', 'density', 'centralization', 'heterogeneity', 'glob_clust', 'avg_clust', 'max_clust', 
                   'skew_deg_dist', 'assortativity_degree' , 'entropy_shannon', 'entropy_rw', "entropy_ks",
                   'node_bc_max', 'node_bc_avg', 'edge_bc_max', 'edge_bc_avg']
    
    ### node-level metrics
    node_degrees = pd.DataFrame(G.degree(), columns=["index", "degree"])
    node_degrees = pd.concat([node_degrees,pd.DataFrame({"index":["file_name"], "degree":[file_name]})], ignore_index=True)
    node_degrees = node_degrees.set_index("index")
    node_degrees = node_degrees.T
    degree_centrality = nx.degree_centrality(G)
    degree_centrality = pd.DataFrame([degree_centrality]).T
    degree_centrality = degree_centrality.reset_index()
    degree_centrality.columns = ["index","degree_centrality"]
    degree_centrality = pd.concat([degree_centrality,pd.DataFrame({"index":["file_name"], "degree":[file_name]})], ignore_index=True)
    degree_centrality = degree_centrality.set_index("index")
    degree_centrality = degree_centrality.T
    clust_df = pd.DataFrame([clustering]).T
    clust_df = clust_df.reset_index()
    clust_df.columns = ["index","clustering"]
    clust_df = pd.concat([clust_df,pd.DataFrame({"index":["file_name"], "degree":[file_name]})], ignore_index=True)
    clust_df = clust_df.set_index("index")
    clust_df = clust_df.T
    node_bc_df = pd.DataFrame([node_bc_dict]).T
    node_bc_df = node_bc_df.reset_index()
    node_bc_df.columns = ["index","node_bc"]
    node_bc_df = pd.concat([node_bc_df,pd.DataFrame({"index":["file_name"], "degree":[file_name]})], ignore_index=True)
    node_bc_df = node_bc_df.set_index("index")
    node_bc_df = node_bc_df.T
    
    
    dat_final = dat
    column_names_final = column_names
    
    #df_metrics = pd.DataFrame([dat_final], columns=column_names_final)
    
    return dat_final, column_names_final, node_degrees, degree_centrality, clust_df, node_bc_df


def run_metric_calculator(file):
    sample = file[8:-14]
    network_df = pd.read_csv("network_data/drosophila_m/"+file)
    G = nx.from_pandas_edgelist(network_df, source="protein1", target="protein2")
    #G.add_nodes_from([(dict(d)["protein1"], {"protein_freq":dict(d)[sample+"_x"], "protein_name":dict(d)["protein1_name"]})
    #                      for n, d in network_df.iterrows()])
    #G.add_nodes_from([(dict(d)["protein2"], {"protein_freq":dict(d)[sample+"_y"], "protein_name":dict(d)["protein2_name"]})
    #                      for n, d in network_df.iterrows()])

    metrics_data, metrics_cols, node_degrees, degree_centrality, clust_df, node_bc_df = metric_calculator(G, file)           
        
    dat_final = metrics_data
    column_names_final = metrics_cols
        
    metrics_df = pd.DataFrame([dat_final], columns=column_names_final)
        
        
    metrics_df.to_csv("metrics_data/drosophila_m/protein_network_metrics_"+sample+".csv")
    degree_centrality.to_csv("metrics_data/drosophila_m/protein_degree_centralities_"+sample+".csv")
    node_degrees.to_csv("metrics_data/drosophila_m/protein_degrees_"+sample+".csv")
    clust_df.to_csv("metrics_data/drosophila_m/protein_clustering_coefs_"+sample+".csv")
    node_bc_df.to_csv("metrics_data/drosophila_m/protein_node_bcs_"+sample+".csv")
    print(file)
    
    
#################################################################################
### run metric calculation
    
files = os.listdir("network_data/drosophila_m/")
files_protein = []
for file in files:
    if "protein" in file:
        files_protein.append(file)
          
with Pool(processes=5) as pool:
    pool.map(run_metric_calculator, files_protein)