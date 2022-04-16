'''
Author : Usman Ghani

This script is used to greedily cluster the models in a directory
based on their pairwise interface RMSDs. The models in the directory
need to have aligned sequence residue numbers for this script to work
'''


from interface_clustering import interface_RMSD
from pymol import cmd
import numpy as np
import shutil
import json
import glob
import sys
import os 

            
def split_up_pdb_into_rec_and_lig(pdb_file, rec_chains, lig_chains, new_dir):
    
    '''Split up the pdb file into receptor and ligand using
    the provided chains.
    
    Args:
        pdb_file (str): Path to the pdb path of interest.
        rec_chains (str): Chains that constitute the receptor.
        lig_chains (str): Chains that constitute the ligand.
        new_dir (str): Path to a directory where the recpetor and
            ligand will be saved.
        
    Returns:
        rec_name (str): Path to the receptor.
        lig_name (str): Path to the ligand.
    
    '''
    
    rec_name = os.path.join(new_dir, os.path.basename(pdb_file[:-4]) + '_rec.pdb')
    lig_name = os.path.join(new_dir, os.path.basename(pdb_file[:-4]) + '_lig.pdb')
    
    cmd.reinitialize()
    cmd.load(pdb_file)
    cmd.create('complex', '*')
    
    cmd.select('rec_sele', f'complex and chain {"+".join(rec_chains)}')
    cmd.select('lig_sele', f'complex and chain {"+".join(lig_chains)}')
    
    cmd.save(rec_name, 'rec_sele')
    cmd.save(lig_name, 'lig_sele')
    
    return(rec_name, lig_name)

def prepare_input_list(file_list, rec_chains, lig_chains):
    
    '''Create an input list in the format required by 
    interface_RMSD.get_pairwise_RMSD.
    
    Args:
        file_list (list): List of pdbs.
        rec_chains (str): Chains that constitute the receptor.
        lig_chains (str): Chains that constitute the ligand.
        
    Returns:
        model_list (list): List of dictionaries. Each entry hasa te path to 
            receptor, the ligand, andd the name.
        temp_dir (str): Path to a temporary directory where the
            receptor and the ligand components for each model are
            saved.
        
    
    '''
    
    model_list = []
    temp_dir = os.path.join(os.path.dirname(file_list[0]), 'temp_parts')
    
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    
    for model in file_list:
        
        pdb = model
        name = os.path.basename(model)[:-4]
        rec_name, lig_name = split_up_pdb_into_rec_and_lig(pdb, rec_chains, lig_chains, temp_dir)
        dict_entry = {'rec':rec_name, 'lig':lig_name, 'name':name}
        
        model_list.append(dict_entry)
    
    return(model_list, temp_dir)
    
    
def get_clusters(RMSD_json, clustering_threshold=5.0):
    
    '''Load the RMSDs from the RMSD json and then performs greedy clustering
    to generate clusters and their centers.
    
    Args:
        RMSD_json (str): Path to JSON with the RMSDs.
        clustering_threshold (float): Threshold below which
            the RMSD indicates that two models are similar.
            
    Return:
        cluster_centers (list): List of names of cluster centers.
        named_clusters (dict of str:list): Maps cluster centers to 
            the members in their cluster. 
    
    '''
    
    
    with open(RMSD_json, 'r') as f:
        RMSD_array = json.load(f)
        
    dists = np.array(RMSD_array['dists'])
    names = RMSD_array['names']
    
    #To keep track of models already assigned to clusters
    done_set = set()
    
    clusters = {}
    while len(done_set) < len(names):
    
        #Dictionary for tracking close contacts in each loop
        temp_close_dict = {}
        for i, row in enumerate(dists):
            
            if i in done_set:
                continue
                
            close_contacts = np.where((row >= 0) & (row <= clustering_threshold))
            
            #Remove close contacts already assigned
            rel_indexes = [x for x in close_contacts[0] if x not in done_set]
            
            if len(rel_indexes) == 0:
                print('Something went wrong here. The cluster should at least contain itself as the cluster center. Adding itself manually')
                rel_indexes = [i]
                
            temp_close_dict[i] = rel_indexes
            
        #Now choose the largest cluster
        poss_clusters = list(temp_close_dict)
        sizes = [len(temp_close_dict[i]) for i in poss_clusters]
        max_size = max(sizes)
        new_cluster_index = poss_clusters[sizes.index(max_size)]
        
        #Add the newly found cluster to the done set
        done_set = done_set | set(temp_close_dict[new_cluster_index]) 
        
        clusters[new_cluster_index] = temp_close_dict[new_cluster_index]
        
    cluster_centers = [names[i] for i in clusters]
    
    #Switch from indexes to actual model names
    named_clusters = {}
    for cluster in clusters:
        name = names[cluster]
        members = [names[i] for i in clusters[cluster]]
        named_clusters[name] = members
    
    return(cluster_centers, named_clusters)
        
def cluster(input_dir, rec_chains, lig_chains, 
            interface_threshold, clustering_threshold, 
            RMSD_json, cluster_json):

    '''Cluster all the models in the input dir. 
    
    All chains other than the rec and lig chains are ignored. 
    Interface threshold defines the distance below which the rec 
    and lig residue atoms are considered to be in contact. 
    Clusutering threshold defines the RMSD below which a model is 
    considerd to be similar to another. The output RMSDs are saved 
    in RMSD_json and the cluster are saved in cluster_json.
    
    Args:
        input_dir (str): Path where the relevant models are 
            saved in pdb format.
        rec_chains (str): Chains that constitute the receptor.
        lig_chains (str): Chains that constitute the ligand.
        interface_threshold (float): Distance threshold below
            which receptor and ligand residues are considered to 
            be in contact and so, at the interface.
        clustering_threshold (float): Threshold below which
            the RMSD indicates that two models are similar.
        RMSD_json (str): Path to JSON where RMSDs are saved.
        cluster_json (str): Path where clustering results are
            saved.
            
    Returns:
        None
    
    '''
        
    file_list = glob.glob(os.path.join(input_dir, '*.pdb'))
    model_list, temp_dir = prepare_input_list(file_list, rec_chains, lig_chains)
    
    if os.path.exists(RMSD_json):
        print(f'Using {RMSD_json} for interface RMSDs...') 
    else:
        print('Calculating interface RMSDs...')
        
        try:
            dists = interface_RMSD.get_pairwise_RMSD(model_list, interface_threshold)
            interface_RMSD.make_result_json(model_list, dists, RMSD_json)
        except:
            print('Interface RMSD calculation failed')
            sys.exit(1)
        
    print(f'Clustering...')   
    cluster_centers, named_clusters = get_clusters(RMSD_json, clustering_threshold=clustering_threshold)
    shutil.rmtree(temp_dir)
    
    print('Cluster centers are:')
    for center in cluster_centers:
        print(center)
    
    with open(cluster_json, 'w') as f:
        json.dump(named_clusters, f)
        
        
        
        
        
        


