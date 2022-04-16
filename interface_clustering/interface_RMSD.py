'''
Author : Usman Ghani

The script is meant for clustering models. The goal is to identify the interacting residues 
for the surfaces and then use the carbon alpha atom distances of the interacting residues. 
The thing of interest here is that this kind of clustering needs structural aligment. 
It uses pymol's align function for this.
'''
from scipy.spatial.distance import cdist
from prody import parsePDB, writePDB
import concurrent.futures
from pymol import cmd
import numpy as np
import prody
import json
import time 
import os

prody.confProDy(verbosity='none')

def get_interface(rec_coords, lig_coords, interface_threshold):

    '''Determine indices of atoms that are at the receptor-ligand
    interface
    
    Args:
        rec_coords (ndarray): Coordinates for each atom in the receptor.
        lig_coords (ndarray): Coordinates for each atom in the ligand.
        interface_threshold (float): Distance threshold below
            which receptor and ligand residues are considered to 
            be in contact and so, at the interface.
            
    Returns:
        rec_coord_indices (ndarray): Indices of the rec_coords where atoms are 
            within the distance threshold to any atom in lig_coords. 
        lig_coord_indices (ndarray): Indices of the lig_coords where atoms are 
            within the distance threshold to any atom in rec_coords.
    
    '''

    dists = cdist(rec_coords, lig_coords, 'sqeuclidean')
    
    lig_coord_indices = np.any(dists < (interface_threshold * interface_threshold), axis=0).nonzero()[0]
    rec_coord_indices = np.any(dists < (interface_threshold * interface_threshold), axis=1).nonzero()[0]
    
    return(lig_coord_indices, rec_coord_indices)
    
def identify_interface(rec_path, lig_path, interface_threshold=10.0):
    
    '''
    Determine the interacting residue numbers.
    
    Args:
        rec_path (str): Path to receptor model in pdb format.
        lig_path (str): Path to ligand model in pdb format. 
        interface_threshold (float): Distance threshold below
            which receptor and ligand residues are considered to 
            be in contact and so, at the interface.
            
    Returns:
        rec_interacting_res (set): Receptor residues at the interface with ligand.
        lig_interacting_res (set): Ligand residues at the interface with receptor.
    '''
    
    interface_dict = {}
    rec = parsePDB(rec_path)
    rec_coords = rec.getCoords()
    lig = parsePDB(lig_path)
    lig_coords = lig.getCoords()
    lig_coord_indices, rec_coord_indices = get_interface(rec_coords, lig_coords, interface_threshold)
    
    #Get interacting residues in the format "{chain}_{residue_number}"
    lig_interacting_res = set([f'{x[0]}_{str(x[1])}' for x in zip(list(lig.getChids()[lig_coord_indices]), list(lig.getResnums()[lig_coord_indices]))])
    rec_interacting_res = set([f'{x[0]}_{str(x[1])}' for x in zip(list(rec.getChids()[rec_coord_indices]), list(rec.getResnums()[rec_coord_indices]))])    
    
    return(rec_interacting_res, lig_interacting_res)
    
    
def get_RMSD(coord_array1, coord_array2):
    
    '''
    Calculate RMSD for the two arrays of coordinates. Ensure
    that the atoms are index matched.
    
    Args:
        coord_array1 (ndarry): Array of atom coordinaates 1.
        coord_array2 (ndarry): Array of atom coordinaates 2.
    
    Returns:
        RMSD (float): RMSD calculated using coord_array1 and
            coord_array2. 
    
    '''
    
    delta = coord_array1 - coord_array2
    delta_square_sum = (np.multiply(delta, delta).sum())/len(delta)
    RMSD = np.sqrt(delta_square_sum)
    
    return(RMSD)
    
    
def process_RMSD(rec1, rec2, lig1, lig2, rec_interface, lig_interface, job_id):

    '''
    Calculate RMSD.
    
    Args:
        rec1 (str): Path to receptor 1 pdb.
        rec2 (str): Path to receptor 2 pdb.
        lig1 (str): Path to ligand 1 pdb.
        lig2 (str): Path to ligand 2 pdb.
        rec_interface (list): Residue numbers for receptor residues at the
            interface with the ligand.
        lig_interface (list): Residue numbers for ligand residues at the
            interface with the receptor.
        job_id (int): Job number.
        
    Returns:
        RMSD (float): RMSD of the interface residues.
        
        
    '''
    
    #Account for the case where there is no interface
    if len(rec_interface) == 0:
        RMSD = -1
        return(RMSD)
        
    temp_rec2, temp_lig2 = build_temp_aligned_models(rec1, rec2, lig2, job_id)
    
    #Get coordinates for everything and combine them
    model1 = np.concatenate((get_interface_CA_cords(rec1, rec_interface), 
                             get_interface_CA_cords(lig1, lig_interface)),
                             axis=0)
                             
    model2 = np.concatenate((get_interface_CA_cords(temp_rec2, rec_interface),
                             get_interface_CA_cords(temp_lig2, lig_interface)), 
                             axis=0)            
    
    assert len(model1) == len(model2), f'CA atoms do not match for:\n {rec1} + {lig1} \n{rec2} + {lig2}\njob_id = {job_id}\nmodel1 = {len(model1)}, model2 = {len(model2)}'
    
    RMSD = get_RMSD(model1, model2)
    
    #Remove the aligned models
    os.remove(temp_rec2)
    os.remove(temp_lig2)
    
    return(RMSD)
    
def build_temp_aligned_models(rec1_path, rec2_path, lig2_path, job_id):

    '''Build structurally aligned models using PyMol.
    
    Args:
        rec1_path (str): Path to receptor 1 pdb.
        rec2_path (str): Path to receptor 2 pdb.
        lig2_path (str): Path to ligand 2 pdb.
        jobs_id (int): Job number.
    
    Returns:
        temp_rec (str): Path where aligned receptor 2 is saved.
        temp_lig (str): Path where aligned ligand 2 is saved.
    
    '''
    
    cmd.reinitialize()
    
    cmd.load(rec1_path, 'rec1')
    cmd.load(rec2_path, 'rec2')
    cmd.load(lig2_path, 'lig2')
    
    cmd.align('rec2', 'rec1')
    cmd.matrix_copy('rec2', 'lig2')
    
    destination = os.path.dirname(rec1_path)
    temp_rec = os.path.join(destination, 'temp_rec2_{}.pdb'.format(job_id))
    temp_lig = os.path.join(destination, 'temp_lig2_{}.pdb'.format(job_id))
    
    cmd.save(temp_rec, 'rec2')
    cmd.save(temp_lig, 'lig2')
    
    return(temp_rec, temp_lig)
    
def get_interface_CA_cords(pdb_path, interface):

    '''Build CA atom array for the interface residues
    of the model specified at the pdb path
    
    Args:
        pdb_path (str): Path to relevant pdb.
        interface (list): Residue numbers for residues at the
            interface.
        
    Returns:
        interface_CAs (ndarray): Atom coordinates of carbon alpha
            atoms of residues at the interface. 
    
    '''

    pdb = parsePDB(pdb_path)
    CAs = pdb.select('name CA')
    
    for i, (chain, res_nums) in enumerate(interface):
    
        interface_rec_string = (' ').join([str(ele) for ele in res_nums])
        interface_chain_CAs = CAs.select(f'chain {chain}').select('resnum {}'.format(interface_rec_string)).getCoords()
        
        if i == 0:
            interface_CAs = interface_chain_CAs
        else:
            interface_CAs = np.concatenate((interface_CAs, interface_chain_CAs), axis=0)
    
    return(interface_CAs)
    
    
def get_pairwise_RMSD(models, interface_threshold):

    '''
    Get pairwise RMSDs.
    
    Models is a list of dicts. Each entry in the list contains the receptor path ("rec"), 
    ligand path ("lig"), and a name ("name") for the model. The name is what is going to 
    be used to keep track of the targets.
    
    Args:
        models (list): List of dictionaries with an entry for each model. 
        interface_threshold (float): Distance threshold below
            which receptor and ligand residues are considered to 
            be in contact and so, at the interface.
            
    Returns:
        dists (ndarray): Array of pairwise RMSDs
    
    '''
    N = len(models)
    dists = np.zeros((N, N), dtype=np.float64)
    
    #Interface dict:
    interface_dict = {}
    for i, model in enumerate(models):
    
        rec = model['rec']
        lig = model['lig']
        rec_interacting_res, lig_interacting_res = identify_interface(rec, lig, interface_threshold=interface_threshold)
        interface_dict[i] = {'rec':rec_interacting_res, 'lig':lig_interacting_res}
    
    rec1_list = []
    lig1_list = []
    rec2_list = []
    lig2_list = []
    rec_interface_list = []
    lig_interface_list = []
    job_list = []
    count = 0
    for i in range(N):
        
        rec1 = models[i]['rec']
        lig1 = models[i]['lig']
        rec1_interacting_res = interface_dict[i]['rec'] 
        lig1_interacting_res = interface_dict[i]['lig'] 
        
        
        #Assume symmetrical distance matrix
        for j in range(i , N):
            
            rec2 = models[j]['rec']
            lig2 = models[j]['lig']
            rec2_interacting_res = interface_dict[j]['rec'] 
            lig2_interacting_res = interface_dict[j]['lig'] 
            
            #Consider the interface residues on both models
            rec_interface = list(rec1_interacting_res | rec2_interacting_res)
            lig_interface = list(lig1_interacting_res | lig2_interacting_res)
            
            #Split these by chain and have the chains and resnums be sorted
            lig_chain_interacting_res = {}
            for res in lig_interface:
                chain = res.split('_')[0]
                res = int(res.split('_')[1])
                if chain not in lig_chain_interacting_res:
                    lig_chain_interacting_res[chain] = []
                lig_chain_interacting_res[chain].append(res)
            
            rec_chain_interacting_res = {}
            for res in rec_interface:
                chain = res.split('_')[0]
                res = int(res.split('_')[1])
                if chain not in rec_chain_interacting_res:
                    rec_chain_interacting_res[chain] = []
                rec_chain_interacting_res[chain].append(res)
                
            #Now turn them into sorted lists
            lig_interactions = []
            for chain in sorted(list(lig_chain_interacting_res)):
                lig_interactions.append((chain, sorted(lig_chain_interacting_res[chain]))) 
            
            rec_interactions = []
            for chain in sorted(list(rec_chain_interacting_res)):
                rec_interactions.append((chain, sorted(rec_chain_interacting_res[chain]))) 
            
            
            #Build lists for each input for the process_RMSD function so that we can use
            #map with concurrent.futures.ProcessPoolExecutor
            rec1_list.append(rec1)
            lig1_list.append(lig1)
            rec2_list.append(rec2)
            lig2_list.append(lig2)
            rec_interface_list.append(rec_interactions)
            lig_interface_list.append(lig_interactions)
            job_list.append(count)
            
            #Count keeps track of the job_id
            count += 1
    with concurrent.futures.ProcessPoolExecutor() as executor:
        RMSDs = executor.map(process_RMSD, rec1_list, rec2_list, 
                             lig1_list, lig2_list, rec_interface_list, 
                             lig_interface_list, job_list)
    
    RMSDs = list(RMSDs)
    count = 0
    #Convert the RMSD list into a distance matrix
    for i in range(N):
        for j in range(i,N):
            
            dist = RMSDs[count]
            dists[i, j] = dist
            dists[j, i] = dist
            count += 1
           
    return(dists)
    
def make_result_json(model_list, dists, results_file_name):

    '''Save the distance matrix and the input model_list to a json for for use
    in clustering. 
    
    Args:
        model_list (list): List of dictionaries with an entry for each model. 
        dists (ndarray): Array of pairwise RMDSDs.
        results_file_name (str): Path to JSON where RMSD array and model information 
            is saved. 
        
    Returns:
        None
    
    '''
    
    model_names = [model['name'] for model in model_list]
    modified_dists = dists.tolist()
    results_dict = {'names':model_names, 'dists':modified_dists}
    with open(results_file_name, 'w') as f:
        json.dump(results_dict, f)
        
        