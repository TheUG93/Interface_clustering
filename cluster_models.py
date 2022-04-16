'''
Author : Usman Ghani

This script uses the scripts in the interface_clustering module
to cluster the models in the provided directory.
'''
from interface_clustering import blast_alignment
from interface_clustering import clustering
import argparse
import shutil
import glob
import sys
import os


def get_models_and_align_them(input_dir):
    
    '''
    Collect all the pdb files in the directory, make a new directory for aligned models,
    and, create aligned models in new directory.
    
    Args:
        input_dir (str): Path to the directory containing the pdbs of interest.
    
    Returns:
        new_dir (str): Path where BLAST aligned versions of the models are stored. 
    
    '''

    models = glob.glob(os.path.join(input_dir, '*.pdb'))
    
    target_model = models[0] #Choosing the first model as the one that all other models should align to
    
    #Make new directory for the aligned models
    new_dir = os.path.join(input_dir, 'BLAST_aligned_models')
    
    if not os.path.exists(new_dir):
        os.mkdir(new_dir)
        
    #targte model does not need to be aligned to itself so 
    #it is directly copied to the aligned model directory
    shutil.copy(target_model, new_dir)
    
    for model in models[1:]:
        blast_alignment.align(model, target_model, 
            os.path.join(new_dir, os.path.basename(model)))
        
    return(new_dir)
    

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', 
        '--input_dir', 
        required=True, 
        help='The directory with all the models that need to be clustered'
        )
    parser.add_argument(
        '-rc', 
        '--rec_chains', 
        required=True, 
        help='Receptor chains'
        )
    parser.add_argument(
        '-lc', 
        '--lig_chains', 
        required=True, 
        help='Ligand chains'
        )
    parser.add_argument(
        '-it', 
        '--interface_threshold', 
        required=True, 
        help='Cutoff distance between receptor and ligand in Angstroms for interface conisderation'
        )
    parser.add_argument(
        '-ct', 
        '--clustering_threshold', 
        required=True, 
        help='RMSD cutoff for models to be considered similar'
        )
    args = parser.parse_args() 
    
    try:
        print('Aligning models using BLAST...')
        BLAST_aligned_models_dir = get_models_and_align_them(args.input_dir)
    except:
        print('Alignment failed.')
        sys.exit(1)
    
    RMSD_json = os.path.join(args.input_dir, 'RMSDs.json')
    cluster_json = os.path.join(args.input_dir, 'clusters.json')
    clustering.cluster(
        BLAST_aligned_models_dir, 
        args.rec_chains, 
        args.lig_chains, 
        float(args.interface_threshold), float(args.clustering_threshold), RMSD_json, cluster_json)
    
    print('Clustering succeeded!')

if __name__ == '__main__':

    main()



