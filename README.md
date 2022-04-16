# Interface_clustering
This repository contains the python scripts that can be used to cluster models
in a directory based on the interfaces between the receptor chains and the ligand chains.  

BLAST aligned models are stored in the input directory in a BLAST_aligned directory.

The RMSDs and clusters are stored in RMSDs.json and clusters.json files respectively in 
the input directory. Clusters.json contains a dictionary in which each cluster center 
is mapped to a list of the members in its cluster. 

Command: 
python cluster_models.py -i {directory containing the models} -rc {receptor chains} -lc {ligand chains} -it {interface threshold in Angstroms} -ct {clustering threshold in Angstroms}

Interface threshold defines the distance under which residue atoms on the receptor and ligand are considered to be in contact. These contacts determine the 
interface. 

Clustering threshold determines the RMSD threshold below which models are considered to be similar. 

Command for provided sample data:

#Low clustering threshold results in smaller, more similar clusters.
python cluster_models.py -i ./sample_models_1 -rc A -lc BC -it 10.0 -ct 0.5 

#High clustering threshold results in bigger, more diverse clusters.
python cluster_models.py -i ./sample_models_2 -rc A -lc B -it 10.0 -ct 5.0


Requirements:  
-Blast  
-Numpy  
-Bio  
-Pymol  
-Prody  
-Scipy  