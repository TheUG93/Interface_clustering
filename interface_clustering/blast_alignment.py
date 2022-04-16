'''
Author : Usman Ghani

This script is used to perform a BLAST sequence alignment on 
all models. It is important for all models to have the same
residue numbering for interface RMSD calculations. 
'''
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from io import StringIO
from Bio import SeqIO
import Bio
import os

def get_sequences(model_file):

    '''Get the sequence and residue numbers for
    each chain in the pdb provided.
    
    Args : 
        model_file: The path to the pdb file of interest.
    
    Returns:
        sequences (dict of str: str): Maps chain names to sequences.
        aligned_to_seq_res_nums_dict (dict of str: list): Maps chain names to a 
            list residue of numbers in the chain. The position of the residue
            number in the list matches the position of the seqeunce.
    '''

    chain_res_numbers = {}
    with open(model_file, 'r') as f:
        
        for line in f:
            
            if line.startswith('ATOM'):
                
                chain = line[21]
                
                if chain not in chain_res_numbers:
                    chain_res_numbers[chain] = {}
                    
                res_num = int(line[22:26])
                
                if res_num in chain_res_numbers[chain]:
                    continue
                
                res_name = line[17:20].strip()
                
                if res_name.lower().capitalize() not in protein_letters_3to1:
                    continue
               
                res_code = protein_letters_3to1[res_name.lower().capitalize()]
                
                chain_res_numbers[chain][res_num] = res_code
    
    
    res_nums = {}
    for chain in chain_res_numbers:
        res_nums[chain] = sorted(list(chain_res_numbers[chain]))
        
        
    sequences = {}
    aligned_to_seq_res_nums_dict = {}    
    for chain, res_num_list in res_nums.items():
        
        sequence = ''   
        aligned_to_seq_res_numbers = [] 
        
        for i, res_num in enumerate(res_num_list):
            
            code = chain_res_numbers[chain][res_num] 
            if i > 0:
                prev_res_num = res_num_list[i-1]
                
                #Check if there is a gap in the residue numbers and 
                #add a "-" to mark it
                if res_num - prev_res_num > 1:
                    code = '-' + code
                    
            sequence += code
            aligned_to_seq_res_numbers.append(res_num)
            
        sequences[chain] = sequence
        aligned_to_seq_res_nums_dict[chain] = aligned_to_seq_res_numbers
        
    return(sequences, aligned_to_seq_res_nums_dict)

def get_model_to_native_mapping(native_lig_sequence, native_aligned_res_numbers, 
                                model_lig_sequence, model_aligned_res_numbers, chain):

    '''Map the model sequence to the native sequence. 
    
    It uses blast to align the sequences, and then builds a dictionary where 
    for each model residue number, it stores the value for the corresponding 
    native residue number.
    
    Args:
        native_lig_sequence (str): Sequence of the native model.
        native_aligned_res_numbers (dict of str: list): list of residue numbers 
            in the native chain aligned to the native chain sequence. 
        model_lig_sequence (str): Sequence of the model.
        model_aligned_res_numbers: list of residue numbers 
            in the model chain aligned to the model chain sequence. 
        chain (str): Name of chain of interest.
        
    Returns:
        actual_model_to_native_mapping (dict of int:int): Residue numbers in the 
            model mapped to the residue numbers in the native.
            
    '''
     
    #Convert sequences into seq objects
    query_seq_obj = SeqRecord(Bio.Seq.Seq(model_lig_sequence), id='query')
    native_seq_obj = SeqRecord(Bio.Seq.Seq(native_lig_sequence), id='native_{}'.format(chain))
    
    #Write seq objects to fasta files
    SeqIO.write(query_seq_obj, 'query.fasta', 'fasta')
    SeqIO.write(native_seq_obj, 'native_{}.fasta'.format(chain), 'fasta')
    
    #Blastp-short has to be used if the sequence is small. Normal Blastp does not produce an alignment for 
    #short sequences
    if len(native_lig_sequence) <= 25:
        output = NcbiblastpCommandline(query="query.fasta", subject="native_{}.fasta".format(chain), 
                                       task='blastp-short', outfmt=5)()[0]
    else:
        output = NcbiblastpCommandline(query="query.fasta", subject="native_{}.fasta".format(chain), 
                                       outfmt=5)()[0]
    
    #Get rid of the fasta files
    os.remove('query.fasta')
    os.remove(f'native_{chain}.fasta')
    
    blast_result_record = NCBIXML.read(StringIO(output))
    
    #Get the highet ranked alignment
    rel_alignment = blast_result_record.alignments[0]
    
    hsp = rel_alignment.hsps[0]
    query_start = hsp.query_start
    native_start = hsp.sbjct_start
    matched_query_seq = hsp.query
    matched_native_seq = hsp.sbjct
    
    counter = query_start - 1 #BLAST starts counting from 1
    query_indxs = []
    for res in matched_query_seq:
        
        if res == '-':
            query_indxs.append(-1)
        else:
            query_indxs.append(counter)
            counter += 1
    
    counter = native_start - 1
    native_indxs = []
    for res in matched_native_seq:
        
        if res == '-':
            native_indxs.append(-1)
        else:
            native_indxs.append(counter)
            counter += 1
    
    #Map the residue numbers for the query to the native. These
    #residue numbers start from 0 as they are relative to their 
    #position in the sequence; not the same as the residue numbers
    #assigned in the model file         
    mapping = {}
    for indx1, indx2 in zip(query_indxs, native_indxs):
        
        if indx1 == -1 or indx2 == -1:
            continue
            
        mapping[indx1] = indx2
    
    #Tranlate the sequence relative indexes to      
    actual_model_to_native_mapping = {}
    for indx1 in mapping:
        
        actual_indx1 = model_aligned_res_numbers[indx1]
        actual_indx2 = native_aligned_res_numbers[mapping[indx1]]
        
        actual_model_to_native_mapping[actual_indx1] = actual_indx2
        
    return(actual_model_to_native_mapping)

def align_models(model1_path, model2_path):

    '''Map residues on the chains in model1 to the relevant residues on 
    the same chains in model2. 
    
    Only keeps chains that are present in both models.
    
    Args:
        model1_path (str): Path to model 1.
        model2_path (str): Path to model 2.
        
    Returns:
        model1_to_model2_mapping (dict of int:int): Residue numbers in model 1
            mapped to residue numbers in model 2.
    
    '''

    model1_sequence, model1_aligned_to_seq_res_nums_dict = get_sequences(model1_path)
    model2_sequence, model2_aligned_to_seq_res_nums_dict = get_sequences(model2_path)
    
    model1_to_model2_mapping = {}
    for chain in (set(model1_sequence) & set(model2_sequence)):
        mapping = get_model_to_native_mapping(model2_sequence[chain], 
                                              model2_aligned_to_seq_res_nums_dict[chain], 
                                              model1_sequence[chain], 
                                              model1_aligned_to_seq_res_nums_dict[chain], 
                                              chain)
        model1_to_model2_mapping[chain] = mapping
        
    assert len(model1_to_model2_mapping) > 0, f'No common chains found for {model1_path} and {model2_path}. Models can not be aligned' 
        
    return(model1_to_model2_mapping)
    
    
def make_new_model(model1_path, model1_to_model2_mapping, output_path):

    '''Make a new model1 file in which the residue numbers match those in
    model2.
    
    Args:
        model1_path (str): Path to model 1.
        model1_to_model2_mapping (dict of int:int): Residue numbers in model 1
            mapped to residue numbers in model 2.
        output_path (str): Path where a new model 1 is saved in which the residue
            numbers correspond to those in model 2.  
    Returns:
        None
    
        
    '''

    new_file = []
    with open(model1_path, 'r') as f:
        
        for line in f:
            
            if not line.startswith('ATOM'):
                new_file.append(line)
                continue
                
            chain = line[21]
                
            res_number = int(line[22:26])
            
            #Ignore residues that were not mapped
            if res_number not in model1_to_model2_mapping[chain]:
                continue
            
            new_res_number = model1_to_model2_mapping[chain][res_number]
            new_res_number = str(new_res_number).rjust(4)
            
            new_line = line[:22] + new_res_number + line[26:]
            new_file.append(new_line)
        
    with open(output_path, 'w') as f:
        f.writelines(new_file)
        

def align(model1_path, model2_path, new_model1_path):
    
    '''Align model1 sequences to model2 sequences and create a new model1 where
    the sequences are aligned.
    
    Args:
        model1_path (str): Path to model 1.
        model2_path (str): Path to model 2.
        new_model1_path (str): Path where a new model 1 is saved in which the residue
            numbers correspond to those in model 2.  

    Returns:
        None
    
    '''
    
    model1_to_model2_mapping = align_models(model1_path, model2_path)
    make_new_model(model1_path, model1_to_model2_mapping, new_model1_path)
    
    
    
    