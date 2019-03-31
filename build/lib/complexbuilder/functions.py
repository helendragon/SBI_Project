import complexbuilder.utilities as utilities
import Bio.PDB
import copy
import os
import sys
import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 

from complexbuilder.classes import * 

def create_dict_pdb(directory):
     parser = Bio.PDB.PDBParser()
     all_files= os.listdir(directory)
     interaction_dict={}     
     files_pdb = [i for i in all_files if i.endswith('.pdb')]
     for pdb in files_pdb:
          count=0
          path1= directory
          struc= parser.get_structure(pdb,os.path.abspath(path1) + "/" + pdb)
          chain = struc.get_chains()    
          for element in chain:
               count +=1
          if count == 2:         
               interaction_dict[pdb]=struc
          else:
               raise IncorrectNumberofChains(count)
     return interaction_dict
    

def checking_files(input_file , parser):
    if input_file == None:
        curr_dir= os.getcwd()         
        number_files= (len([name for name in os.listdir(curr_dir) if name.endswith(".pdb")]))
        if number_files == 0:
            raise NoPDBFiles(number_files)
        if parser:
            s= "%s PDB files found\n" %(number_files)
            sys.stderr.write(s)
        return create_dict_pdb(curr_dir)                   
    elif os.path.isdir(input_file) == True:
        number_files= (len([name for name in os.listdir(input_file) if name.endswith(".pdb")]))
        if number_files == 0:
            raise NoPDBFiles(number_files)
        if parser:
            s= "%s PDB files found\n" %(number_files)
            sys.stderr.write(s)        
        return create_dict_pdb(input_file)


def get_fasta_dictionary(interaction_dict):
    """DOCU"""
    fasta_dict={}
    for name, structure in interaction_dict.items():
        chains =[]
        sequences = get_residues_sequence(structure)
        if len(sequences) == 2: #there are 2 chains
            A = sequences[0]
            B = sequences [1]
            for model in structure:
                for chain in model:
                    chains.append(chain.id)
            name_fastaA = name + "_" + chains[0]
            name_fastaB = name + "_" + chains[1]
            fasta_dict[name_fastaA] = A
            fasta_dict[name_fastaB] = B
        elif len(sequences) ==1: #there is only 1 chain there
            A = A = sequences[0]
            for model in structure:
                for chain in model:
                    chains.append(chain.id)
            name_fastaA = name + "_" + chains[0]
            fasta_dict[name_fastaA] = A
        else: #there are no chains
            continue
    return fasta_dict


def preparing(fasta_list, pdb_dict):
     for item1 in fasta_list:
          matchObj = re.search( '^(.*)_([a-zA-Z0-9])$', item1[0])
          fasta1= item1[1]
          if matchObj:
               original_name1= matchObj.group(1)
               original_structure1=pdb_dict[original_name1]
               chain_1= matchObj.group(2)               
               yield fasta1, [original_structure1, chain_1]    


def pairwise(fasta_list, pdb_dict, cutoff=0.98):
     interaction_setlist=[]
     for a in preparing(fasta_list, pdb_dict):
          for b in preparing(fasta_list, pdb_dict):
               merge= a[1] + b[1]
               setito = list(merge)
               if setito not in interaction_setlist:
                    alignments=pairwise2.align.globalxx(a[0], b[0], score_only=1)                      
                    len_max = min(len(a[0]), len(b[0]))
                    cutoff_chain = alignments / len_max
                    if cutoff_chain > cutoff:
                         interaction_setlist.append(setito)
     return interaction_setlist



def get_CA_atoms (residue_list):
    ''' Docu'''
    list_atoms=[]
    for residue in residue_list:
        atoms_minilist = list(residue.get_atoms())
        for atom in atoms_minilist:
            if atom.id =="CA" or atom.id == 'P':
                list_atoms.append(atom)
    return list_atoms


def get_common_residues (chainA, chainB):
    '''DOCU'''
   
    chainA_residues = list(chainA.get_residues())
    chainB_residues = list(chainB.get_residues())
    
    chainB_residues_ids=[x.id[1] for x in chainB_residues if ('CA' in x or 'P' in x)]

    same_res= [r for r in chainA_residues if(r.id[1] in chainB_residues_ids and ('CA' in r or 'P' in r))]

    return same_res


def get_residues_sequence (structure):
    """Function that takes a structure and returns a LIST of strings with the fasta sequences of AA"""
    list_fastas=[]
    residues_dict = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
         'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
         'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    rna_dict={ '  G': 'G', '  C': 'C', '  A': 'A', '  T': 'T', '  U': 'U'}
    dna_dict={ ' DG': 'G', ' DC': 'C', ' DA': 'A', ' DT': 'T', ' DU': 'U'}

    for model in structure:
        for chain in model:
            string=""
            for residue in chain:
                for atom in residue:
                    if atom.id == "CA":
                        if residue.resname == 'UNK':
                            continue                        
                        else:
                            residue_name = residues_dict[residue.resname]
                            string = string + residue_name
                    elif atom.id == "P":
                        if residue.resname.startswith(" D"):
                            residue_name= dna_dict[residue.resname]
                        else:
                            residue_name= rna_dict[residue.resname]
                        string= string + residue_name
            if string != "":
                list_fastas.append(string)
    return list_fastas
   





def get_clashes (fixed_struc, moving_struc, minimum_clash_distance):
    """DOCU"""
    clash_set=set()
    
    NS = Bio.PDB.NeighborSearch(list(fixed_struc.get_atoms()))
    clashes =0
    for atoms in moving_struc.get_atoms():
        close = NS.search(atoms.get_coord(), minimum_clash_distance)
        if len(close) > 0:
            for item in close:
                clash_set.add(item.get_parent().get_parent().id)
                clashes +=1
    return clashes        


def build_complex (seed, interaction_setlist1, number_of_chlashes, clash_distance, unique, parser):
    """DOCU"""
    if not unique:
        interaction_setlist=copy.deepcopy(interaction_setlist1)
    else:
        interaction_setlist = interaction_setlist1
    superimposer = Bio.PDB.Superimposer()
    seed_copy = copy.deepcopy(seed)
    complex1 = seed_copy
    c=0
    added_chains = {}
    for chains in complex1[0]:
        added_chains[chains.id] = get_parent_struc (chains, added_chains)
    for fixed_chain in complex1[0]: 
        #We prepare the fixed chain parameters to superimpose
        fixed_id = get_chain_name(fixed_chain)
        structure_fixed= get_parent_struc (fixed_chain, added_chains)
        if parser:
            print(structure_fixed.id + "fixed")
        fixed = list([structure_fixed, fixed_id])
        if fixed_id not in added_chains.keys():
            added_chains[fixed_id] =  structure_fixed
        for interaction in interaction_setlist:
            #check which interaction the fixed chain has
                if interaction[0:2] == fixed: 
                    #we prepare the moving chain to superimpose
                    moving = interaction[2:4]
                    moving_struc = moving[0]
                    id_chainM = moving[1]
                    if parser:
                        print(moving_struc.id + "moving")

                    chain_moving = moving_struc[0][id_chainM]
                    #get the common atoms from both chains
                    common_1and2 = get_common_residues(fixed_chain,chain_moving)
                    common_2and1 = get_common_residues(chain_moving, fixed_chain)
                    #get the CA from both chains
                    atoms_1and2 = get_CA_atoms(common_1and2)
                    atoms_2and1 = get_CA_atoms(common_2and1) 
                    
                    #We make sure that the lengths are the same and not 0 before superimpsoe
                    if len(atoms_1and2) == len(atoms_2and1) and len(common_1and2) >0 and len(atoms_1and2) != 0:
                        superimposer.set_atoms(atoms_1and2, atoms_2and1)
                        #we copy the moving structure to detach a chain without modifying the real one
                        moving_struc_copy = copy.deepcopy(moving_struc)
                        #We detach the moving structure to check the clashes between the complex and the other structure
                        for model in moving_struc_copy:
                            model.detach_child(id_chainM)
                        superimposer.apply(list(moving_struc_copy.get_atoms()))
                        #We check the clashes 
                        if get_clashes(complex1, moving_struc_copy, clash_distance) < number_of_chlashes:
                            changes = change_names(complex1, moving_struc_copy, c)
                            complex1 = changes[0]
                            moving_struc_copy = changes[1]
                            c = changes[2]
                            chainAdded= get_chains(moving_struc_copy)
                            complex1[0].add(moving_struc_copy[0][chainAdded[0].id])
                            added_chains[chainAdded[0].id] = moving_struc
                            if parser:
                                print ("Chain added")
                            if not unique:
                                for i,x in enumerate(interaction_setlist):                             
                                    if len(x) != 0: #ta mal 
                                        if x[2:4] == moving:
                                            interaction_setlist[i] = []
                                        if x == [moving, fixed]:
                                            interaction_setlist[i] = []  
                        else:
                            if parser:
                                print("Clashes")
                    else:
                        continue
    return complex1

                        

def get_chain_name (chain):
    """DOCU"""
    if "-" in chain.id:
        id_chain=chain.id[-1]
    else:
        id_chain=chain.id
    return id_chain

def get_parent_struc (chain, added_chains_dict):
    """DOCU"""
    if "-" in chain.id:

        structure = added_chains_dict[chain.id]
    else:
        structure = chain.get_parent().get_parent()
    return structure

def change_names (fixed_structure, moving_structure, index):
    """Function that changes the name of the new chains to be added to the complex"""
    for chain in moving_structure[0]:
        chain.id = utilities.merged_list[index]+"-"+chain.id
        index +=1
    return (fixed_structure, moving_structure, index)

def get_chains (structure):
    """DOCU"""
    chains=[]
    for chain in structure[0]:
        chains.append(chain)
    return chains

def save_complex (complex, outputpath ,count):
    """DOCU"""
    c=0
    for chain in complex[0]:
        chain.id = utilities.merged_list[c]
        c+=1
    if os.path.isdir(outputpath): # The file is saved in this path with default name
        if len(complex[0]) < 52:
            format1= ".pdb"
            io = Bio.PDB.PDBIO()          
        else:
            format1= ".cif"
            io = Bio.PDB.MMCIFIO()
        io.set_structure(complex)
        io.save(outputpath + "Complex"+ "_" + str(count) + format1)  
    else: #if the output is just the file name, it is saved in the current directory 
        if len(complex[0]) < 52:
            io = Bio.PDB.PDBIO()
            if not ".pdb" or not ".cif" in outputpath: #the output name does not have the correct format
                format1 = ".pdb"
            elif outputpath.endswith (".cif"):
                raise IncorrectFileFormat(outputpath) 
            else:
                format1= ".pdb"
                outputpath= outputpath[:-4]
                             
        else:
            if not ".pdb" or not ".cif" in outputpath:
                format1 = ".cif"
            elif outputpath.endswith (".pdb"):
                raise IncorrectFileFormat(outputpath)
            else:
                format1= ".cif"
                ouputpath= outputpath[:-4]                            
            io = Bio.PDB.MMCIFIO()
        io.set_structure(complex)
        io.save(outputpath + "_" + str(count) + format1)

    return "Model Saved"


    

