import patata.utilities as utilities
import Bio.PDB
import copy
import os
import sys
import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment 

from patata.classes import * 

def create_dict_pdb(directory):
    """ This functions looks into the directory given by the user, all the files that end with .pdb; then counts how many chains are in each one. If the file does not contain exactly 2 chains an error is raised. This function returns a dictionary with the structure as a value and the key is the pdb file name."""
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
    """ This function gets the input of the user. If there argument "-i" is None, then it looks in the current directory for files that end with .pdb. If a directory is given and no PDB files are found, an error is raised. The function calls the create_dict function that will create a dictionary with the PDB files."""
    if input_file == None:
        curr_dir= os.getcwd()         
        number_files= (len([name for name in os.listdir(curr_dir) if name.endswith(".pdb")]))
        if number_files == 0:
            raise NoPDBFiles(number_files)
        if parser:
            s= "%s PDB files found\n" %(number_files)
            sys.stdout.write(s)
        return create_dict_pdb(curr_dir)                   
    elif os.path.isdir(input_file) == True:
        number_files= (len([name for name in os.listdir(input_file) if name.endswith(".pdb")]))
        if number_files == 0:
            raise NoPDBFiles(number_files)
        if parser:
            s= "%s PDB files found\n" %(number_files)
            sys.stdout.write(s)        
        return create_dict_pdb(input_file)


def get_fasta_dictionary(interaction_dict):
    """This function gets the name and structure of the pdb files given by the user saved in the dictionary created by the create_pdb_dict function. This functions calls the get_residues_sequence to obtain a list with the sequences of the chains of the structure. If the structure has 2 chains it will save the name of the file + the id of the chain  and its corresponding seqeunce in a dictionary. It will do the same if the structure only has 1 chain. However, it will do nothing if no chains are found in a structure. This function returns a dictionary where the key is the name of the structure + the chain ID and the value is the sequence of the chain. """
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
    """ This is a generatior function that takes the a fasta list with the contents of the fasta dictionary as (key, value) tuples and the pdb dictionary with the structure as value and the name of the sequence as key.
    The function first looks at the first item of the tuple (the name of the pdb file + the chain that it has the sequence) and extracts the name of the file and the sequence name. This function returns a tuple with the sequence and a list with the structure object of the chain and the chain id."""
    for item1 in fasta_list:
        matchObj = re.search( '^(.*)_([a-zA-Z0-9])$', item1[0])
        fasta1= item1[1]
        if matchObj:
            original_name1= matchObj.group(1)
            original_structure1=pdb_dict[original_name1]
            chain_1= matchObj.group(2)               
            yield fasta1, [original_structure1, chain_1]    


def pairwise(fasta_list, pdb_dict, cutoff=0.98):
    """ This function takes the fasta list of tuples of (structure+chain name, sequence) and the dictionary with the structure object as value and the pdb file name as key. This function iterates through the generator created by the preparing function and aligns each of the chains. If the alignment has a score higher than the cuttof, it is added to a list of interactions. This list of interactions is a list of list, each interaction consists of a list that has the structure objects of the 2 stuctures aligned and the ID of the chain that was aligned like so: [struc1, chain1, sctruc2, chain2]. This function returns the list of interactions."""
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
    ''' This functions returns a list of CA or P atoms from a list of residues. First is gets the atoms from the resdiue object list and looks for CA or P ids. If it finds one, it appends it to the list of atoms. '''
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
    """Function that takes a structure and returns a LIST of strings with the fasta sequences of AA or nucleotides (DNA or RNA) depending on the residuename of the CA or P atom. If the functions finds a residue that is UNK, it will skip it. If a residue is not found in any of the dictionaries (DNA, RNA or protein) it will raise an error depending on where it is not found."""
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
                            if residue.resname not in residues_dict:
                                raise NoExistsResidue(residue.resname)                            
                            residue_name = residues_dict[residue.resname]                                
                            string = string + residue_name
                    elif atom.id == "P":
                        if residue.resname.startswith(" D"):
                            residue_name= dna_dict[residue.resname]
                            if residue.resname not in dna_dict:
                                raise NoExistsDesoxiribonucleotide(residue.resname)
                        else:
                            residue_name= rna_dict[residue.resname]
                            if residue.resname not in rna_dict:
                                raise NoExistsRibonucleotide(residue.resname)
                        string= string + residue_name
            if string != "":
                list_fastas.append(string)
    return list_fastas

   





def get_clashes (fixed_struc, moving_struc, minimum_clash_distance):
    """This function checks the clashes between the fixed structure (the complex) and the moving structure before being added. The minimum clash distance is defined by the user, but the default is 1.5 A. It returns the number of steric clashes found"""
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


def build_complex (seed, interaction_setlist1, number_of_chlashes, clash_distance, unique, parser, n_max_chains):
    """This functions takes as arguments:
    - seed: structure that has more interactions and it's used to start the complex.
    - interaction_setlist: list of interactions created by the pairwise function. 
    - number_of_clashes: minimum number of clashes between 2 structures, this can be changed by the user, the default value is 5.
    - clash_distance: this is the minimum clash distance used by the get_clashes function. The default value is 1.5 A, but can be changed by the user.
    - unique: True or False value that indicates that the input provided by the user contains non redundant interactions. By deafult this is False. 
    - parser: True or False value, depends on the -verbose statement. 

    This function iterates through all the chains of the complex and finds on the interaction list the chains that interact with it. Once a chain that interacts with the complex chain is found, the common residues are found and they are superimposed. Once they are superimposed, the steric clashes are found, and if there are found less than the minimum given, the structure is added to the complex. 

    Once the structue is added, the name of the chain is changed (to aboid problems with chains comming from different structures). The name of the chain is saved into a dictionary, where the parent structure is saved as well.
    If the unique function is not activated, the interactions that add the chain that has just been added are removed and the interactions that are the oposite of what was just added.

    This function returns the complex
    """
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
        if len(complex1[0]) < n_max_chains:
            fixed_id = get_chain_name(fixed_chain)
            structure_fixed= get_parent_struc (fixed_chain, added_chains)
            if parser:
                sys.stdout.write("The structure " +structure_fixed.id + " is fixed\n")
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
                            sys.stdout.write("The movile structure is "+moving_struc.id + "\n")

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
                                    sys.stdout.write ("The Chain "+ chainAdded[0].id + "has been added to the complex\n")
                                if not unique:
                                    for i,x in enumerate(interaction_setlist):                             
                                        if len(x) != 0: #ta mal 
                                            if x[2:4] == moving:
                                                interaction_setlist[i] = []
                                            if x == [moving, fixed]:
                                                interaction_setlist[i] = []  
                            else:
                                if parser:
                                    sys.stdout.write("Too many clashes were found, going to the next structure\n")
                        else:
                            continue
    return complex1

                        

def get_chain_name (chain):
    """Returns the ID of the chain object. It does not matter if the id has been changes after adding the chain the to complex. It returns the "real" id of the chain object """
    if "-" in chain.id:
        id_chain=chain.id[-1]
    else:
        id_chain=chain.id
    return id_chain

def get_parent_struc (chain, added_chains_dict):
    """Function that gets the parent structure from a chain. It needs as an input the chain object and the dictionary of added chains to the complex.  
    The function checks if the ID of the chain has been modified (has a "-") and if it finds one it looks for this id in the added dictionary.
    This function returns the parent structure object  """
    if "-" in chain.id:
        structure = added_chains_dict[chain.id]
    else:
        structure = chain.get_parent().get_parent()
    return structure

def change_names (fixed_structure, moving_structure, index):
    """Function that changes the name of the new chains to be added to the complex, it adds a letter from the alphabet before the real ID separated by a "-" """
    for chain in moving_structure[0]:
        chain.id = utilities.merged_list[index]+"-"+chain.id
        index +=1
    return (fixed_structure, moving_structure, index)

def get_chains (structure):
    """This functions returns a list with chain objects of a given structure"""
    chains=[]
    for chain in structure[0]:
        chains.append(chain)
    return chains

def save_complex (complex, outputpath ,count):
    """
    This function takes the complex and saves it into a pdb or cif file (depending on how many chains it has). One of the inputs is the output path - the function checks if it's a directory and if the name of the complex is given (it will be saved in the current directory). The count is necessary to help the user keep track of the models that are saved. 
    """
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
                outputpath= outputpath[:-4]                            
            io = Bio.PDB.MMCIFIO()
        io.set_structure(complex)
        io.save(outputpath + "_" + str(count) + format1)

    return "Model Saved"


    

