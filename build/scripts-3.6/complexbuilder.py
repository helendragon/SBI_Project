import Bio.PDB

from Bio import BiopythonWarning
import copy
import sys
import os
import re

import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

import complexbuilder.functions 



from Bio import pairwise2
from Bio.pairwise2 import format_alignment 


import argparse

parser = argparse.ArgumentParser(description="This program does BLA BLA BLA")

parser.add_argument('-i', '--input',
    dest = "infile",
    action = "store",
    default = None,
    help = "Input FASTA formatted file")
parser.add_argument('-o', '--output',
    dest = "outfile",
    action = "store",
    default = "Complex",
    help = "outputfile")
parser.add_argument('-c', '--clashes',
    dest = "minimum_clash_distance",
    action = "store",
    default = 1.5,
    type = int,
    help = "Integer specifying")
parser.add_argument('-n', '--number',
    dest= "number_of_clashes",
    action= "store",
    default= 5,
    type= int,
    help= "Integer specifying")
parser.add_argument( '-u', '--unique',
    dest = "interactions",
    action= "store",
    default= False,
    help= "interact")
parser.add_argument('-v', '--verbose',
    dest = "verbose",
    action = "store_true",
    default = False,
    help = "Print log in stderr")

options = parser.parse_args()
options.infile
options.verbose



#The path given in the input is a path of a directory
#input_dir= sys.argv[1]

#input_miau="/home/helena/Documentos/Master/SBI/Project/hemoglobin_fail"

pdb_dict= functions.checking_files(options.infile)
#pdb_dict= functions.checking_files(input_miau)



#tenim un dicionari amb el nom del pdb de key i la estructura com a value

fasta_dict = functions.get_fasta_dictionary(pdb_dict)


sd = sorted(fasta_dict.items())
fasta_list=[]
for k,v in sd:
    tupla= (k, v)
    fasta_list.append(tupla)

interaction_setlist = functions.pairwise(fasta_list, pdb_dict)


# Ara tenim una llista de tuplas amb cadenes que son iguales. 



Dict={}
# Now we superimpose
for element in interaction_setlist:
     Id= element[0].id
     ch= element[1]
     key= Id
     if key not in Dict:
          Dict[key] = 1
     else:
          Dict[key] += 1

superimposer = Bio.PDB.Superimposer()

# for the seed we choose the structures with more interactions - pyt it in a list

#max_interactions = max(Dict, key=Dict.get)
maxValue = max(Dict.values())
possible_seeds = ([key for key in Dict.keys() if Dict[key]==maxValue])

#nameChain = max_interactions[-1]
#nameFile = max_interactions [:-1]
seeds = []
for possible in possible_seeds:
     for item in interaction_setlist:
          if item[0].id == possible:
               seeds.append(item[0])
               break

outputpath= options.outfile
count = 1
maximum_models=5
if options.verbose:
     print (str(len(seeds)) + " models can be generated, your maximum is " + str(maximum_models))

model_list=[]
for seed in seeds:
     while (count <= maximum_models):
          complex=functions.build_complex(seed, interaction_setlist, options.number_of_clashes, options.minimum_clash_distance, options.interactions)
          model_list.append(complex)
          functions.save_complex(complex, outputpath, count)
          if options.verbose:
               print ("Model " + str(count) + " generated")
          count += 1
          break




               