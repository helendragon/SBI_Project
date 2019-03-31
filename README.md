# SBI-project

A python package for macrocomplex modeling from protein interaction pairs.

## Introduction

This is a usage tutorial of the Complexbuilder command-line application.

The objective of this project is to reconstruct a complete macrocomplex from a given set of interacting pairs of proteins. The package returns a PDB of CIF file with the possible models that were built.

Developers: Núria Olvera Ocaña and Helena Catena Sánchez, from the UPF MSc of Bioinformatics for Health Sciences.


## Installation

To install this package you must clone the it from GitHub:

    git clone https://github.com/helendragon/SBI_Project.git

    cd complexbuilder

And run the next command:

    sudo python3 setup.py install

Now the package can be ran from any part of your system, you don’t have to be in the complexbuilder folder to use it.

## How does it work?

### Input and Output Files

The input is a set of __PDB files__ holding pairs of proteins interacting. These inputs can be either redundant (all the possible interactions of the macrocomplex) or unique interactions (only those necessary to build the complex). <br />
This program __does not__ accept compressed files.

It is important to consider that the output files will depend on the command-line options and arguments the user determines while performing the analysis.
The output will be saved in the directory defined by the user or in the current directory. <br />
Depending on the size of the macrocomplex, it will be saved in .pdb or .cif format.

### Python Modules

* __complexbuilder:__ This is the _main_ module or program created to reconstruct a macrocomplex given interaction pairs of protein-protein or protein - RNA.
In addition, this module has the Argument Parser object, from the argparse module from python. This module it's used to give to the program command-line options and arguments. 

* __functions.py:__ This module is formed by the different __functions__ used to go from the set of pairs of interactions to the final macrocomplex. This module is imported to the _main_ module so that the functions can be used.

* __classes.py:__ This module is composed by the different __class__ errors that can be raised during the execution of this package.

* __utilities.py:__ This module is composed by different variables needed for the functions and __main__ module during the analysis.

### Command-line Arguments

* __-i --input__: Input a directory containing the pdb files with the interactions pairs that the user wants to process. If no input is given, the package will try to look for the files in the current directory.

* __-o --output__: Output directory where the macrocomplexes will be saved as PDB files or CIF files. This can also be the name that the user wants the complexes to have. If no output is defined, it will be saved by the default name in the current directory.

* __-c --clashes__: Minimum clash distance between 2 atoms. The default minimum is 1.5 A.

* __-n --number__: Maximum number of clahses permited between 2 structures. The default maximum is 5.

* __-C --chains__: Maximum number of chains (iterations) the user wants to add to the complex. This is useful when only one interaction pair is given and it could be added infinite times. The default value is 100.

* __-m --models__: Maximum number of files (complexes) the user wants to obtain. If this is not specified, the maximum number of complexed will be 5.

* __-u --unique__: Gives the user the option to tell the program that the interactions that are given are not redundant. By default this is set to false.

* __-v --verbose__: indicates if the progression log has to be printed. By default this is et to false.

### Algorithm of the Package

The problem that is presented to us is to reconstruct a macrocomplex using the interaction pairs given by the user. <br /> <br />
The way that we have aproached this issue is by first of all identifying those sequences that are the same (that have a pairwise sequence identity of > 98%). We assumed that, if they have a very similar sequence, there are a lot of chances that they are the same chain, therefore they can be superimposed. <br /><br />
Next, in order to start building the complex, we thought that is was important that the first sequence to be added in the complex was one that had the maximum number of interactions. By doing this we would aboid starting with structures that would have very limited interactions and could not create the whole complex from them. Since it is possible that a many structures have the maximum number of interactions, all of them are used as the first structure (seed structure) to be added to the complex (by default the maximum models is 5). <br /><br />
The next step is the construction of the macrocomplex. First we look for all the interactions that a given chain has from the macrocomplex has, and we try to superimpose the structures that interact with it. Once we superimpsoe, we apply the rotation matrix to the chain that has to be added and the clashes between the atoms are checked. If clashes are found, the chain is not added and the program passes to the next interaction. If no clashes are found the structure is added to the complex and the programm passes to the next interaction. <br />
The program will continue to check if chains can be added until there are no more chains to look for in the complex. In other words, all the chains in the complex have been checked and no more structures can be added because there are clahses. <br />
This process is repeated for all the seed structures. <br /><br />
Finally, the models are saved into PDB or CIF format, depending on how many chains the final complex has. 

### System Requirements

In order to run this package with all its functionalities the user must have several programs:
* Python 3
* Chimera or Pymol
* Python modules:
    * biopython
    * copy
    * argparse
    * re
    * sys
    * os



