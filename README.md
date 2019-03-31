# SBI-project

A python package for macrocomplex modeling from protein interaction pairs.

## Introduction

This is a usage tutorial of the Complexbuilder command-line application. 

The objective of this project is to reconstruct a complete macrocomplex from a given set of interacting pairs of proteins. The package returns a PDB of CIF file with the possible models that were built.

Developers: Núria Olvera Ocaña and Helena Catena Sánchez, from the UPF MSc of Bioinformatics for Health Sciences.


## Installation

To install this package you must clone the it from GitHub:

    git clone MYLINK
    cd complexbuilder

And run the next comand:

    sudo python3 setup.py install


## How does it work?

### Input and output files

The input is a set of __PDB files__ holding pairs of proteins interacting. These inputs can be either redundant (all the possible interactions of the macrocomplex) or unique interactions (only those necesary to build the complex).

It is important to consider that the output files will depend on the command-line options and arguments the user determines while performing the analysis. 
The output will be saved in the directory defined by the user or in the current directory. 
Depending on the size of the macrocomplex, it will be saved in .pdb or .cif format.

### Python modules

* __complexbuilder:__ This is the _main_ module or program created to reconstruct a macrocomplex given interaction pairs of protein-protein or protein - RNA.
In addition, this module has the Argument Parser objext, from the argparse module from python. This module it's wused to give to the program comand-line options and arguments.  

* __functions.py:__ 

is module is composed by a set of different functions to solve biological and technical problems during the analysis. Thus, it is imported into the other modules in order to use the defined functions.
* __classes.py:__
* __utilities.py:__ 
