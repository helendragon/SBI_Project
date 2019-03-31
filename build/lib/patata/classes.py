
class IncorrectNumberofChains(ValueError):
    """Checks if the input structure has 2 chains. This program only takes pairs of interaction"""
    def __init__(self, n_chains):
        self.n_chains= n_chains        
    def __str__(self):
        return "You are providing %s chains and the program only accepts pairs of chains" %(self.n_chains)    

class IncorrectFileFormat(ValueError):
    """ Checks if the output extension provided by the user matches with what the program can do. If a complex has more than 52 chains, it will not be able to save it into a pdb file --> an error will be raised if the user tries to save it as pdb """
    def __init__(self, outputpath):
        self.outputpath= outputpath        
    def __str__(self):
        return "The format provided is not correct"

class NoPDBFiles(ValueError):
    """ Checks if there are pdb files in the input given by the user. If there are not this error will be raised. """
    def __init__(self, number_files):
        self.number_files= number_files        
    def __str__(self):
        return "No PDB files have been found in this directory"

class NoExistsResidue(ValueError):
    def __init__(self, residue):
        self.residue= residue        
    def __str__(self):
        return "%s is not found in the aminoacids dictionary" %(self.residue)
class NoExistsDesoxiribonucleotide(ValueError):
    def __init__(self, desoxi):
        self.desoxi= desoxi       
    def __str__(self):
        return "%s is not found in the desoxiribonucleotide dictionary" %(self.desoxi)
class NoExistsRibonucleotide(ValueError):
    def __init__(self, ribo):
        self.ribo= ribo       
    def __str__(self):
        return "%s is not found in the ribonucleotide dictionary" %(self.ribo)
