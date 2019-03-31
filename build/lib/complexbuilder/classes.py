
class IncorrectNumberofChains(ValueError):
    def __init__(self, n_chains):
        self.n_chains= n_chains        
    def __str__(self):
        return "You are providing %s chains and the program only accepts pairs of chains" %(self.n_chains)    

class IncorrectFileFormat(ValueError):
    def __init__(self, outputpath):
        self.outputpath= outputpath        
    def __str__(self):
        return "The format provided is not correct"

class NoPDBFiles(ValueError):
    def __init__(self, number_files):
        self.number_files= number_files        
    def __str__(self):
        return "No PDB files have been found in this directory"

