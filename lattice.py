import numpy as np
class Lattice:
    """Lattice class
    This class defines lattice, sublattices, and their hoppings. 
    This class returns primitive unit cell (or supercell)        
    """
    
    def __init__ (self, a1, a2=None, a3=None):
        """ 
        Construct a lattice object
        
        Parameters
        -------------------
        a1, a2, a3 : list of float 
            lattice vectors
        -------------------
        Examples:
        2D square lattice
        a1 = [1, 0, 0]
        a2 = [0, 1, 0]
        """
        
        
        if a1 != None:
            self.a1 = a1
        if a2 != None:
            self.a2 = a2
        else:
            self.a2 = None
            
        if a3 != None:
            self.a3 = a3
        else:
            self.a3 = None
        
        #number of dimensions
        self.NDIM = 3

        #start with empty sublattices
        self.sublattices =[]
        
        #atomic positions will be saved in numpy array
        self.positions = np.zeros((0))
        
        
    def add_one_sublattice(self, name, position, energy=0):
        """
        This methods add one sublattice to the cell
        
        Parameters
        ------------
        name : str 
            atom type/name 
        position: list 
            fractional coordinates of sublattice(s) respect to the primitive (or supercell) vectors
        energy: float
            onsite energy of the atom 
        """
        self.sublattices.append((name, position, energy))
        
    def add_sublattices(self, *sublattices):
        for sub in sublattices:
            self.add_one_sublattice(*sub)
            
            
    def getFromDataFile(self, datafile):
        """
        This method takes data from a data file (MD simulations, or experiments).
        
        -------------
        Data input format: 
        atom_id  x y z onsite_energy
        x y z are in Angstrom unit
        
        
        Example:
        Pb 0.0000 0.0000 0.0000 (if known, otherwise set 0)
        S  2.0000 2.0000 1.0000 (if known, otherwise set 0)
        ..
        ..
        ..
        
        1. The method wil transform the data to fractional coordinate relative to unit
        cell lattice vectors
        2. Add each sublattice to the cell
        
        """
        
        #data = np.loadtxt(datafile)
        f = open(datafile)
        data = f.readlines()
        
        self.positions = np.zeros((len(data), self.NDIM))
        
        for line in range(len(data)):
            oneline = data[line].split()
            
            name = oneline[0]
            position = [float(oneline[1]), float(oneline[2]), float(oneline[3])]
            
            #for MD/experimental data, onsite energy is usually unknown
            #for now we set onsite energy =0
            energy = float(oneline[4]) 
            
   
            
            self.sublattices.append((name, position, energy))
            
            #We will do extensive array manipulations we need numpy for better performance, 
            #save positions as numpy 
            
            self.positions[line] = position
            
            
