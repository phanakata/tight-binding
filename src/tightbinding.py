from lattice import Lattice
import math 
import numpy as np
import time, sys

"""Tight binding class
This class is used to model electronic systems based on their atomic positions and their interactions. 

"""
class TightBinding(Lattice):
    def __init__(self, lattice):
        """
        Parameters
        --------------
        lattice : lattice object 
        """
        super().__init__(lattice)
        self.lattice = lattice
        
        #create an empty neighbor list
        self.nlist = []
        self.N = len(self.lattice.positions)
        self.H = np.zeros((self.N, self.N))
        
        
    
    #Main methods 
    def createHamiltonian(self):
        """A main method
        1. Find neighbor(s) 
        2. Construct Hamiltonian
        """
       
        self.findNeighbors()
        self.hamiltonian()
        print("..ALL jobs are DONE!!..")
    
    
    
    #Helper functions
    def hamiltonian(self):
        """Function to construct Hamiltonian"""
        for i in range(self.N):
            
            time.sleep(0.1)
            self.update_progress("Constructing Hamiltonian", i/float(self.N))
            
            #onsite energy
            H[i][i] = self.lattice.sublattices[i][2]

            neighborList = len(self.nlist[i])
            for neighbor in range(neighborList):
                #Hij
                self.H[i][self.nlist[i][neighbor]]=self.tij(self.lattice.positions[i], self.lattice.positions[self.nlist[i][neighbor]])
        
        self.update_progress("Constructing Hamiltonian", 1) 
    
    def findNeighbors(self):
        """Function to find (nearest) neighbor(s)"""
        
        for i in range(self.N):
            self.nlist.append([])
    
        for i in range(self.N):
            time.sleep(0.1)
            self.update_progress("Finding neighbor(s)", i/float(self.N))
            for j in range(self.N):
           
            
                xi = self.lattice.positions[i]
                xj = self.lattice.positions[j]
                if self.lattice.pbc is True:
                    if self.distance2(xj, xi)<cut*cut and i !=j:
                        self.nlist[i].append(j)
                    else:
                        if abs(xi[0]-xj[0])>Lx/2:
                            xj[0] = xj[0] - Lx * (xj[0]-xi[0])/abs(xj[0]-xi[0])
                
                        if abs(xi[1]-xj[1])>Ly/2:
                            xj[1] = xj[1] -  Ly * (xj[1]-xi[1])/abs(xj[1]-xi[1])
                
                        if self.distance2(xj, xi)<cut*cut and i !=j:
                            self.nlist[i].append(j)
                
                else:
                    if self.distance2(xj, xi)<cut*cut and i !=j:
                        self.nlist[i].append(j)
            
            
        
        self.update_progress("Finding neighbor(s)", 1)
                
            
        
    
    
          
    @staticmethod    
    def tij(xi, xj):
        """Helper function to calculate tij
        Form simplicity tij is assumed to be proportional to 1/r^3 Harrison's rule
        Different functionals may be considered 
        """
        
        rij2 = self.distance2(xj, xi)
        if  rij2 <cut*cut:
            return t * (d_cc/math.sqrt(rij2))**3
        else:
            if self.lattice.pbc is True:
                if abs(xj[0]-xi[0])>Lx/2:
                    xj[0] = xj[0] - Lx * (xj[0]-xi[0])/abs(xj[0]-xi[0])
                if abs(xj[1]-xi[1])>Ly/2:
                    xj[1] = xj[1] -  Ly * (xj[1]-xi[1])/abs(xj[1]-xi[1])
                
                rij2 = self.distance2(xj, xi)
                if  rij2 <cut*cut:
                    return t * (d_cc/math.sqrt(rij2))**3



    @staticmethod            
    def distance2(xj, xi):
        """Helper function to calculate squared distance"""
        dx = xj[0]-xi[0]
        dy = xj[1]-xi[1]
        dz = xj[2]-xi[2]
        
        return dx*dx + dy*dy + dz*dz
    
    @staticmethod            
    def update_progress(job_title, progress):
        """Helper function to calculate calculation progress"""
        length = 20 # modify this to change the length
        block = int(round(length*progress))
        msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block), round(progress*100, 2))
        if progress >= 1: msg += " DONE\r\n"
        sys.stdout.write(msg)
        sys.stdout.flush()
    
 
                
            
        
