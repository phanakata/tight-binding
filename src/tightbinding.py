from lattice import Lattice
import math 
import numpy as np
import sys #,time
from time import time


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
        t0 = time()
        self.findNearestNeighbors2()
        print ("find NearestNeighbors2 time:", round(time()-t0, 3), "s")
        
        self.hamiltonian()
        print("..ALL jobs are DONE!!..")
    
    
    def benchmarkAlgo(self):
        """
        For development
        Method to time and benchmark differenrt algo 
        """

        t1 = time()
        self.findNearestNeighbors_z()
        print ("find NearestNeighbors time:", round(time()-t1, 3), "s")

        t1 = time()
        self.findNearestNeighbors()
        print ("find NearestNeighbors time:", round(time()-t1, 3), "s")
        
        #t1 = time()
        #self.findNearestNeighbors_np()
        #print ("find NearestNeighbors_np time:", round(time()-t1, 3), "s")
        
        
        #t1 = time()
        #self.findNearestNeighbors()
        #print ("find NearestNeighbors time:", round(time()-t1, 3), "s")
        
        #print(self.nlist)
        
        #self.nlist= []
        
        t0 = time()
        self.findNearestNeighbors2()
        print ("find NearestNeighbors2 time:", round(time()-t0, 3), "s")
        
  
      



    
    #Helper functions
    def hamiltonian(self):
        """Function to construct Hamiltonian"""
        for i in range(self.N):
            
            #time.sleep(0.1)
            self.update_progress("Constructing Hamiltonian", i/float(self.N))
            
            #onsite energy
            self.H[i][i] = self.lattice.sublattices[i][2]

            neighborList = len(self.nlist[i])
            for neighbor in range(neighborList):
                #Hij
                self.H[i][self.nlist[i][neighbor]]=self.tij_nn(self.lattice.positions[i], self.lattice.positions[self.nlist[i][neighbor]])
        
        self.update_progress("Constructing Hamiltonian", 1) 
    
    def findNearestNeighbors(self):
        """Function to find (nearest) neighbor(s)"""
        
        #to speed caclulation time
        #go to next loop whne nnearest = 3
        cut = []
        
        #append nearest neighbor cut(s) parameters
        for nn in range(len(self.lattice.parameters)):
            cut.append(self.lattice.parameters[nn][2])

        for i in range(self.N):
            self.nlist.append([])
    
        for i in range(self.N):
            #time.sleep(0.1)
            self.update_progress("Finding neighbor(s)", i/float(self.N))
            
            xi = self.lattice.positions[i]
            j = 0
            while j<self.N and len(self.nlist[i])<len(cut):
                
                xj = self.lattice.positions[j]
                nn = len(self.nlist[i])   #nearest neighbor index
                
                if i !=j and j not in self.nlist[i]:
                    if self.lattice.pbc is False:
                        if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                            self.nlist[i].append(j)
                            #this also means i is nearest for jth
                            self.nlist[j].append(i)
                    elif self.lattice.pbc is True:
                        if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                            self.nlist[i].append(j)
                            self.nlist[j].append(i)
                        else:
                            #pbc image
                            if abs(xi[0]-xj[0])>Lx/2:
                                xj[0] = xj[0] - Lx * (xj[0]-xi[0])/abs(xj[0]-xi[0])
                            if abs(xi[1]-xj[1])>Ly/2:
                                xj[1] = xj[1] -  Ly * (xj[1]-xi[1])/abs(xj[1]-xi[1])
                            
                            if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                                self.nlist[i].append(j)
                                self.nlist[j].append(i)
                                
                j = j +1  
        
        self.update_progress("Finding neighbor(s)", 1)
        
        
    def findNearestNeighbors_np(self):
        """Function to find (nearest) neighbor(s)"""
        
        #to speed caclulation time
        #go to next loop whne nnearest = 3
        cut = []
        
        
        
        #append nearest neighbor cut(s) parameters
        for nn in range(len(self.lattice.parameters)):
            cut.append(self.lattice.parameters[nn][2])

        self.nlist_np = np.zeros((self.N, len(self.lattice.parameters)))
    
        for i in range(self.N):
            #time.sleep(0.1)
            self.update_progress("Finding neighbor(s)", i/float(self.N))
            
            xi = self.lattice.positions[i]
            nn = 0
            j = 0
            while j<self.N and nn<len(cut):
                
                xj = self.lattice.positions[j]
                #nearest neighbor index
                
                if i !=j and j not in self.nlist_np[i]:
                    if self.lattice.pbc is False:
                        if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                            self.nlist_np[i][nn] = j 
                            #this also means i is nearest for jth
                            self.nlist_np[j][nn] = i
                            nn = nn + 1
                    elif self.lattice.pbc is True:
                        if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                            self.nlist_np[i][nn] = j 
                            #this also means i is nearest for jth
                            self.nlist_np[j][nn] = i
                            nn = nn + 1
                        else:
                            #pbc image
                            if abs(xi[0]-xj[0])>Lx/2:
                                xj[0] = xj[0] - Lx * (xj[0]-xi[0])/abs(xj[0]-xi[0])
                            if abs(xi[1]-xj[1])>Ly/2:
                                xj[1] = xj[1] -  Ly * (xj[1]-xi[1])/abs(xj[1]-xi[1])
                            
                            if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                                self.nlist_np[i][nn] = j 
                                #this also means i is nearest for jth
                                self.nlist_np[j][nn] = i
                                nn = nn + 1
                                
                                
                j = j +1  
        
        self.update_progress("Finding neighbor(s)", 1)    
        

    def findNearestNeighbors_z(self):
        """Function to find (nearest) neighbor(s)"""
        
        #to speed caclulation time
        #go to next loop whne nnearest = 3
        cut = []
        
        
        
        #append nearest neighbor cut(s) parameters
        for nn in range(len(self.lattice.parameters)):
            cut.append(self.lattice.parameters[nn][2])

        self.nlist_np = np.zeros((self.N, len(self.lattice.parameters)))
    
        for i in range(self.N):
            #time.sleep(0.1)
            self.update_progress("Finding neighbor(s)", i/float(self.N))
            
            xi = self.lattice.positions[i]
            nn = 0
            j = 0
            while j<self.N and nn<len(cut):
                
                xj = self.lattice.positions[j]
                #nearest neighbor index

                #dx = xj[0] - xi[0]
                #dy = xj[1] - xi[1]
                #dz = xj[2] - xi[2] 
                if self.lattice.pbc is False:
                    if i==j or abs(xj[0] - xi[0]) > cut or abs(xj[1] - xi[1])>cut or abs(xj[2] - xi[2])>cut:
                        continue 
                    elif self.distance2(xj, xi)<cut[nn]*cut[nn] and nn<len(cut):
                        self.nlist_np[i][nn] = j 
                        nn = nn + 1

                elif self.lattice.pbc is True:
                    if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                        self.nlist_np[i][nn] = j
                        #this also means i is nearest for jth
                        self.nlist_np[j][nn] = i
                        nn = nn + 1
                    else:
                        #pbc image
                        if abs(xi[0]-xj[0])>Lx/2:
                            xj[0] = xj[0] - Lx * (xj[0]-xi[0])/abs(xj[0]-xi[0])
                        if abs(xi[1]-xj[1])>Ly/2:
                            xj[1] = xj[1] -  Ly * (xj[1]-xi[1])/abs(xj[1]-xi[1])
                            
                        if self.distance2(xj, xi)<cut[nn]*cut[nn]:
                            self.nlist_np[i][nn] = j 
                            #this also means i is nearest for jth
                            self.nlist_np[j][nn] = i
                            nn = nn + 1
                                
                                
                j = j +1  
        
        self.update_progress("Finding neighbor(s)", 1)    
        
        
    def findNearestNeighbors2(self):
        """Function to find (nearest) neighbor(s)"""
        
        #to speed caclulation time
        #go to next loop whne nnearest = 3
        #nearestMax = 3 

        cut = self.lattice.parameters[0][2]
        
        for i in range(self.N):
            self.nlist.append([])
    
        for i in range(self.N):
            self.update_progress("Finding neighbor(s)", i/float(self.N))
            
            xi = self.lattice.positions[i]
            for j in range(self.N):
                
                xj = self.lattice.positions[j]
                
                if self.lattice.pbc is False:
                    if self.distance2(xj, xi)<cut*cut and i !=j:
                        self.nlist[i].append(j)
                
                elif self.lattice.pbc is True:
                    if self.distance2(xj, xi)<cut*cut and i !=j:
                        self.nlist[i].append(j)
                    else:
                        if abs(xi[0]-xj[0])>Lx/2:
                            xj[0] = xj[0] - Lx * (xj[0]-xi[0])/abs(xj[0]-xi[0])
                
                        if abs(xi[1]-xj[1])>Ly/2:
                            xj[1] = xj[1] -  Ly * (xj[1]-xi[1])/abs(xj[1]-xi[1])
                
                        if self.distance2(xj, xi)<cut*cut and i !=j:
                            self.nlist[i].append(j)
                
            
        
        self.update_progress("Finding neighbor(s)", 1)
            
        
    
    
          
    @staticmethod    
    def tij_nn(xi, xj):
        """Helper function to calculate NEAREST NEIGHBOR tij 
        Form simplicity tij is assumed to be proportional to 1/r^3 Harrison's rule
        Different functionals may be considered 
        """
        
        #Assign nearest neighbor values, index 0 mean first nearest neighbor
        d_cc = self.lattice.parameters[0][1]
        cut = self.lattice.parameters[0][2]
        t = self.lattice.parameters[0][3]

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
    
 
