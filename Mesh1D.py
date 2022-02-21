import numpy as np
import sys
import Quadrature as Q
from math import sqrt

class Mesh1d:
    """Class defining a one dimensional mesh :
    - Points
    - Cells
    - Faces
    - Boundary Faces
    """

    def __init__(self):
        self.Points=[]
        self.Cells=[]
        self.NPoints=0
        self.NCells=0
        self.a=0
        self.b=0

    def __init__(self,a,b,N):
        if ( a > b ):
            print ("Error in Mesh1d: a should be lower than b")
            print ("a = ",a,"\t b = ",b)
            sys.exit()
        if ( N <=0 ):
            print ("Error in Mesh1d: N must be a positive integer")
            print ("N = ",N)
            sys.exit()
        self.a=a
        self.b=b
        self.NCells = N
        self.NPoints = N+1
        self.Points=np.linspace(a,b,N+1)
        self.Cells=[]
        for i in range(N):
            cell=[]
            cell.append(i)
            cell.append(i+1)
            (self.Cells).append(cell)

    def Init(self,funcInit):
        uh=np.zeros(self.NPoints)
        for iPoint in range(0,self.NPoints):
            x=self.Points[iPoint]
            uh[iPoint]=funcInit(x)
        return uh

    def Exact(self,funcExact,t):
        uh=np.zeros(self.NPoints)
        for iPoint in range(0,self.NPoints):
            x=self.Points[iPoint]
            uh[iPoint]=funcExact(t,x)
        return uh
