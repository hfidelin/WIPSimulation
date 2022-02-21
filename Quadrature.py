import numpy as np
import sys
from math import sqrt

class Quadrature:
    """Class defining the quadrature formula :
    """

    def __init__(self):
        self.nQuad=0
        self.weights=[]
        self.points=[]
        self.exactitudeDegree=0
        self.order=0

    def __str__(self):
        toPrint=f"Quadrature formula {self.__class__.__name__},\t npoints = {self.npoints},\t weights = {self.weights}, \t points = {self.points}"
        return toPrint

class RectangleLeft(Quadrature):
    def __init__(self):
        self.nQuad=1
        self.exactitudeDegree=0
        self.order=1
        self.weights=[1.]
        self.points=[0.]
        

class RectangleRight(Quadrature):
    def __init__(self):
        self.nQuad=1
        self.exactitudeDegree=0
        self.order=1
        self.weights=[1.]
        self.points=[1.]

class MiddlePoint(Quadrature):
    def __init__(self):
        self.nQuad=1
        self.exactitudeDegree=1
        self.order=2
        self.weights=[1.]
        self.points=[0.5]

class Trapeze(Quadrature):
    def __init__(self):
        self.nQuad=2
        self.exactitudeDegree=1
        self.order=2
        self.weights=[0.5, 0.5]
        self.points=[0., 1.]

class Simpson(Quadrature):
    def __init__(self):
        self.nQuad=3
        self.exactitudeDegree=3
        self.order=4
        self.weights=[1./6., 4./6., 1./6.]
        self.points=[0., 0.5, 1.]

class Gauss2(Quadrature):
    def __init__(self):
        self.nQuad=2
        self.exactitudeDegree=3
        self.order=4
        self.weights=[0.5, 0.5]
        self.points=[0.5-0.5/sqrt(3), 0.5+0.5/sqrt(3)]