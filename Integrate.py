from math import sqrt
#import Quadrature as Q
#import Mesh1D as M

def integrate(func,quad,mesh):
        intVal = 0.
        for iCell in range(mesh.NCells):
                cell = mesh.Cells[ iCell ]
                xL=mesh.Points[ cell[0] ]
                xR=mesh.Points[ cell[1] ]
                for iqf in range(quad.nQuad):
                        xqf = xL + (xR-xL) * quad.points[iqf]
                        intVal += (xR-xL) * quad.weights[iqf] * func( xqf )
        return intVal

def ErrorInf(uh,funcExact,quad,mesh):
        error = 0.
        for iPoint in range(mesh.NPoints):
            x = mesh.Points[iPoint]
            error = max(error,abs(uh[iPoint]-funcExact(x)))
        return error

def ErrorL1(uh,funcExact,quad,mesh):
    error = 0.
    #print(quad)
    for iCell in range(mesh.NCells):
        cell = mesh.Cells[ iCell ]
        xL=mesh.Points[ cell[0] ]
        xR=mesh.Points[ cell[1] ]
        uL=uh[cell[0]]
        uR=uh[cell[1]]
        for iqf in range(quad.nQuad):
            xqf = xL + (xR-xL) * quad.points[iqf]
            uxqf = (uR - uL)/(xR-xL) * (xqf - xL) + uL
            error += (xR-xL) * quad.weights[iqf] * abs( funcExact( xqf ) - uxqf )
    return error

def ErrorL2(uh,funcExact,quad,mesh):
    error = 0.
    for iCell in range(mesh.NCells):
        cell = mesh.Cells[ iCell ]
        xL=mesh.Points[ cell[0] ]
        xR=mesh.Points[ cell[1] ]
        uL=uh[cell[0]]
        uR=uh[cell[1]]
        for iqf in range(quad.nQuad):
            xqf = xL + (xR-xL) * quad.points[iqf]
            uxqf = (uR - uL)/(xR-xL) * (xqf - xL) + uL
            error += (xR-xL) * quad.weights[iqf] * ( funcExact( xqf ) - uxqf )**2
    return sqrt(error)

def ErrorInfT(uh,funcExact,t,quad,mesh):
    error = 0.
    for iPoint in range(mesh.NPoints):
        x = mesh.Points[iPoint]
        error = max(error,abs(uh[iPoint]-funcExact(t,x)))
    return error

    
def ErrorL1T(uh,funcExact,t,quad,mesh):
    error = 0.
    for iCell in range(mesh.NCells):
        cell = mesh.Cells[ iCell ]
        xL=mesh.Points[ cell[0] ]
        xR=mesh.Points[ cell[1] ]
        uL=uh[cell[0]]
        uR=uh[cell[1]]
        for iqf in range(quad.nQuad):
            xqf = xL + (xR-xL) * quad.points[iqf]
            uxqf = (uR - uL)/(xR-xL) * (xqf - xL) + uL
            error += (xR-xL) * quad.weights[iqf] * abs( funcExact( t ,xqf) - uxqf )
    return error

                          
def ErrorL2T(uh,funcExact,t,quad,mesh):
    error = 0.
    for iCell in range(mesh.NCells):
        cell = mesh.Cells[ iCell ]
        xL=mesh.Points[ cell[0] ]
        xR=mesh.Points[ cell[1] ]
        uL=uh[cell[0]]
        uR=uh[cell[1]]
        for iqf in range(quad.nQuad):
            xqf = xL + (xR-xL) * quad.points[iqf]
            uxqf = (uR - uL)/(xR-xL) * (xqf - xL) + uL
            error += (xR-xL) * quad.weights[iqf] * ( funcExact( t , xqf) - uxqf )**2
    return sqrt(error)