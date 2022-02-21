import numpy as np
from matplotlib import pyplot as plt
import Mesh1D as M
import Quadrature as Q
import Integrate as I

###FUNCTION###

def u0(x):
    return np.exp(-((x-x0)/d)**2)

def uex(x,t):
    y=0
    if x-a*t >= 0:
        y=u0(x-a*t)
    else :
        y=0
    return y

def initA(mesh,CFL):
    M=np.zeros((mesh.NPoints,mesh.NPoints))
    M[0,0]=1
    M[0,1]=0
    M[mesh.NPoints-1,mesh.NPoints-2]=0
    M[mesh.NPoints-1,mesh.NPoints-1]=1
    for i in range(1,mesh.NPoints-1):
        M[i,i]=1
    for j in range(1,mesh.NPoints-1):
        if i==j-1:
            M[i,j]=-CFL*0.5
        if i==j+1:
            M[i,j]=CFL*0.5
    return M

def initb(mesh,un):
    b=np.zeros(mesh.NPoints)
    b[0]=0
    b[mesh.NPoints-1]=0
    for i in range(1,mesh.NPoints-1):
        b[i]=un[i]
    return b

def SolveEulerImp(mesh,uinit,T,dt,CFL):
    un=np.zeros(mesh.NPoints)
    un1=np.zeros(mesh.NPoints)
    tt=np.linspace(mesh.a,mesh.b,int(T/dt))
    A=initA(mesh,CFL)
    b=initb(mesh,uinit)
    for n in tt:
        un1=np.linalg.solve(A,b)
        b=initb(mesh,un1)
        un=un1
    return un


###MAIN SCRIPT###
N=200
a=1
x0=0.2
d=0.05
T=0.5
mesh=M.Mesh1d(0,1,N-1)
quad=Q.Gauss2()
dx=(mesh.b-mesh.a)/mesh.NPoints
dt=0.5*dx
CFL1=a*(dt/dx)
CFL2=a*(dt/dx)*4
uinit=np.zeros(mesh.NPoints)
U1=np.zeros(mesh.NPoints)

for i in range(mesh.NPoints):
    uinit[i]=u0(mesh.Points[i])
    U1[i]=uex(mesh.Points[i],T)

ySol1=SolveEulerImp(mesh,uinit,dt,dt,CFL1)
ySol2=SolveEulerImp(mesh,uinit,dt,dt,CFL2)
###PLOT###
print(CFL1)
print(CFL2)
plt.plot(mesh.Points,uinit,c="red",label="Initial condition")
plt.plot(mesh.Points,U1, ls=":",label="Exact solution")
plt.plot(mesh.Points,ySol1,label="Approx solution CFL=0.05")
plt.grid()
plt.legend()
plt.show()