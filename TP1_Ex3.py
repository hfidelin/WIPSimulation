import numpy as np
from matplotlib import pyplot as plt
import Mesh1D as M
import Quadrature as Q
import Integrate as I

###FUNCTIONS###

def c01(x):
    return np.exp(-((x-15)**2)/s1)


def c02(x):
    return np.exp(-((x-15)**2)/s2)

def u(x):
    return 1+U0*np.sin(np.pi*x)

def F(x):
    R=1.111*(x-0.7)*(1-x)
    if x < 0.7 :
        R=0
    return R

def initA (mesh):
    A=np.zeros((mesh.NPoints,mesh.NPoints))
    A[0,0]=1
    A[-1,-1]=1
    for i in range(1,mesh.NPoints-1):
        A[i,i-1]=-U[i]*l+U[i]*m
        A[i,i  ]=1+2*U[i]*l
        A[i,i+1]=-U[i]*l-U[i]*m
    return A

def initb(mesh,cn):
    b=np.zeros(mesh.NPoints)
    for i in range(1,mesh.NPoints-1):
        b[i]=(U[i]*l-U[i]*m)*cn[i-1]+(1-2*U[i]*l)*cn[i]+(U[i]*l+U[i]*m)*cn[i+1]+Fvec[i]
    return b
    
def SolveADR(mesh,Cinit,T,dt):
    cn=np.zeros(mesh.NPoints)
    cn1=np.zeros(mesh.NPoints)
    A=initA(mesh)
    b=initb(mesh,Cinit)
    for n in range(int(T/dt)):
        cn1=np.linalg.solve(A,b)
        b=initb(mesh,cn1)
        cn[...]=cn1[...]
    return cn
###MAIN SCRIPT###

s1=8
s2=2
b=1
g=0.1
U0=0.92
L=30
T=25
N=10
mesh=M.Mesh1d(0,L,N-1)
quad=Q.Gauss2()
dx=(mesh.b-mesh.a)/mesh.NPoints
dt=0.1
l=(g*dt)/(2*b*dx*dx)
m=1/(4*dx)
Cinit1=mesh.Init(c01)
Cinit2=mesh.Init(c02)
U=np.zeros(mesh.NPoints)
Fvec=np.zeros(mesh.NPoints)
B=initb(mesh,Cinit1)
M=initA(mesh)

for i in range(mesh.NPoints):
    U[i]=u(mesh.Points[i])
    Fvec[i]=F(mesh.Points[i])

Sol1=SolveADR(mesh,Cinit1,dt,dt)

###PLOT###

plt.plot(mesh.Points,Cinit1, label="Initial condition with sigma =8")
#plt.plot(mesh.Points,Cinit2, label="Initial condition with sigma =2")
print("Initial solution vector size",Cinit1.size)

plt.plot(mesh.Points,Sol1,label="Approx solution")

print("Time vector size : ",U.size)
print("F vector size :",Fvec.size)
print("Doit afficher 4:",B.size)
print(U)
plt.legend()
plt.grid()
plt.show()