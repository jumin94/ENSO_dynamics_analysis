import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt       

      


    
    
# --------------------
# Function definitions
# --------------------
def ecuaciones(v, dv):
    x=v[0]
    y=v[1]
    z=v[2]
    eps1,A, eps2,omega=0.2, 7.5,1,3.5
    dv[0]=y
    dv[1]=x-y-x*x*x+x*y+eps1+A*np.cos(z)+(eps2)*x*x
    dv[2]=omega
    return dv
def rk4(dv,v,n,t,dt):
    v1=[]
    k1=[]
    k2=[]
    k3=[]
    k4=[]
    for x in range(0, n):
        v1.append(x)
        k1.append(x)
        k2.append(x)
        k3.append(x)
        k4.append(x)

    dt2=dt/2.0
    dt6=dt/6.0
    for x in range(0, n):
        v1[x]=v[x]
    dv(v1, k1)
    for x in range(0, n):
        v1[x]=v[x]+dt2*k1[x]
    dv(v1, k2)     
    for x in range(0, n):
        v1[x]=v[x]+dt2*k2[x]
    dv(v1, k3)
    for x in range(0, n):
        v1[x]=v[x]+dt*k3[x]
    dv(v1, k4)
    for x in range(0, n):
        v1[x]=v[x]+dt*k4[x]        
    for x in range(0, n):
        v[x]=v[x]+dt6*(2.0*(k2[x]+k3[x])+k1[x]+k4[x])
    return v
n=3 #Cantidad de variables   
v=[]
for x in range(0, n):
    v.append(x)
v[0]=1.
v[1]=3.
v[2]=2.
dt=0.01
t=0.0
t_pre=0.0
t_max=100.0
x=[]
y=[]
z=[]
print ("Ahí va el gráfico")
#print(v[0], v[1], v[2])
mpl.rcParams['legend.fontsize'] = 10
fig=plt.figure()
cont=0
f=open("Nino.dat","w")
while t<t_max:
    rk4(ecuaciones,v,n,t,dt)
    t+=dt
    x.append(cont)  #ACÁ ARMO LOS ARREGLOS DE X Y Z CON LOS RESULTADOS QUE VA LARGANDO "V"
    y.append(cont)
    z.append(cont)
    x[cont]=v[0]
    y[cont]=v[1]
    z[cont]=v[2]
    cont=cont+1
    f.write(str(v[0])+"\t"+str(v[1])+"\t"+str(v[2])+"\n")
#    f.write(str(t)+"\n")
f.close()

x1d=np.zeros(len(x))
x2d=np.zeros(len(x))
x0d=np.zeros(len(x))
for i in range(30,len(x)):
    x1d[i]=-x[i-15]
    x2d[i]=-x[i-30]
    x0d[i]=-x[i]

ax=fig.add_subplot(111, projection= '3d')

ax.plot(x0d,x1d, x2d, label="física")

plt.show()


    
    
    
    
    
    