# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 11:43:02 2016

@author: Iván
"""

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
    dv[0]=10.0*(y-x)
    dv[1]=28.0*x-y-x*z
    dv[2]=x*y-8.0*z/3.0
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
v[0]=10.0
v[1]=10.0
v[2]=40.0
dt=0.001
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
f=open("Lorenz_new.dat","w")
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
ax=fig.add_subplot(111, projection= '3d')

ax.plot(x,y, z, label="física")

plt.show()

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    