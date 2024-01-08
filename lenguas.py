#Imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy import signal
plt.rcParams['axes.formatter.useoffset'] = False
import numpy as np
from scipy.optimize import curve_fit
import scipy as sp
from scipy import stats
import glob, os
import random

#Functions
from scipy import fftpack
def plot_espectro(senal,xmax,ymax):
    X = fftpack.fft(senal)
    plt.plot(fftpack.fftfreq(len(t[:]),1./365), np.abs(X))
    plt.ylim(0, ymax)
    plt.xlim(0, xmax)
    plt.xlabel('freq [años$^{-1}$]')
    plt.show()

from scipy.signal import filtfilt
def PasaBajos(senal):
  fs = 365 #pienso en años, como la unidad es el año, mi frec de muestreo es 365 muestras por año (señal en días)
  frec_c = 1.2 #Tres muestras por año o sea frecuencia de corte es un trimestre (2 sería 6 meses y es demasiado)
  nyq = .5*fs #Nysquits para pasarle al filtro
  frec = frec_c / nyq
  orden = 6
  b,a = scipy.signal.butter(orden,frec,'low',analog=False,output='ba') #Filtro butterworth
  y = scipy.signal.filtfilt(b,a,senal,axis=0)
  return y

def PasaAltos(senal):
  fs = 365 #pienso en años, como la unidad es el año, mi frec de muestreo es 365 muestras por año (señal en días)
  frec_c = 0.001 #Tres muestras por año o sea frecuencia de corte es un trimestre (2 sería 6 meses y es demasiado)
  nyq = .5*fs #Nysquits para pasarle al filtro
  frec = frec_c / nyq
  orden = 6
  b,a = scipy.signal.butter(orden,frec,'high',analog=False,output='ba') #Filtro butterworth
  y = scipy.signal.filtfilt(b,a,senal,axis=0)
  return y

def Polo(senal):
  fs = 365 #pienso en años, como la unidad es el año, mi frec de muestreo es 365 muestras por año (señal en días)
  nyq = .5*fs #Nysquits para pasarle al filtro
  z, p, k = signal.butter(0, 0.1/nyq, output='zpk', fs=fs)
  y = signal.filtfilt(z,p,senal)
  return y


def close_returns(senal):
    cr =  np.zeros((len(senal),5000)) #cr de close returns
    eps = (np.max(senal)-np.min(senal))*0.005
    for j in range(5000):
      for i in range(len(senal)-5000):
        if (np.abs(senal[i+j] - senal[i]) < eps):
          cr[i,j] = 1
    return cr

def embedding(tau,dato):
     w = 3
     tau = tau
     embedding = dato[(np.arange(w)*(tau+1))+ np.arange(np.max(dato.shape[0] - (w-1)*(tau+1), 0)).reshape(-1,1)]
     dim1 = embedding[:,0]
     dim2 = embedding[:,1]
     dim3 = embedding[:,2]
     return dim1, dim2, dim3

def remover_ma(dato,ventana):
  numbers_series = pd.Series(dato)
  window_size = ventana
  windows = numbers_series.rolling(window_size)
  moving_averages = windows.mean()
  moving_averages_list = moving_averages.tolist()
  without_nans = moving_averages_list[window_size - 1:]
  ones = np.ones(int((len(dato)-len(without_nans)+1)/2))*np.mean(dato)
  moving_average = np.concatenate([ones,without_nans,ones])
  correccion = np.ones(len(moving_average))*np.mean(moving_average) - moving_average
  corregido = dato + correccion[:-1]
  return corregido, moving_average


def plot_orbita(df,inicio,fin,color='-b'):
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_subplot(111, projection='3d')
  ax.plot(df[0][inicio:fin], df[1][inicio:fin], df[2][inicio:fin],color,linewidth=3)
  #ax.set_xlim(23,30)
  #ax.set_ylim(23,30)
  #ax.set_zlim3d(23,30)
  return ax

def plot_orbita(df,inicio,fin,color='-b'):
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_subplot(111, projection='3d')
  ax.plot(df[0][inicio:fin], df[1][inicio:fin], df[2][inicio:fin],color,linewidth=3)
  #ax.set_xlim(23,30)
  #ax.set_ylim(23,30)
  #ax.set_zlim3d(23,30)
  return ax

def find_repeated_fragments(time_series, tolerancia=0.005, lag=25):
    fragmentos_repetidos = []
    series_length = len(time_series)

    for i in range(series_length - lag):
        fragmento = time_series[i:i + lag + 1]
        for j in range(i + lag + 1, series_length - lag):
            candidato = time_series[j:j + lag + 1]
            if len(fragmento) == len(candidato) and all(abs(fragmento[k] - candidato[k]) <= tolerancia for k in range(len(fragmento))):
                fragmentos_repetidos.append((i, j))

    return fragmentos_repetidos

def find_repeated_fragments_barato(time_series, T, tolerancia=0.005, lag=25):
    fragmentos_repetidos = []
    series_length = len(time_series)
    length_orbit = []
    cont = 0
    while cont <= 10:
        for i in range(50 - lag):
            fragmento = time_series[i:i + lag + 1]
            for j in range(i + lag + 1, series_length - lag):
                candidato = time_series[j:j + lag + 1]
                if len(fragmento) == len(candidato) and all(abs(fragmento[k] - candidato[k]) <= tolerancia for k in range(len(fragmento))) and (round((j - i)/36.5,2) < 10):
                    len_orbit = round((j - i)/T,0)
                    length_orbit.append(len_orbit)
                    cont += 1
                    print(cont)

    return length_orbit

    len_orbit = round((orbita[1] - orbita[0])/36.5,2)
    if ((np.abs(orbita[0] - orbita0[0]) > 30) and (np.abs(orbita[1] - orbita0[1]) > 30) and (len_orbit<10)):
        orbitas_posta.append(orbita)
        orbita0 = orbita


def rk4(dv,v,n,t,dt,e1,A,e2,omega):
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
    dv(v1, k1,e1,A,e2,omega)
    for x in range(0, n):
        v1[x]=v[x]+dt2*k1[x]
    dv(v1, k2,e1,A,e2,omega)
    for x in range(0, n):
        v1[x]=v[x]+dt2*k2[x]
    dv(v1, k3,e1,A,e2,omega)
    for x in range(0, n):
        v1[x]=v[x]+dt*k3[x]
    dv(v1, k4,e1,A,e2,omega)
    for x in range(0, n):
        v1[x]=v[x]+dt*k4[x]
    for x in range(0, n):
        v[x]=v[x]+dt6*(2.0*(k2[x]+k3[x])+k1[x]+k4[x])
    return v


def ecuaciones(v, dv,eps1,A,eps2,omega):
    """Defino la ecuación diferencial que quiero ajustar"""
    x=v[0]
    y=v[1]
    z=v[2]
    #A, eps2,omega = 7.5,1,4.7 #Cerca de la lengua de arnold de periodo 4, pero donde hay caos
    #A=A*(1+4.0*random.normalvariate(0,0.5)) #Agrego ruido a los parámetros
    dv[0]=y
    dv[1]=x-y-x*x*x+x*y+eps1+(eps2)*x*x+A*np.cos(z)
    dv[2]=omega
    return dv

def integracion(eps1=0.25,A=7.5,eps2=1,omega=4.7):
    for m in range(1):
        n=3 #Cantidad de variables   
        v=[]
        for x in range(0, n):
            v.append(x)

        v[0]=1.*(1+0.5*random.normalvariate(0,0.5))
        v[1]=3.*(1+0.5*random.normalvariate(0,0.5))
        v[2]=2.*(1+0.5*random.normalvariate(0,0.5))
        dt=2*np.pi / 365 / omega
        t=0.0
        t_pre=0.0
        t_max=1000.0
        x=[]
        y=[]
        z=[]
        cont=0
        while t<t_max:
            #eps1 = 0.25
            rk4(ecuaciones,v,n,t,dt,eps1,A,eps2,omega)
            t+=dt
            x.append(cont)  #ACÁ ARMO LOS ARREGLOS DE X Y Z CON LOS RESULTADOS QUE VA LARGANDO "V"
            y.append(cont)
            z.append(cont)
            x[cont]=v[0]
            y[cont]=v[1]
            z[cont]=v[2]
            cont=cont+1

        x = x[10000:40000]
        return x


def count_orbit(A_array,omega_array,file_path):
    matrix_A_omega = np.zeros((len(A_array),len(omega_array)))
    with open(file_path,'w') as f:
        for i,A in enumerate(A_array):
            for j, omega in enumerate(omega_array):
                int_x = np.array(integracion(0.2,A,1,omega))
                orbit_lengths = []
                orbitas = find_repeated_fragments(int_x[:3000:5])
                if len(orbitas) == 0:
                    matrix_A_omega[i,j] = 0
                else:
                    orbitas = find_repeated_fragments(int_x[::5])
                    orbita0 = orbitas[0]
                    for orbita in orbitas:
                        len_orbit = round((orbita[1] - orbita[0])/(365/5),2)
                        if ((np.abs(orbita[0] - orbita0[0]) > 100) and (np.abs(orbita[1] - orbita0[1]) > 100) and ((len_orbit < 10) and (len_orbit>0.8))):
                            orbit_lengths.append(len_orbit)
                            orbita0 = orbita

                    orbit_lengths = np.array(orbit_lengths)
                    #list(np.around(orbit_lengths,0)).count(4)
                    orbit_numbers=set(np.around(orbit_lengths,0))
                    list_numbers = list(orbit_numbers)
                    if len(list_numbers) == 1:
                        matrix_A_omega[i,j] = list_numbers[0]
                        print('save',i,' de ',len(A_array)*len(omega_array))

            np.savetxt(f, matrix_A_omega, fmt='%.2f')

    return matrix_A_omega

A_array = np.arange(7,12,0.01)
omega_array = np.arange(3,6,0.01)
print(len(A_array)*len(omega_array))
path_here = os.getcwd()
print('sale!!!!')
path_save = path_here+'/lengua_'+str(A_array[0])+'_'+str(A_array[-1])+'_'+str(omega_array[0])+'_'+str(omega_array[-1])+'.txt'
matrix_salida = count_orbit(A_array,omega_array,path_save)
df = pd.DataFrame(matrix_salida)
df.to_csv(path_here+'/lengua_'+str(A_array[0])+'_'+str(A_array[-1])+'_'+str(omega_array[0])+'_'+str(omega_array[-1])+'.csv', index=False)


