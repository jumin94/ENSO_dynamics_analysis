#@title
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import scipy
from scipy import signal
from scipy.optimize import curve_fit
import scipy as sp
from scipy import stats
plt.rcParams['axes.formatter.useoffset'] = False
import glob
import os

#FILTRO-----------------------------------------
from scipy.signal import filtfilt
def PasaBajos(senal):
    """Aplica un filtro butterworth
    entrada: senial (np.array)
    salida: senial filtrada (np.array)
    """
    fs = 365 #pienso en años, como la unidad es el año, mi frec de muestreo es 365 muestras por año (señal en días)
    frec_c = 1.2 #Tres muestras por año o sea frecuencia de corte es un trimestre (2 sería 6 meses y es demasiado)
    nyq = .5*fs #Nysquits para pasarle al filtro
    frec = frec_c / nyq 
    orden = 6
    b,a = scipy.signal.butter(orden,frec,'low',analog=False,output='ba') #Filtro butterworth
    y = scipy.signal.filtfilt(b,a,senal,axis=0)
    return y

def remover_ma(dato,ventana):
    """Aplica una ventana movil y la remueve a la senial
    entrada: senial (np.array), largo de la ventana (escalar)
    salida: senial filtrada (np.array), ventana movil (componente baja frec)
    """
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

# Plot espectro
from scipy import fftpack
def plot_espectro(senal):
    """Plot del espectro de la senial
    entrada: senial (np.array)
    salida: figura
    """
    N = len(senal)
    rate = 365
    L = N/rate
    t = np.arange(N) / rate
    X = fftpack.fft(senal)
    fig = plt.figure()
    plt.plot(fftpack.fftfreq(len(t[:]),1./365), np.abs(X))
    plt.ylim(0, 30000)
    plt.xlim(0, 6)
    plt.xlabel('freq [años$^{-1}$]')
    return fig
    
#CLOSE RETURNS-----------------------------

def close_returns(senal):
    """Calcula el mapa de close returns
    entrada: senial (np.array 1D)
    salida: matriz (i,p) con los close returns (np.array 2D)
    """
    cr =  np.zeros((len(senal),2000)) #cr de close returns
    eps = (np.max(senal)-np.min(senal))*0.001
    for j in range(2000):
      for i in range(len(senal)-2000):
        if (np.abs(senal[i+j] - senal[i]) < eps) and (np.abs(senal[i+j+4] - senal[i+4]) < eps):
          cr[i,j] = 1
    return cr

def find_indices(matrix):
    """Genera los indices x,y de las posiciones de los unos en la matriz
    de close returns
    entrada: matriz (np.array 2D)
    salida: lista con la posicion (x,y) de cada uno (list)
    """
    indices = np.where(matrix == 1)
    return list(zip(indices[0], indices[1]))

# Plot close returns-----------------------------
def close_returns_for_plot(senal):
    """Overkill mal pero genera arrays a partir de la lista de (x,y)
    entrada: lista (x,y) de los close returns (list)
    salida: vectores x, y (dos np.array)
    """
    cr = close_returns(senal)
    result = find_indices(cr)
    x = []
    y = []
    for m in result:
        x.append(m[0])
        y.append(m[1])
        
    return x, y


def figure_cr(x,y,ylim=1600,xlim=15000,hist_lim=1000):
    
    # definitions for the axes
    left, width = 0.1, 2
    bottom, height = 0.1, 0.2
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # start with a rectangular Figure
    fig = plt.figure(figsize=(8, 8))

    ax_scatter = plt.axes(rect_scatter)
    ax_scatter.tick_params(direction='in', top=True, right=True)    
    ax_histy = plt.axes(rect_histy)
    ax_histy.tick_params(direction='in', labelleft=False)

    # the scatter plot:
    ax_scatter.scatter(x, y,s=0.8)

    # now determine nice limits by hand:
    binwidth = 0.25
    lim = np.ceil(np.abs([x, y]).max() / binwidth) * binwidth
    ax_scatter.set_xlim((0, xlim))
    ax_scatter.set_ylim((0, ylim))

    # Format the x-axis labels
    ax_scatter.set_yticks(np.arange(0,ylim,365))
    xtick_ylim = np.ceil(ylim/365)
    ax_scatter.set_yticklabels(np.arange(0,xtick_ylim,1).astype(int))
    ax_scatter.set_xticks(np.arange(0,xlim,365))
    xtick_xlim = np.ceil(xlim/365)
    ax_scatter.set_xticklabels(np.arange(0,xtick_xlim,1).astype(int))
    # Rotate the x-axis labels for better readability (optional)

    # Add labels and title as needed
    ax_scatter.set_xlabel('i [T = 365 days]',fontsize=18)
    ax_scatter.set_ylabel('j [T = 365 days]',fontsize=18)

    bins = np.arange(0, ylim + 30, 30)
    #ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal') #
    #ax_histx.set_xlim(ax_scatter.get_xlim())
    ax_histy.set_ylim(ax_scatter.get_ylim())
    ax_histy.set_xlim(1,hist_lim)
    ax_scatter.xaxis.set_ticks_position('top')
    ax_scatter.tick_params(axis='both', labelsize= 18)
    ax_histy.tick_params(axis='both', labelsize= 18)
    return fig       

#COMPLEJIDAD -------------------------
from itertools import permutations
import numpy as np

#Funcion para evaluar la permutacion a la que corresponde unna ventana
def evaluate_order(values):
    """Evalua el orden de los valores (permutacion)
    por ejemplo si los valores son x1 < x2 < x3, devuelve [0,1,2]
    entrada: valores a ordenar (np.array o lista)
    salida: lista de orden (list)
    """
    # Ordeno la lista
    sorted_values = sorted(values)
    # Mapeo los valores ordenados con su indice
    order_dict = {value: index + 1 for index, value in enumerate(sorted_values)}
    # Genero lista con permutacion 
    order = [order_dict[value] for value in values]
    return order

def entropia_y_complejidad(serie_tiempo,n):
    """Calcula la entropia y luego la complejidad de la serie
    entrada: senial, orden para evaluar entropia (n) (np.array,escalar)
    salida: entropia, complejidad (escalar, escalar)
    """
    #Defino largo de la ventana
    #n = 4
    #Permutaciones posibles 
    permutaciones_posibles = list(permutations(np.arange(1,n+1,1)))
    #Selecciono ventanas de largo n en la serie de tiempo
    ventanas = [serie_tiempo[i:i+n] for i in range(len(serie_tiempo) - n + 1)]
    #Evaluo la permutacion de cada ventana
    permutaciones = [evaluate_order(ventana) for ventana in ventanas]

    #Evaluo las probabilidades de cada permutacion posible y la guardo en un vector p_i
    probabilidades = np.zeros(len(permutaciones_posibles))
    for i,lista_a_contar in enumerate(permutaciones_posibles):
        # Inicializar el contador para la permutacion i
        contador = 0
        # Recorrer la lista de permutaciones de las ventanitas
        for sublista in permutaciones:
            if tuple(sublista) == lista_a_contar: 
                contador += 1      #Suma al contador cada vez que una ventanita corresponde a la permutacion i
        
        probabilidades[i] = contador/len(permutaciones) #Calcula la probabilidad de la permutacion i


    #Evaluo H normalizada, desequilibrio y complejidad
    p = np.array(probabilidades)
    p_sin_ceros = p[p!=0]
    H = -sum(p_sin_ceros*np.log2(p_sin_ceros))/np.log2(np.math.factorial(n))
    c = H*(sum(p-(1/len(p))**2))
    return H, c

def entropia_y_complejidad_wootters(serie_tiempo,n):
    """Calcula la entropia y luego la complejidad de la serie
    entrada: senial, orden para evaluar entropia (n) (np.array,escalar)
    salida: entropia, complejidad (escalar, escalar)
    
    La correccion de wooters consiste en evaluar la distancia a 1/N en el 
    espacio de probabilidades y por lo tanto con una formula mas compleja que
    la distancia Euclidea
    """
    #Defino largo de la ventana
    #n = 4
    #Permutaciones posibles 
    permutaciones_posibles = list(permutations(np.arange(1,n+1,1)))
    #Selecciono ventanas de largo n en la serie de tiempo
    ventanas = [serie_tiempo[i:i+n] for i in range(len(serie_tiempo) - n + 1)]
    #Evaluo la permutacion de cada ventana
    permutaciones = [evaluate_order(ventana) for ventana in ventanas]

    #Evaluo las probabilidades de cada permutacion posible y la guardo en un vector p_i
    probabilidades = np.zeros(len(permutaciones_posibles))
    for i,lista_a_contar in enumerate(permutaciones_posibles):
        # Inicializar el contador para la permutacion i
        contador = 0
        # Recorrer la lista de permutaciones de las ventanitas
        for sublista in permutaciones:
            if tuple(sublista) == lista_a_contar: 
                contador += 1      #Suma al contador cada vez que una ventanita corresponde a la permutacion i
        
        probabilidades[i] = contador/len(permutaciones) #Calcula la probabilidad de la permutacion i

    #Evaluo H normalizada, desequilibrio y complejidad
    p = np.array(probabilidades)
    p_sin_ceros = p[p!=0]
    H = -sum(p_sin_ceros*np.log2(p_sin_ceros))/np.log2(np.math.factorial(n))
    qw = (1/np.arccos((1/np.math.factorial(n))**(1/2)))*np.arccos(sum(p**(1/2)*(1/np.math.factorial(n))**(1/2)))
    c = H*qw
    return H, c



#EMBEDDING-----
def embedding(tau,dato):
     w = 3
     tau = tau
     embedding = dato[(np.arange(w)*(tau+1))+ np.arange(np.max(dato.shape[0] - (w-1)*(tau+1), 0)).reshape(-1,1)]
     dim1 = embedding[:,0]
     dim2 = embedding[:,1]
     dim3 = embedding[:,2]
     return dim1, dim2, dim3
 
 
def plot_orbita(df,inicio,fin,color='-g'):
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_subplot(111, projection='3d')
  ax.plot(df[0][inicio:fin], df[1][inicio:fin], df[2][inicio:fin],color,linewidth=3)
  #ax.set_xlim(23,30)
  #ax.set_ylim(23,30)
  #ax.set_zlim3d(23,30)
  return ax

#Integrar sistema dinamico-----------------------
import random
class Interar:
    def __init__(self):
        'Integro'
    
    def A_w_eps1(self,A,w,eps1,eps2,ruido1,ruido2,dt):
        self.A = A
        self.w = w
        self.eps1 = eps1
        self.eps2 = eps2
        self.ruido1 = 0.5
        self.ruido2 = 0.5
        self.dt = dt
        
    
    def ecuaciones(self,v, dv,eps1,eps2):
        """Defino la ecuación diferencial que quiero ajustar"""
        x=v[0]
        y=v[1]
        z=v[2]
        eps1, A, eps2,omega = eps1,self.A,eps2,self.w
        dv[0]=y
        dv[1]=x-y-x*x*x+x*y+eps1+(eps2)*x*x+A*np.cos(z)
        dv[2]=omega
        return dv
    
    def rk4(self,dv,v,n,t,dt,e1,e2):
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
        dv(v1, k1,e1,e2)
        for x in range(0, n):
            v1[x]=v[x]+dt2*k1[x]
        dv(v1, k2,e1,e2)     
        for x in range(0, n):
            v1[x]=v[x]+dt2*k2[x]
        dv(v1, k3,e1,e2)
        for x in range(0, n):
            v1[x]=v[x]+dt*k3[x]
        dv(v1, k4,e1,e2)
        for x in range(0, n):
            v1[x]=v[x]+dt*k4[x]        
        for x in range(0, n):
            v[x]=v[x]+dt6*(2.0*(k2[x]+k3[x])+k1[x]+k4[x])
        return v  
    
    def integracion(self):
        
        n=3 #Cantidad de variables   
        v=[]
        for x in range(0, n):
            v.append(x)

        v[0]=1.#*(1+self.ruido*random.normalvariate(0,0.5))
        v[1]=3.#*(1+self.ruido*random.normalvariate(0,0.5))
        v[2]=2.#*(1+self.ruido*random.normalvariate(0,0.5))
        dt=self.dt
        t=0.0
        t_pre=0.0
        t_max=1000.0
        x=[]
        y=[]
        z=[]
        cont=0
        while t<t_max:
            eps1 = self.eps1*(1+self.ruido1*random.normalvariate(0,0.5))
            eps2 = self.eps2*(1+self.ruido2*random.normalvariate(0,0.5))
            self.rk4(self.ecuaciones,v,n,t,dt,eps1,eps2)
            t+=dt
            x.append(cont)  #ACÁ ARMO LOS ARREGLOS DE X Y Z CON LOS RESULTADOS QUE VA LARGANDO "V"
            y.append(cont)
            z.append(cont)
            x[cont]=v[0]
            y[cont]=v[1]
            z[cont]=v[2]
            cont=cont+1

        x = x[:60000]
        return x
    
    
def vecinos(y, dim):
    #y contiene los puntos del embedding [[P1], [P2], [P3],...,[PN]]
    #dim = la dimension del embedding
    #Definimos el numero de puntos
    N = len(y)
    #Transformamos el vector y en un array
    y = np.array(y)
    #Inicializamos la matriz de distancia cuadrática entre puntos i j
    Dij = []
    #Calculamos las distancias
    for k in range(len(y)-1):
        #Armamos una matriz de "N-k-1" filas y "dim" columnas en las que se
        #repite el punto Pk al que le vamos a calcular la distancia con el resto
        My = np.full((N-k-1, dim), y[k])
        if dim ==1:
            #Si la dimension es 1, la distancia es restar (Pk - Pj)^2
            Dktodos = np.power((My.transpose() - y[k+1:])[0], 2, dtype=float)
        else:
            #Si la dimension es >1, la distancia es ((xk-xj)^2+(yk-yj)^2+...)
            #por eso hacemos la resta de componentes, elevamos al cuadrado
            #y despues sumamos
            Dktodos = np.sum(np.power(My-y[k+1:],2, dtype=float), 1, dtype=float)
        #La fila de la matriz va a ser ceros hasta el lugar k y despues empiezan
        #las distancias que calculamos
        fila = list(np.concatenate((np.zeros(k+1), Dktodos)))
        Dij.append(fila)
        #print("Porcentaje :", 100*k/N)
    #Al ultimo punto no hace falta calcularle las distancias. Agregamos una fila
    #de ceros
    Dij.append(list(np.zeros(N)))
    Dij = np.array(Dij)
    #Hasta aca tenemos una matriz triangular superior, con ceros en la diagonal
    #La matriz debe ser simetrica. Y para que al buscar la minima distancia
    #no nos agarre el punto consigmo mismo; le sumamos algo grande a la diagonal    
    Dij = Dij + Dij.transpose() + np.eye(N) *100
    #Ya teniendo la matriz de distancias cuadraticas, buscamos los indices
    #en los que esta el minimo de cada fila
    indice_min = np.argmin(Dij, 1)
    #Teniendo los indices, guardamos la minima distancia
    dist_min = []
    for k in range(N):
        dist_min.append(np.sqrt(Dij[k][indice_min[k]]))
    #Por utlimo, devolvemos los indices de minimo y el valor de la distancia
    return indice_min, dist_min



def porcentaje_falsos_vecinos(x):
    porcentaje_falsos_vecinos = []
    dim = 1

    T =15
    indice_min, dist_min = vecinos(x, dim)
    #Determino cuales son falsos vecinos para 1d
    R_crecimiento = []
    falsos_vecinos = 0
    puntos_noanalizados = 0
    for k in range(len(indice_min)-T):
        if indice_min[k]+dim*T < len(x):
            R_aux = np.abs( x[k+dim*T] - x[indice_min[k]+dim*T] ) / dist_min[k]
            R_crecimiento.append(R_aux)
            if R_aux >= 10:
                falsos_vecinos+=1
        else:
            puntos_noanalizados+=1

    porcentaje_falsos_vecinos.append(falsos_vecinos / (len(x)-T))
    
    y_emb_2 = []

    for k in range(len(x)-1*T):
        y_emb_2.append([x[k], x[k+T]])

    dim = 2
    indice_min, dist_min = vecinos(y_emb_2, dim)

    #Determino cuales son falsos vecinos para 2d
    R_crecimiento = []
    falsos_vecinos = 0
    puntos_noanalizados = 0
    for k in range(len(indice_min)-T):
        if indice_min[k]+dim*T < len(x):
            R_aux = np.abs( x[k+dim*T] - x[indice_min[k]+dim*T] ) / dist_min[k]
            R_crecimiento.append(R_aux)
            if R_aux >= 10:
                falsos_vecinos+=1
        else:
            puntos_noanalizados+=1
    print("El porcentaje de falsos vecinos es:", falsos_vecinos / (len(x)-T))
    print(puntos_noanalizados)
    porcentaje_falsos_vecinos.append(falsos_vecinos / (len(x)-T))

    y_emb_3 = []

    for k in range(len(x)-2*T):
        y_emb_3.append([x[k], x[k+T], x[k+2*T]])

    dim = 3
    indice_min, dist_min = vecinos(y_emb_3, dim)
    #Determino cuales son falsos vecinos para 3d
    R_crecimiento = []
    falsos_vecinos = 0
    puntos_noanalizados = 0
    for k in range(len(indice_min)-T):
        if indice_min[k]+dim*T < len(x):
            R_aux = np.abs( x[k+dim*T] - x[indice_min[k]+dim*T] ) / dist_min[k]
            R_crecimiento.append(R_aux)
            if R_aux >= 10:
                falsos_vecinos+=1
        else:
            puntos_noanalizados+=1
    print("El porcentaje de falsos vecinos es:", falsos_vecinos / (len(x)-T))
    print(puntos_noanalizados)
    porcentaje_falsos_vecinos.append(falsos_vecinos / (len(x)-T))


    y_emb_4 = []

    for k in range(len(x)-3*T):
        y_emb_4.append([x[k], x[k+T], x[k+2*T], x[k+3*T]])

    dim = 4
    indice_min, dist_min = vecinos(y_emb_4, dim)

    #Determino cuales son falsos vecinos para 4d
    R_crecimiento = []
    falsos_vecinos = 0
    puntos_noanalizados = 0
    for k in range(len(indice_min)-T):
        if indice_min[k]+dim*T < len(x):
            R_aux = np.abs( x[k+dim*T] - x[indice_min[k]+dim*T] ) / dist_min[k]
            R_crecimiento.append(R_aux)
            if R_aux >= 10:
                falsos_vecinos+=1
        else:
            puntos_noanalizados+=1
            
    porcentaje_falsos_vecinos.append(falsos_vecinos / (len(x)-T))

    y_emb_5 = []

    for k in range(len(x)-4*T):
        y_emb_5.append([x[k], x[k+T], x[k+2*T], x[k+3*T], x[k+4*T]])

    dim = 5
    indice_min, dist_min = vecinos(y_emb_5, dim)
    #Determino cuales son falsos vecinos para 4d
    R_crecimiento = []
    falsos_vecinos = 0
    puntos_noanalizados = 0
    for k in range(len(indice_min)-T):
        if indice_min[k]+dim*T < len(x):
            R_aux = np.abs( x[k+dim*T] - x[indice_min[k]+dim*T] ) / dist_min[k]
            R_crecimiento.append(R_aux)
            if R_aux >= 10:
                falsos_vecinos+=1
        else:
            puntos_noanalizados+=1

    porcentaje_falsos_vecinos.append(falsos_vecinos / (len(x)-T))
    return porcentaje_falsos_vecinos


    