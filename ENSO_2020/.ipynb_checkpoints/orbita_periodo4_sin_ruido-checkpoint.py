def ecuaciones(v, dv,eps1):
    """Defino la ecuación diferencial que quiero ajustar"""
    x=v[0]
    y=v[1]
    z=v[2]
    eps1, A, eps2,omega = 0.2,10,1,4.7
    #A=A*(1+5.0*random.normalvariate(0,0.5))
    dv[0]=y
    dv[1]=x-y-x*x*x+x*y+eps1+(eps2)*x*x+A*np.cos(z)
    dv[2]=omega
    return dv

for m in range(5):
    n=3 #Cantidad de variables   
    v=[]
    for x in range(0, n):
        v.append(x)

    v[0]=1.*(1+0.5*random.normalvariate(0,0.5))
    v[1]=3.*(1+0.5*random.normalvariate(0,0.5))
    v[2]=2.*(1+0.5*random.normalvariate(0,0.5))
    dt=0.0036
    t=0.0
    t_pre=0.0
    t_max=1000.0
    x=[]
    y=[]
    z=[]
    print ("Ahí va el gráfico")
    #print(v[0], v[1], v[2])
    mpl.rcParams['legend.fontsize'] = 10
    fig=plt.figure()
    cont=0
    f=open("Nino_eps1_0.25+ruido_.dat","w")
    while t<t_max:
        eps1 = 0.25
        rk4(ecuaciones,v,n,t,dt,eps1)
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
    x = x[:60000]
    np.savetxt('nino_periodo4_sin_ruido_'+str(m)+'.txt', x)

