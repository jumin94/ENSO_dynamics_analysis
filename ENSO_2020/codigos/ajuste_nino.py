#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:22:22 2019

@author: gabrielmindlin
"""

import numpy as np
from sklearn import linear_model

xyz=np.loadtxt('nino.dat')

x=xyz[:,0]
y=xyz[:,1]
z=xyz[:,2]
xr=np.roll(x,1)
yr=np.roll(y,1)
zr=np.roll(z,1)

dxdt=1000*(xr-x)
dydt=1000*(yr-y)
dzdt=1000*(zr-z)


theta=np.zeros((xyz.shape[0],12))
theta[:,0]=np.ones_like(x)
theta[:,1]=x
theta[:,2]=y
theta[:,3]=z
theta[:,4]=x*x
theta[:,5]=x*y
theta[:,6]=y*y
theta[:,7]=x*x*x
theta[:,8]=y*y*x
theta[:,9]=x*x*y
theta[:,10]=y*y*y
theta[:,11]=np.cos(z)

print('ajuste de x \n')
clf=linear_model.Lasso(alpha=0.1,max_iter=100000,fit_intercept=False,normalize=False)
clf.fit(-theta,dxdt)
print(clf.coef_)
print(clf.intercept_)
print('\n')

print('ajuste de y \n')
clf=linear_model.Lasso(alpha=0.1,max_iter=100000,fit_intercept=False,normalize=False)
clf.fit(-theta,dydt)
print(clf.coef_)
print(clf.intercept_)
print('\n')

print('ajuste de z \n')
clf=linear_model.Lasso(alpha=0.1,max_iter=100000,fit_intercept=False,normalize=False)
clf.fit(-theta,dzdt)
print(clf.coef_)
print(clf.intercept_)


