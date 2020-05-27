# -*- coding: utf-8 -*-
"""
Created on Sat Nov 16 14:16:35 2019

@author: jeff
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time
from numba import jit


#D = 1
#rho = 0.5


init_H = {(-1, 0): [1, 3], (0, 0.5): [4, 0, 2, 2], (1, 0): [3, 1], (-1, 1): [3.0, 1], (0, 1.5): [2.0], (1, 1): [1, 3], (2, 0.5): [2], (-2, 0.5): [2]}
init_U = {(-1, 0): np.array([[np.sqrt(0.5), np.sqrt(0.5)],
       [ -np.sqrt(0.5), np.sqrt(0.5)]]), (0, 0.5): np.array([[ 5.00000000e-01, -3.53553391e-01, -6.12372436e-01,
         5.00000000e-01],
       [-5.00000000e-01, -3.53553391e-01, -6.12372436e-01,
        -5.00000000e-01],
       [-np.sqrt(0.5),  1.87763707e-04, -1.08405427e-04,
         np.sqrt(0.5)],
       [ 0, -8.66025404e-01,  5e-01, 0]]), (1, 0): np.array([[ np.sqrt(0.5),  np.sqrt(0.5)],
       [-np.sqrt(0.5),  np.sqrt(0.5)]]), (-1, 1): np.array([[ np.sqrt(0.5),  np.sqrt(0.5)],
       [-np.sqrt(0.5),  np.sqrt(0.5)]]), (0, 1.5): np.array([[1.]]), (1, 1): np.array([[-np.sqrt(0.5), -np.sqrt(0.5)],
       [ np.sqrt(0.5), -np.sqrt(0.5)]]), (2, 0.5): np.array([[1.]]), (-2, 0.5): np.array([[1]])}


#Tk = (J*rho)**0.5*np.exp(-1/(2*D*J*rho))

E_init = init_H
U = init_U
N_cut = 100
pair1 = [(-1, 0), (0, 0.5), (1, 0), (-1, 1), (0, 1.5), (1, 1), (2, 0.5), (-2, 0.5)]
pair2 = label(pair1)
constant1 = J*(delta**(-0.5))
constant2 = J*(delta**(-0.5))
S0 = [entropy(T[0],E_init,0)]
C0 = [heat_capacity(beta,E_init,0)]
sus0 = [susceptibity(beta,E_init,0)]


index0 = {(-1, 0): {1: [1, 3]}, (0, 0.5): {1: [1, 2, 3, 4]}, (1, 0): {1: [3, 4]}, (-1, 1): {1: [1, 2]}, (0, 1.5): {1: [2]}, (1, 1): {1: [2, 4]}, (2, 0.5): {1: [4]}, (-2, 0.5): {1: [1]}}
dim_init = 8
t1 = time.time()
for i in range(1,M):
    print(i)
    dim = dim_init
    if dim < N_cut:
        H , index0, pair1 = H_N(pair1 , pair2 ,E_init , U, index0, delta )
        pair2 = label(pair1)
        E_init,U = diagnalization(H,delta)
        
        TN = T[i]
        sus0 .append(susceptibity(beta,E_init,i))
        S0.append(entropy(beta,E_init,i))
        C0.append(heat_capacity(beta,E_init,i))
    elif dim >= N_cut:
        E_init , pair1 ,U = truncation( N_cut, E_init , U,index0)
        pair2 = label(pair1)
        H , index0, pair1 = H_N(pair1 , pair2 ,E_init , U, index0 ,delta)
        E_init,U = diagnalization(H,delta)
        TN = T[i]
        sus0 .append(susceptibity(beta,E_init,i))
        C0.append(heat_capacity(beta,E_init,i))
        S0.append(entropy(beta,E_init,i))
    dim_init = 0
    for pair in E_init:#E = {(Q,S):[...],...}
        dim_init += len(E_init[pair])*(2*pair[1]+1)
    
E_init , pair1 ,U = truncation( N_cut, E_init , U,index0)
t2 = time.time()
print(t2-t1) 

