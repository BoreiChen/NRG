# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 22:38:07 2019

@author: jeff
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import time
from numba import jit
from numba import vectorize , int64
M = 
J= 0.04
D = 1
rho = 0.5
delta = 4
cc = 4*rho/(1 + 1/delta)
beta = 0.5


#J = J/cc
#init_H =anderson0(U,delta,D,gamma,ed)
init_H = initH0(J,delta,rho)
init_U = initU()
N_cut = 100
#
#Tk = D*(abs(J)*rho)**0.5*np.exp(-1/(2*abs(J)*rho))
def energy_flow(E,N):
    E_flow = {(0,0):[],(1,0.5):[],(-1,0.5):[],(0,1):[]}
    for i in range(N):
        if i%2 == 0:
            E_flow[(0,0)].append(E[(0,0)])
    


def entropy(b,E,N):
    e_avg = 0
    partition = 0
    for pair in E:
        q,s = pair
        for e in E[pair]:
            e_avg += (2*s+1)*(e)*np.exp(-e*b)
            partition += (2*s+1)*np.exp(-e*b)
    S = (e_avg*beta/partition) + np.log(partition)
    return S

def heat_capacity(b,E,N):
    esqrt_avg = 0
    e_avg = 0
    partition = 0
    for pair in E:
        q,s = pair
        for e in E[pair]:
            esqrt_avg += (2*s+1)*((e)**2)*np.exp(-e*b)
            e_avg += (2*s+1)*(e)*np.exp(-e*b)
            partition += (2*s+1)*np.exp(-e*b)
    
    C = b**2*((esqrt_avg/partition) - (e_avg/partition)**2)
    return C
    

def susceptibity(b,E,N):
    msqrt_avg = 0
    m_avg = 0
    partition = 0
    for pair in E:
        q,s = pair
        for e in E[pair]:
            msqrt_avg += ((2*s+1)*((2*s+1)**2 -1)/12)*np.exp(-e*b)
            m_avg += 0
            partition += (2*s+1)*np.exp(-e*b)
    M_total = ((msqrt_avg - m_avg**2)/partition)
    return M_total


e1 = []
E = init_H
U = init_U

pair1 = [(0,0),(0,1),(1,0.5),(-1,0.5)]
pair2 = label(pair1)
constant1 = J*(delta**(-0.5))
constant2 = J*(delta**(-0.5))
C = []
chi = []
T = []
index0 = {(0,0):{1:[3]},(0,1):{1:[2]},(1,0.5):{1:[4]}
                                  ,(-1,0.5):{1:[1]}}
#append susepbility
S = []
sus = []
T = []
lnT = []
#append energy flow
e_flow = []
n_flow = []
dim_init = 8
t1 = time.time()
for i in range(M):
#    t3 = time.time()
    print(i)
    dim = dim_init
    if dim < N_cut:
        H , index0, pair1 = H_N(pair1 , pair2 ,E , U, index0,delta )
        pair2 = label(pair1)
        
        E,U = diagnalization(H,delta)
        TN = D*0.5*(1 + 1/delta)*delta**(-0.5*(i-1))/beta
        T.append(TN)
        lnT.append(np.log10(TN))
        sus .append(susceptibity(beta,E,i))
        S.append(entropy(beta,E,i))
        C.append(heat_capacity(beta,E,i))
#        e1 = []
#        for pair in E:
#            for k in E[pair]:
#                e1.append(k)
#        e1 = sorted(e1)
#        for j in range(10):
#            e_flow.append(e1[j])
#            n_flow.append(i)

    elif dim >= N_cut:
        E , pair1 ,U = truncation( N_cut, E , U,index0)
        pair2 = label(pair1)
        H , index0, pair1 = H_N(pair1 , pair2 ,E , U, index0 ,delta)
        E,U = diagnalization(H,delta)
        #construct energy flow diagram
#        e1 = []
#        for pair in E:
#            for k in E[pair]:
#                e1.append(k)
#        e1 = sorted(e1)
#        for j in range(10):
#            e_flow.append(e1[j])
#            n_flow.append(i)
        
        TN = D*0.5*(1 + 1/delta)*delta**(-0.5*(i-1))/beta
        T.append(TN)
        lnT.append(np.log10(TN))
        sus .append(susceptibity(beta,E,i))
        S.append(entropy(beta,E,i))
        C.append(heat_capacity(beta,E,i))
    dim_init = 0
    for pair in E:#E = {(Q,S):[...],...}
        dim_init += len(E[pair])*(2*pair[1]+1)
#    t4 = time.time()
#    print(t4-t3)
E , pair1 ,U = truncation( N_cut, E , U,index0)    
t2 = time.time()
print(t2-t1) 

    
    
