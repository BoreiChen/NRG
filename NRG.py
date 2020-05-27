# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 19:55:38 2019

@author: jeff
"""
import numpy as np
import math
from scipy import linalg
from scipy.linalg import *

def diagnalization(matrix,delta):
    #delta coupling*delta
    E = {}
    U = {}
    e_collect = []
    for pair in matrix:
        eigvalue,eigvector = eig(matrix[pair])
        e_collect.append(min(eigvalue))
        E[pair] = eigvalue.real
        U[pair] = np.transpose(eigvector).real
    E_min = min(e_collect)
    for pair in E:
#        #E[pair] = [E[pair][i] - E_min for i in range(len(E[pair]))]
        E[pair] =  [(E[pair][i]-E_min).real for i in range(len(E[pair]))]
    return E,U
      
def truncation( N_cut, e_value, U, index):
    #U = {(Q,Sz):[...],...}
    index_r = {}
    e_value1 = {}
    for pair in e_value:
        pair_0, pair_1 = pair
        r = 1
        dim = 0
        dimm = len(index[pair][r])
        for i in range(len(e_value[pair])):
            e_value1[(pair_0,pair_1,i)] = e_value[pair][i]
            index_r[(pair_0,pair_1,i)] = (r,index[pair][r][i-dim])
            if i < dimm-1:
                continue
            else:
                dim += len(index[pair][r])
                r += 1
                if r in index[pair]:
                    dimm += len(index[pair][r])
    e_value_new = sorted(e_value1.items(),key=lambda x:x[1])
    e_value_new = [pair[0] for pair in e_value_new]
    """
        starting truncate
    """
    """
        parameter:
            
    """
    index_rr = {}
    energy = {}#contain all energies after truncation
    "ordering labels (Q,S,i)"
    num = 0
    ordering = 0
    while num < N_cut:
        num += 2*e_value_new[ordering][1]+1
        energy[e_value_new[ordering]] = e_value1[e_value_new[ordering]]
        ordering += 1
    #energy = {(Q,Sz,r): E,.... }
    pairs_new = []
    for pair in energy:
        pair_0, pair_1,pair_2 = pair
        pairs_new.append((pair_0,pair_1))
    #pair_new = [(Q,Sz),....]--------repeated
    label_new = []
    for i in pairs_new:
        if not i in label_new:
            label_new.append(i)
    #new index
    for pair in energy:
        index_rr[(pair[0],pair[1])] = {}
    for pair in energy:
        if index_r[pair][0] in index_rr[(pair[0],pair[1])]:
            index_rr[(pair[0],pair[1])][index_r[pair][0]].append(index_r[pair][1])
        else:
            index_rr[(pair[0],pair[1])][index_r[pair][0]] = [index_r[pair][1]]
    
    #label_new = [(Q,Sz),...]--------remove repeated elements
    U_new = {}
    energy_new = {}
    for pair in label_new:#(Q,Sz)
        pair_0, pair_1 = pair
        energy_new[pair] = []
        U_new[pair] = []
        for pair1 in energy:
            pair1_0, pair1_1, pair1_2 = pair1
            if pair_0 == pair1_0 and pair_1 == pair1_1:
                energy_new[pair].append(energy[pair1])
                U_new[pair].append(U[pair][pair1_2])
    #U_new = {[Q,Sz]:[.......] }
    #energy_new = {[Q,Sz]:[.....]}
    return energy_new, label_new ,U_new 
        
        
    
        
        
        
        