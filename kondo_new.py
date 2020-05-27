# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 23:49:52 2019

@author: jeff
"""
import numpy as np
import math
from numba import jit
from numba import vectorize, int64
def same(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    iset = set1.intersection(set2)
    return list(iset)
  
def label(pair):
    "collect all posible pairs in each iteration"
    """
    input:
        pair:(Q,S)
    """
    y = []
    for i in pair:
        i_0, i_1 = i
        y.append((i_0-1,i_1))
        y.append((i_0,i_1+0.5))
        if i_1 >= 0.5:
            y.append((i_0,i_1-0.5))
        y.append((i_0+1,i_1))
    re_y = []
    for pairs in y:
        if pairs not in re_y:
            re_y.append(pairs)
    return re_y
    


def H_N(pair_old , pair_new ,E_old , U_old, index_old,delta ):
    """
    this function is used to construct the hamiltonian , 
    parameters:
        
    pair_old:
        all possible combination for H(N+1)
        [[Q,S] ....]
    pair_new:
        [[Q,S,r,i]....]
    order:
        {(Q,S):{(r,i):1.....}......}
    N :the number of site
    energy: eigen energy of H(N)
    order_old : {(Q,S):{(r,i)=1.....}......}
    U_old: {(Q,S):[]}
        
    H_new = {(Q,S)= np.array([....]).....}
    """
    
    H_new = {}
    E_new = {}  
    index = {}
    for pair in pair_new:
        E_new[pair] = {}
        H_new[pair] = {}
        index[pair] = {}
    "1, diagonal part"
    for j in pair_old:
        old_0,old_1 = j
        for l in range(len(E_old[j])):
            index[(old_0-1,old_1)][l+1] = []
            index[(old_0,old_1+0.5)][l+1] = []
            E_new[(old_0-1,old_1)][l+1] = []
            E_new[(old_0,old_1+0.5)][l+1] = []
            if old_1 >= 0.5:
                index[(old_0,old_1-0.5)][l+1] = []
                E_new[(old_0,old_1-0.5)][l+1] = []
            index[(old_0+1,old_1)][l+1] = []
            E_new[(old_0+1,old_1)][l+1] = []
    for j in pair_old:
        old_0,old_1 = j
        for k in range(len(E_old[j])):
            E_new[(old_0-1,old_1)][k+1].append(E_old[j][k])
            index[(old_0-1,old_1)][k+1].append(1)
            E_new[(old_0,old_1+0.5)][k+1].append(E_old[j][k])
            index[(old_0,old_1+0.5)][k+1].append(2)
            if old_1 >= 0.5:
                E_new[(old_0 , old_1-0.5)][k+1].append(E_old[j][k])
                index[(old_0,old_1-0.5)][k+1].append(3)
            E_new[(old_0+1 , old_1)][k+1].append(E_old[j][k])
            index[(old_0+1 , old_1)][k+1].append(4)
       #E_new = {(Q,Sz):{(r,i):E,....}......}
       #E_old = {(Q,Sz):[........],.....}
       #order = {(Q,Sz):{(r,i):1,...}...}
       #H_new = {(Q,Sz):[...],....}
    index1 = index
    for i in index:
        for j in index[i]:
            index[i][j] = sorted(index[i][j])
    
    
    for i in index:
        dim_H = 0
        for j in index[i]:
            dim_H += len(index[i][j])
        H_new[i] = np.zeros((dim_H,dim_H))
    
    "fill in the diagonal element"
    for pair in E_new:
        dim = 0
        for ind in E_new[pair]:
            if 1 in index[pair][ind]:
                e = index1[pair][ind]
                H_new[pair][dim + e.index(1)][dim + e.index(1)] =  E_old[(pair[0]+1,pair[1])][ind-1]*delta**0.5
            if 2 in index[pair][ind]:
                e = index1[pair][ind]
                H_new[pair][dim + e.index(2)][dim + e.index(2)] =  E_old[(pair[0],pair[1]-0.5)][ind-1]*delta**0.5
            if 3 in index[pair][ind]:
                e = index1[pair][ind]
                H_new[pair][dim + e.index(3)][dim + e.index(3)] =  E_old[(pair[0],pair[1]+0.5)][ind-1]*delta**0.5
            if 4 in index[pair][ind]:
                e = index1[pair][ind]
                H_new[pair][dim + e.index(4)][dim + e.index(4)] =  E_old[(pair[0]-1,pair[1])][ind-1]*delta**0.5
            dim += len(E_new[pair][ind])
    
    "2, off diagonal part"
    
    for pair in E_new:
        dim1 = 0
        for r1 in index[pair]:
            dim2 = 0
            for r2 in index[pair]:
                if 1 in index[pair][r1] and 2 in index[pair][r2]:
                    i1 = index_old[(pair[0]+1,pair[1])]
                    i2 = index_old[(pair[0],pair[1]-0.5)]
                    V12 = 0
                    dim_r31 = 0
                    dim_r32 = 0
                    i12 = same(list(i1), list(i2))
                    
                    for r3 in i12:
                        if 2 in i1[r3] and 1 in i2[r3] :
                            V12 += (U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(2)]*
                                      U_old[(pair[0],pair[1]-0.5)][r2-1][dim_r32 + i2[r3].index(1)])
                        if  4 in i1[r3] and 3 in i2[r3]:    
                            V12 += ((2*pair[1]/(2*pair[1] + 1))**0.5*
                                   U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(4)]*
                                      U_old[(pair[0],pair[1]-0.5)][r2-1][dim_r32 + i2[r3].index(3)]) 
#                        if 3 in i1[r3] and 1 in i2[r3] :
#                            V12 += (U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(3)]*
#                                      U_old[(pair[0],pair[1]-0.5)][r2-1][dim_r32 + i2[r3].index(1)])
#                        if 4 in i1[r3] and 2 in i2[r3]:       
#                            V12   -= ((2*pair[1]/(2*pair[1] - 1))**0.5*
#                                   U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(4)]*
#                                      U_old[(pair[0],pair[1]-0.5)][r2-1][dim_r32 + i2[r3].index(2)]) 
                        dim_r31 += len(i1[r3])
                        dim_r32 += len(i2[r3])
                    s1 = index[pair][r1].index(1)
                    s2 = index[pair][r2].index(2)
                    H_new[pair][dim1 + s1][dim2 + s2] = V12
                    H_new[pair][dim2 + s2][dim1 + s1] = np.conjugate(V12)
                if 1 in index[pair][r1] and 3 in index[pair][r2]:
                    i1 = index_old[(pair[0]+1,pair[1])]
                    i2 = index_old[(pair[0],pair[1]+0.5)]
                    V13 = 0
                    dim_r31 = 0
                    dim_r32 = 0
                    i12 = same(list(i1), list(i2))
                    for r3 in i12:
#                        print(i1)
#                        print(i2)
#                        if 2 in i1[r3] and 1 in i2[r3]:
#                            print('yes')
#                            V13 +=  (U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(2)]*
#                                      U_old[(pair[0],pair[1]+0.5)][r2-1][dim_r32 + i2[r3].index(1)])
#                        if 4 in i1[r3] and 3 in i2[r3]:    
#                            V13 += (((2*pair[1] + 2)/(2*pair[1] + 3))**0.5*
#                                   U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(4)]*
#                                      U_old[(pair[0],pair[1]+0.5)][r2-1][dim_r32 + i2[r3].index(3)]) 
                        if 3 in i1[r3] and 1 in i2[r3]:
                            
                            V13 +=  (U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(3)]*
                                      U_old[(pair[0],pair[1]+0.5)][r2-1][dim_r32 + i2[r3].index(1)])
                        if 4 in i1[r3] and 2 in i2[r3]:          
                            V13 -= (((2*pair[1] + 2)/(2*pair[1] + 1))**0.5*
                                   U_old[(pair[0]+1,pair[1])][r1-1][dim_r31 + i1[r3].index(4)]*
                                      U_old[(pair[0],pair[1]+0.5)][r2-1][dim_r32 + i2[r3].index(2)])
                        dim_r31 += len(i1[r3])
                        dim_r32 += len(i2[r3])
                    s1 = index[pair][r1].index(1)
                    s3 = index[pair][r2].index(3)
                    H_new[pair][dim1 + s1][dim2 + s3] = V13
                    H_new[pair][dim2 + s3][dim1 + s1] = np.conjugate(V13)
                if 2 in index[pair][r1] and 4 in index[pair][r2]:
                    i1 = index_old[(pair[0],pair[1]-0.5)]
                    i2 = index_old[(pair[0]-1,pair[1])]
                    V24 = 0
                    dim_r31 = 0
                    dim_r32 = 0
                    i12 = same(list(i1), list(i2))
                    for r3 in i12:
#                        if 2 in i1[r3] and 1 in i2[r3]:                          
#                            V24 +=  (U_old[(pair[0],pair[1]-0.5)][r1-1][dim_r31 + i1[r3].index(2)]*
#                                      U_old[(pair[0]-1,pair[1])][r2-1][dim_r32 + i2[r3].index(1)])
#                        if 4 in i1[r3] and 3 in i2[r3]: 
#                            V24 += (((2*pair[1]+1)/(2*pair[1] + 2))**0.5*
#                                   U_old[(pair[0],pair[1]-0.5)][r1-1][dim_r31 + i1[r3].index(4)]*
#                                      U_old[(pair[0]-1 ,pair[1])][r2-1][dim_r32 + i2[r3].index(3)]) 
                        if 3 in i1[r3] and 1 in i2[r3]:  
                            V24 += (U_old[(pair[0],pair[1]-0.5)][r1-1][dim_r31 + i1[r3].index(3)]*
                                      U_old[(pair[0]-1,pair[1])][r2-1][dim_r32 + i2[r3].index(1)])
                        if 4 in i1[r3] and 2 in i2[r3]:
                            V24 -= (((2*pair[1]+1)/(2*pair[1]))**0.5*
                                   U_old[(pair[0],pair[1]-0.5)][r1-1][dim_r31 + i1[r3].index(4)]*
                                      U_old[(pair[0]-1,pair[1])][r2-1][dim_r32 + i2[r3].index(2)])
                        dim_r31 += len(i1[r3])
                        dim_r32 += len(i2[r3])
                    s2 = index[pair][r1].index(2)
                    s4 = index[pair][r2].index(4)
                    H_new[pair][dim1 + s2][dim2 + s4] = V24*(2*pair[1]/(2*pair[1]+1))**0.5
                    H_new[pair][dim2 + s4][dim1 + s2] = np.conjugate(V24*(2*pair[1]/(2*pair[1]+1))**0.5)
                if 3 in index[pair][r1] and 4 in index[pair][r2]:
                    i1 = index_old[(pair[0],pair[1]+0.5)]
                    i2 = index_old[(pair[0]-1,pair[1])]
                    dim_r31 = 0
                    dim_r32 = 0
                    V34 = 0
                    i12 = same(list(i1), list(i2))
                    for r3 in i12:
          
                        if 2 in i1[r3] and 1 in i2[r3]:  
                            V34 +=  (U_old[(pair[0],pair[1]+0.5)][r1-1][dim_r31 + i1[r3].index(2)]*
                                      U_old[(pair[0]-1,pair[1])][r2-1][dim_r32 + i2[r3].index(1)])
                        if 4 in i1[r3] and 3 in i2[r3]: 
                            V34 += (((2*pair[1]+1)/(2*pair[1] + 2))**0.5*
                                   U_old[(pair[0],pair[1]+0.5)][r1-1][dim_r31 + i1[r3].index(4)]*
                                      U_old[(pair[0]-1,pair[1])][r2-1][dim_r32 + i2[r3].index(3)] )
                            
#                        if 3 in i1[r3] and 1 in i2[r3]: 
#                            print('yes')
#                            V34 += ( U_old[(pair[0],pair[1]+0.5)][r1-1][dim_r31 + i1[r3].index(3)]*
#                                      U_old[(pair[0]-1 ,pair[1])][r2-1][dim_r32 + i2[r3].index(1)])
#                        if 4 in i1[r3] and 2 in i2[r3]:
#                            V34 -= (((2*pair[1]+1)/(2*pair[1]))**0.5*
#                                   U_old[(pair[0],pair[1]+0.5)][r1-1][dim_r31 + i1[r3].index(4)]*
#                                      U_old[(pair[0]-1,pair[1])][r2-1][dim_r32 + i2[r3].index(2)] )
                        dim_r31 += len(i1[r3])
                        dim_r32 += len(i2[r3])
                        
                    s3 = index[pair][r1].index(3)
                    s4 = index[pair][r2].index(4)
                    H_new[pair][dim1 + s3][dim2 + s4] = -V34*((2*pair[1]+2)/(2*pair[1]+1))**0.5
                    H_new[pair][dim2 + s4][dim1 + s3] = np.conjugate(-V34*((2*pair[1]+2)/(2*pair[1]+1))**0.5)
                dim2 += len(index[pair][r2])
            dim1 += len(index[pair][r1])
                
    return H_new , index ,pair_new
            
    

    
    
    
    
    




