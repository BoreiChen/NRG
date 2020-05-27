# -*- coding: utf-8 -*-
"""
Created on Sat Aug 24 09:33:28 2019

@author: jeff
"""
"In this code I will show the NRG process with SU(2)xU(1) symmetry"
import numpy as np
import math
from numba import jit

        
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
    
    
@jit
def H_next( pair_old , pair_new ,E_old , U_old, order_old ):
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
    order = {}
    for pair in pair_new:
        E_new[pair] = {}
        order[pair] = {}
    "1, diagonal part"
    for j in pair_old:
        old_0,old_1 = j
        for k in range(len(E_old[j])):
            E_new[(old_0-1,old_1)][(k+1,1)] = E_old[j][k]
            E_new[(old_0,old_1+0.5)][(k+1,2)] = E_old[j][k]
            if old_1 >= 0.5:
                E_new[(old_0 , old_1-0.5)][(k+1,3)] = E_old[j][k]        
            E_new[(old_0+1 , old_1)][(k+1,4)] = E_old[j][k]    
       #E_new = {(Q,Sz):{(r,i):E,....}......}
       #E_old = {(Q,Sz):[........],.....}
       #order = {(Q,Sz):{(r,i):1,...}...}
       #H_new = {(Q,Sz):[...],....}
    """
    asign the order
    """
    
    for pair in order:
        label_number = 0
        for pair2 in E_new[pair]:#(r,1)
            order[pair][pair2] = label_number
            label_number += 1
    for pair in pair_new:#(Q,S)
        H_new[pair] = np.zeros((len(order[pair]),len(order[pair])))
        
    for pair in order :
        for pair1 in order[pair]:
            H_new[pair][order[pair][pair1]][order[pair][pair1]] = E_new[pair][pair1]
    H_diag = H_new
    "2, off-diagonal part"
    
    
        #order_old : {(Q,S):{(r,i)=1.....}......}
    for pair in pair_new:#(Q,S)
        pair_0,pair_1 = pair 
        for pair1 in order[pair]:#{(k,r):1,....}
            pair1_0,pair1_1 = pair1
            for pair2 in order[pair]:#{(k,r):1,....}
                pair2_0,pair2_1 = pair2
                if (pair1_1,pair2_1) == (1,2):
                    rlist1 = [*order_old[(pair_0 +1,pair_1)]]
                    rlist1 = [rlist1[i][0] for i in range(len(rlist1))]
                    rlist2 = [*order_old[(pair_0 ,pair_1-0.5)]]
                    rlist2 = [rlist2[i][0] for i in range(len(rlist2))]
                    rnew = [l for l in rlist1 if l in rlist2]
                    for pair3 in order_old[(pair_0 +1,pair_1)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0,pair_1-0.5)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3_1,pair4_1) == (3,1):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0 +1,pair_1)][pair3]
                                b2 = order_old[(pair_0,pair_1-0.5)][pair4]
                                H_new[pair][a1][a2] +=( 
                                U_old[(pair_0 +1,pair_1)][pair1_0-1][b1]*U_old[(pair_0,pair_1-0.5)][pair2_0-1][b2])
                                
                            if (pair3_1,pair4_1) == (4,2):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0 +1,pair_1)][pair3]
                                b2 = order_old[(pair_0,pair_1-0.5)][pair4]
                                H_new[pair][a1][a2] -=( np.sqrt((2*pair_1)/(2*pair_1 +1))*U_old[(pair_0 +1,pair_1)][pair1_0-1][b1]*U_old[(pair_0,pair_1-0.5)][pair2_0-1][b2])
                                
                if (pair1_1,pair2_1) == (3,4):
                    
                    
                    
                    
                    for pair3 in order_old[(pair_0,pair_1 +0.5)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0 -1,pair_1)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3[1],pair4[1]) == (3,1):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0,pair_1 +0.5)][pair3]
                                b2 = order_old[(pair_0 -1,pair_1)][pair4]
                                H_new[pair][a1][a2] +=( 
                                U_old[(pair_0,pair_1 +0.5)][pair1_0-1][b1]*U_old[(pair_0 -1,pair_1)][pair2_0-1][b2])
                                
                            if (pair3[1],pair4[1]) == (4,2):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0,pair_1 +0.5)][pair3]
                                b2 = order_old[(pair_0 -1,pair_1)][pair4]
                                H_new[pair][a1][a2] -=( np.sqrt((2*pair_1+1)/(2*pair_1 +2))*U_old[(pair_0,pair_1 +0.5)][pair1_0-1][b1]*U_old[(pair_0 -1,pair_1)][pair2_0-1][b2])
                                
                if (pair1_1,pair2_1) == (1,3):
                    
                    
                    
                    
                    
                    
                    for pair3 in order_old[(pair_0 +1,pair_1)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0,pair_1 +0.5)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3_1,pair4_1) == (2,1):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0 +1,pair_1)][pair3]
                                b2 = order_old[(pair_0,pair_1+ 0.5)][pair4]
                                H_new[pair][a1][a2] +=( 
                                U_old[(pair_0 +1,pair_1)][pair1_0-1][b1]*U_old[(pair_0,pair_1 +0.5)][pair2_0-1][b2])
                                
                            if (pair3_1,pair4_1) == (4,3):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0 +1,pair_1)][pair3]
                                b2 = order_old[(pair_0,pair_1+ 0.5)][pair4]
                                H_new[pair][a1][a2] +=( np.sqrt((2*pair_1+2)/(2*pair_1 +1))*U_old[(pair_0 +1,pair_1)][pair1_0-1][b1]*U_old[(pair_0,pair_1 +0.5)][pair2_0-1][b2])
                                
                if (pair1_1,pair2_1) == (2,4):
                    
                    
                    
                    
                    
                    
                    for pair3 in order_old[(pair_0,pair_1 -0.5)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0 -1,pair_1)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3_1,pair4_1) == (2,1):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0,pair_1 -0.5)][pair3]
                                b2 = order_old[(pair_0 -1,pair_1)][pair4]
                                H_new[pair][a1][a2] +=( 
                                U_old[(pair_0,pair_1 -0.5)][pair1_0-1][b1]*U_old[(pair_0-1,pair_1)][pair2_0-1][b2])
                                
                            if (pair3_1,pair4_1) == (4,3):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0,pair_1 -0.5)][pair3]
                                b2 = order_old[(pair_0 -1,pair_1)][pair4]
                                H_new[pair][a1][a2] +=( np.sqrt((2*pair_1 +1)/(2*pair_1 ))*U_old[(pair_0,pair_1 -0.5)][pair1_0-1][b1]*U_old[(pair_0-1,pair_1)][pair2_0-1][b2])                     
    for pair in order:#(Q,S)
        pair_0,pair_1 = pair
        for pair1 in order[pair]:
            pair1_0,pair1_1 = pair1
            for pair2 in order[pair]:
                pair2_0,pair2_1 = pair2
                if (pair1_1,pair2_1) == (4,2):
                    a1 = order[pair][pair1]
                    a2 = order[pair][pair2]
                    H_new[pair][a1][a2] *= np.sqrt(2*pair_1/(2*pair_1+1))
                    
                elif (pair1_1,pair2_1) == (4,3):
                    a1 = order[pair][pair1]
                    a2 = order[pair][pair2]
                    H_new[pair][a1][a2] *= -np.sqrt(2*pair_1/(2*pair_1+1))
    for pair in H_new:
        H_trans = np.transpose(np.conjugate(H_new[pair]))
        H_new[pair] = (H_new[pair] +  H_trans) - H_diag[pair]
                    
                
    pair_final = []
    for pair in order:
        pair_final.append(pair)
    
    return H_new  , pair_final, order











