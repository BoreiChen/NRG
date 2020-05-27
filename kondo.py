# -*- coding: utf-8 -*-
"""
Created on Sun Jul 14 22:36:19 2019

@author: jeff
"""
"In this code I will show the NRG process with U(1)xU(1) symmetry"

import numpy as np
import math

class kondo:
    def __init__(self,delta = 1,couple = 1,dim = 1,matrix = np.eye(2),u1 = np.eye(8),u2 = np.eye(8)):
        "initialize the parameters"
        self.delta = delta
        self.couple = couple
        self.dim = dim
        self.matrix = matrix
        self.u1 = u1
        self.u2 = u2
        
        init_label = [(0,0),(0,1),(0,-1),(1,0.5),(-1,0.5),(1,-0.5),(-1,-0.5)]
    def initH1(self):
        """
        1,create initial hamiltonian in block form
        2,the pair denotes (Q,Sz)
        """
        c = 1
        x = ({(0,0):c*np.array([-3,1]) ,(0,1):c*np.array([1]) ,(0,-1):c*np.array([1])
                ,(1,0.5):np.array([0]) ,(-1,0.5):np.array([0]) ,(1,-0.5):np.array([0])
                ,(-1,-0.5):np.array([0])})
        
        return x
    def initU1(self):
            u = {(0,0):np.array([[-1/np.sqrt(2),1/np.sqrt(2)],[1/np.sqrt(2),1/np.sqrt(2)]]),(0,1):np.array([[1]]) ,(0,-1):np.array([[1]])
                ,(1,0.5):np.array([[1]]) ,(-1,0.5):np.array([[1]]) ,(1,-0.5):np.array([[1]])
                ,(-1,-0.5):np.array([[1]])}
            return u
        
        
def label(pair):
    "collect all posible pairs in each iteration"
    """
    input:
        pair:(Q,Sz)
        
    """
    y = []
    for i in pair:
        i_0, i_1 = i
        y.append((i_0-1,i_1))
        y.append((i_0,i_1+0.5))
        y.append((i_0,i_1-0.5))
        y.append((i_0+1,i_1))
        
    return y
    
    
    
def H_next( pair_old , pair_new ,E_old , U_old, order_old ):
    """
    this function is used to construct the hamiltonian , 
    parameters:
        
    pair_old:
        all possible combination for H(N+1)
        [[Q,Sz] ....]
    pair_new:
        [[Q,Sz,r,i]....]
    order:
        {(Q,Sz):{(r,i):1.....}......}
    N :the number of site
    energy: eigen energy of H(N)
    order_old : {(Q,Sz):{(r,i)=1.....}......}
    U_old: {(Q,Sz):[]}
        
    H_new = {(Q,Sz)= np.array([....]).....}
    """
    "firstly, I will perform mode counting to determime spin multipicity"
    H_new = {}
        
    E_new = {}  
    order = {}
    H_off = {}
    for pair in pair_new:
        E_new[pair] = {}
        order[pair] = {}
    "1, diagonal part"
    for j in pair_old:
        old_0 , old_1 = j
        for k in range(len(E_old[j])):
            E_new[(old_0-1 , old_1)][(k+1,1)] = E_old[j][k]
            E_new[(old_0 , old_1+0.5)][(k+1,2)] = E_old[j][k]
            E_new[(old_0 , old_1-0.5)][(k+1,3)] = E_old[j][k]        
            E_new[(old_0+1 , old_1)][(k+1,4)] = E_old[j][k]    
       #E_new = {(Q,Sz):{(k,r):E,....}......}
       #E_old = {(Q,Sz):[........],.....}
       #order = {(Q,Sz):{(k,r):1,...}...}
       #H_new = {(Q,Sz):[...],....}
    """
    asign the order
    """
    
    for pair in order:
        label_number = 0
        for pair2 in E_new[pair]:#(k,r)
            order[pair][pair2] = label_number
            label_number += 1
    for pair in pair_new:#(Q,Sz)
        H_new[pair] = np.zeros((len(order[pair]),len(order[pair])))
        
    for pair in order :
        for pair1 in order[pair]:
            H_new[pair][order[pair][pair1]][order[pair][pair1]] = E_new[pair][pair1]
        
        "2, off-diagonal part"
    
    
        #order_old : {(Q,Sz):{(r,i)=1.....}......}
    for pair in pair_new:#(Q,Sz)
        pair_0,pair_1 = pair 
        for pair1 in order[pair]:#{(k,r):1,....}
            pair1_0,pair1_1 = pair1
            for pair2 in order[pair]:#{(k,r):1,....}
                pair2_0,pair2_1 = pair2
                if (pair1_1,pair2_1) == (1,2):
                    for pair3 in order_old[(pair_0 +1,pair_1)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0,pair_1-0.5)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3_1,pair4_1) == (2,1):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0 +1,pair_1)][pair3]
                                b2 = order_old[(pair_0,pair_1-0.5)][pair4]
                                H_new[pair][a1][a2] +=( 
                                U_old[(pair_0 +1,pair_1)][pair1_0-1][b1]*U_old[(pair_0,pair_1-0.5)][pair2_0-1][b2])
                                H_new[pair][a2][a1] = H_new[pair][a1][a2]
                                
                if (pair1_1,pair2_1) == (3,4):
                    for pair3 in order_old[(pair_0,pair_1 +0.5)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0 -1,pair_1)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3[1],pair4[1]) == (4,3):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0,pair_1 +0.5)][pair3]
                                b2 = order_old[(pair_0 -1,pair_1)][pair4]
                                H_new[pair][a1][a2] +=( 
                                U_old[(pair_0,pair_1 +0.5)][pair1_0-1][b1]*U_old[(pair_0 -1,pair_1)][pair2_0-1][b2])
                                H_new[pair][a2][a1] = H_new[pair][a1][a2]
                                
                if (pair1_1,pair2_1) == (1,3):
                    for pair3 in order_old[(pair_0 +1,pair_1)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0,pair_1 +0.5)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3_1,pair4_1) == (3,1):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0 +1,pair_1)][pair3]
                                b2 = order_old[(pair_0,pair_1+ 0.5)][pair4]
                                H_new[pair][a1][a2] +=( 
                                U_old[(pair_0 +1,pair_1)][pair1_0-1][b1]*U_old[(pair_0,pair_1 +0.5)][pair2_0-1][b2])
                                H_new[pair][a2][a1] = H_new[pair][a1][a2]
                                
                if (pair1_1,pair2_1) == (2,4):
                    for pair3 in order_old[(pair_0,pair_1 -0.5)]:
                        pair3_0,pair3_1 = pair3
                        for pair4 in order_old[(pair_0 -1,pair_1)]:
                            pair4_0,pair4_1 = pair4
                            if (pair3_1,pair4_1) == (4,2):
                                a1 = order[pair][pair1]
                                a2 = order[pair][pair2]
                                b1 = order_old[(pair_0,pair_1 -0.5)][pair3]
                                b2 = order_old[(pair_0 -1,pair_1)][pair4]
                                H_new[pair][a1][a2] -=( 
                                U_old[(pair_0,pair_1 -0.5)][pair1_0-1][b1]*U_old[(pair_0-1,pair_1)][pair2_0-1][b2])
                                H_new[pair][a2][a1] = H_new[pair][a1][a2]
                                
    for pair in order:#(Q,Sz)
        for pair1 in order[pair]:
            pair1_0,pair1_1 = pair1
            for pair2 in order[pair]:
                pair2_0,pair2_1 = pair2
                if ((pair1_1,pair2_1) == (1,2) 
                or (pair1_1,pair2_1) == (1,3) or (pair1_1,pair2_1) == (2,4)):
                    a1 = order[pair][pair1]
                    a2 = order[pair][pair2]
                    H_new[pair][a1][a2] *= 1
                    H_new[pair][a1][a2] *= 1
                elif (pair1_1,pair2_1) == (3,4):
                    a1 = order[pair][pair1]
                    a2 = order[pair][pair2]
                    H_new[pair][a1][a2] *= -1
                    H_new[pair][a1][a2] *= -1
                
    pair_final = []
    for pair in order:
        pair_final.append(pair)
    return H_new  , pair_final, order

        
        
        
        
        
        
        
        
        
        
        