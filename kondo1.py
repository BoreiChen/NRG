# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 18:10:05 2019

@author: jeff
"""


def initH0(couple,delta,rho):
    """
    1,create initial hamiltonian in block form
    2,the pair denotes (Q,S)
    """
    
    couple1 = 4*rho*couple/(1+1/delta)
    c = 0.5*couple1*delta**(-0.5)
#    c = -2
    x = {(0,0):c*np.array([-3]), (0,1):c*np.array([1]),(1,0.5):c*np.array([0]),(-1,0.5):c*np.array([0])}
    return x

def anderson0(U,delta,D,gamma,ed):
    """
    1,create initial hamiltonian in block form
    2,the pair denotes (Q,S)
    """
    
    deltad = (2/(1+1/delta))*(ed + 0.5*U)/D
    U1 = (2/(1+1/delta))*U/(2*D)
    gamma1 = (2/(1+1/delta))*(2*gamma)/(np.pi*D)
    c = delta*(-0.5)
    x = {(-2,0):c*np.array([U1]),
         (-1,-0.5):c*np.array([[U1,gamma1**0.5],[gamma1**0.5,deltad]]),
         (-1,0.5):c*np.array([[U1,gamma1**0.5],[gamma1**0.5,deltad]]),
         (0,-1):c*np.array([deltad]),
         (0,0):c*np.array([[U1, -gamma1**0.5, gamma1**0.5,0],[-gamma1**0.5,deltad,0,0],
         [gamma1**0.5,0,deltad,0],[0,0,0,U1+2*deltad]]), 
         (0,1):c*np.array([deltad]),
         (1,0.5):c*np.array([[U1,-gamma1**0.5],[-gamma1**0.5,deltad]]),
         (1,-0.5):c*np.array([[U1,-gamma1**0.5],[-gamma1**0.5,deltad]]),
         (2,0):c*np.array([U1+2*deltad])}
    return x


def initH00():
    """
    1,create initial hamiltonian in block form
    2,the pair denotes (Q,S)
    """
    
   
    x = {(0,0):np.array([0]), (0,1):np.array([0]),(1,0.5):np.array([0]),(-1,0.5):np.array([0])}
    return x

def initH1(couple,delta):
    """
    1,create initial hamiltonian in block form
    2,the pair denotes (Q,S)
    """
    rho = 0.3
    couple1 = 4*rho*couple/(1+delta**0.5)
    c = 0.5*couple1*delta**(-0.5)
    x = {(-2,0.5):c*np.array([[0]]) , (-1,0):np.array([[-3,-1],[-1,0]]) , (-1,1):np.array([[1,1],[1,0]])
         ,(0,0.5):np.array([[0,-1,-1.5**0.5,0],[-1,-3,0,-0.5*2**0.5],[-1.5**0.5,0,1,-1.5**0.5],[0,-0.5*2**0.5,-1.5**0.5,0]])
         ,(0,1.5):np.array([[1]]), (1,0): np.array([[0,-0.5*2**0.5],[-0.5*2**0.5,-3]])
         ,(1,1):np.array([[0,-1],[-1,1]]), (2,0.5): np.array([[0]])}
        
    return x
def initU():
    u = {(0,0):np.array([[-1]]),(0,1):np.array([[1]]) 
        ,(1,0.5):np.array([[1]]) ,(-1,0.5):np.array([[1]])}
    return u
        