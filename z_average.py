# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:16:39 2020

@author: jeff
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import random

sus_imp1z = []
sus_imp2z = []
sus_imp3z = []
sus_imp4z = []
sus_imp5z = []
C_imp1z = []
C_imp2z = []
C_imp3z = []
C_imp4z = []
C_imp5z = []
S_imp1z = []
S_imp2z = []
S_imp3z = []
S_imp4z = []
S_imp5z = []
for i in range(1,len(T_odd)-1):
    sus_z1 = 0.5*(sus_imp1[i]+sus_imp1[i-1]+
                 (sus_imp1[i+1]-sus_imp1[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    sus_z2 = 0.5*(sus_imp2[i]+sus_imp2[i-1]+
                 (sus_imp2[i+1]-sus_imp2[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    sus_z3 = 0.5*(sus_imp3[i]+sus_imp3[i-1]+
                 (sus_imp3[i+1]-sus_imp3[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    sus_z4 = 0.5*(sus_imp4[i]+sus_imp4[i-1]+
                 (sus_imp4[i+1]-sus_imp4[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    sus_z5 = 0.5*(sus_imp5[i]+sus_imp5[i-1]+
                 (sus_imp5[i+1]-sus_imp5[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    c_z1 = 0.5*(C_imp1[i]+C_imp1[i-1]+
                 (C_imp1[i+1]-C_imp1[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    c_z2 = 0.5*(C_imp2[i]+C_imp2[i-1]+
                 (C_imp2[i+1]-C_imp2[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    c_z3 = 0.5*(C_imp3[i]+C_imp3[i-1]+
                 (C_imp3[i+1]-C_imp3[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    c_z4 = 0.5*(C_imp4[i]+C_imp4[i-1]+
                 (C_imp4[i+1]-C_imp4[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    c_z5 = 0.5*(C_imp5[i]+C_imp5[i-1]+
                 (C_imp5[i+1]-C_imp5[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    s_z1 = 0.5*(S_imp1[i]+S_imp1[i-1]+
                 (S_imp1[i+1]-S_imp1[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    s_z2 = 0.5*(S_imp2[i]+S_imp2[i-1]+
                 (S_imp2[i+1]-S_imp2[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    s_z3 = 0.5*(S_imp3[i]+S_imp3[i-1]+
                 (S_imp3[i+1]-S_imp3[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    s_z4 = 0.5*(S_imp4[i]+S_imp4[i-1]+
                 (S_imp4[i+1]-S_imp4[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    s_z5 = 0.5*(S_imp5[i]+S_imp5[i-1]+
                 (S_imp5[i+1]-S_imp5[i-1])*(T_odd[i]-T_odd[i-1])/(T_odd[i+1]-T_odd[i-1]))
    sus_imp1z.append(sus_z1)
    sus_imp2z.append(sus_z2)
    sus_imp3z.append(sus_z3)
    sus_imp4z.append(sus_z4)
    sus_imp5z.append(sus_z5)
    C_imp1z.append(c_z1)
    C_imp2z.append(c_z2)
    C_imp3z.append(c_z3)
    C_imp4z.append(c_z4)
    C_imp5z.append(c_z5)
    S_imp1z.append(s_z1)
    S_imp2z.append(s_z2)
    S_imp3z.append(s_z3)
    S_imp4z.append(s_z4)
    S_imp5z.append(s_z5)
Tn = [T_odd[i] for i in range(1,len(T_odd)-1)]












