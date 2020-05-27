# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 01:51:23 2019

@author: jeff
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from kondo import kondo
#plt.plot(T,c0)
#plt.plot(T,c4)
#plt.plot(T,c1)
#plt.plot(T,c2)
#plt.plot(T,c3)

#susceptibility
#plt.plot(Tn,sus_imp1z)
#plt.plot(Tn,sus_imp2z)
#plt.plot(Tn,sus_imp3z)
#plt.plot(Tn,sus_imp4z)
#plt.plot(Tn,sus_imp5z)

#specific heat
#plt.plot(Tn,C_imp1z)
#plt.plot(Tn,C_imp2z)
#plt.plot(Tn,C_imp3z)
#plt.plot(Tn,C_imp4z)
#plt.plot(Tn,C_imp5z)
#entropy
#plt.plot(Tn,S_imp1z)
#plt.plot(Tn,S_imp2z)
#plt.plot(Tn,S_imp3z)
#plt.plot(Tn,S_imp4z)
#plt.plot(Tn,S_imp5z)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.figure(figsize=(9,6))
plt.title(r'Susceptibility ',fontsize = 24)
plt.plot(Tn,sus_imp1z,'o-',label = 'J = 0.05')
plt.plot(Tn,sus_imp2z,'o-',label = 'J = 0.03')
plt.plot(Tn,sus_imp3z,'o-',label = 'J = 0.07')
plt.plot(Tn,sus_imp4z,'o-',label = 'J = 0.06')
plt.plot(Tn,sus_imp5z,'o-',label = 'J = 0.04')
plt.axis([-25, 0 , 0, 0.3])
plt.legend(fontsize = 16)
plt.xlabel('lnT',fontsize = 16)
plt.ylabel(r'$\chi_imp$',fontsize = 16)

plt.figure(figsize=(9,6))
plt.title(r'Specific heat  ',fontsize = 24)
plt.plot(Tn,C_imp1z,'o-',label = 'J = 0.05')
plt.plot(Tn,C_imp2z,'o-',label = 'J = 0.03')
plt.plot(Tn,C_imp3z,'o-',label = 'J = 0.07')
plt.plot(Tn,C_imp4z,'o-',label = 'J = 0.06')
plt.plot(Tn,C_imp5z,'o-',label = 'J = 0.04')
plt.axis([-25, 0 , 0, 0.25])
plt.legend(fontsize = 16)
plt.xlabel('lnT',fontsize = 16)
plt.ylabel(r'$C_imp$',fontsize = 16)

plt.figure(figsize=(9,6))
plt.title(r'Entropy  ',fontsize = 24)
plt.plot(Tn,S_imp1z,'o-',label = 'J = 0.05')
plt.plot(Tn,S_imp2z,'o-',label = 'J = 0.03')
plt.plot(Tn,S_imp3z,'o-',label = 'J = 0.07')
plt.plot(Tn,S_imp4z,'o-',label = 'J = 0.06')
plt.plot(Tn,S_imp5z,'o-',label = 'J = 0.04')
plt.legend(fontsize = 16)
plt.axis([-25, 0 , 0, 1.2])
plt.xlabel('lnT',fontsize = 16)
plt.ylabel(r'$S_imp/ln2$',fontsize = 16)

plt.figure(figsize=(9,6))
plt.title(r'J = 0.05  ',fontsize = 24)
plt.plot(Tn,S_imp1z,'o-',label = 'entropy')
plt.plot(Tn,C_imp1z,'o-',label = 'specific heat')
plt.plot(Tn,sus_imp1z,'o-',label = 'susceptibility')
s = [0 + i*0.1 for i in range(11)]
t = [-12.546349804879155 for i in range(len(s))]
plt.plot(t,s,'--,r')
plt.legend(fontsize = 16)
plt.text(-14, 0.8, r'$T_{k}$',color='red',fontsize = 27)
plt.xlabel('lnT',fontsize = 16)
plt.ylabel(r'magnititude',fontsize = 16)

