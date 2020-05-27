# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:35:35 2020

@author: jeff
"""

sus_imp = []
C_imp = []
T_odd = []
S_imp = []
wilson = []
for i in range(2,len(sus)):
    if i > 0:
        sus_imp.append(sus[i] - sus0[i]+0.25)
        C_imp.append(C[i] - C0[i])
        S_imp.append((S[i]-S0[i])/(np.log(2))+1)
        T_odd.append(lnT[i])
#for i in range(len(sus_imp)):
#    wilson.append((4*np.pi**2/3)*(sus_imp[i]/C_imp[i]))
#for i in range(2,len(sus)):
#    sus_imp.append(sus[i] - sus0[i] +0.25)
#    C_imp.append(C[i] - C0[i])
#    T_odd.append(lnT[i])
    
    
plt.plot(T_odd,sus_imp)
plt.plot(T_odd,S_imp)
plt.plot(T_odd,C_imp)