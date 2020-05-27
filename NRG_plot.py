# -*- coding: utf-8 -*-
"""
Created on Sun Aug 18 10:45:05 2019

@author: jeff
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from kondo import kondo


plt.figure(figsize=(9,6))
plt.title(r'Suscepbility ',fontsize = 24)
plt.plot(T_odd,sus_imp2,'o',label = 'J = -1')
plt.plot(T_odd,sus_imp3,'o',label = 'J = -3')
plt.plot(T_odd,sus_imp1,'o',label = 'J = -5')

plt.legend(fontsize = 16)
plt.xlabel('lnT',fontsize = 16)
plt.ylabel('$\chi_{imp}$',fontsize = 16)



