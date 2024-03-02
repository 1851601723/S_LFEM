# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:24:49 2023

@author: XinLi
"""
import BendCruvedHMSRod3 as inp1
from Method.NewmarkMethod import NewmarkMethod
import matplotlib.pyplot as plt
import base
import numpy as np
a = NewmarkMethod(inp1,Iplt = 1)

s = inp1.nodes[-1]
u = []
w = []
t = []


for i in range(inp1.numTimeStep2):
    t.append(i*0.5/(inp1.numTimeStep2))
    u.append((s.R[1+i][0,0] - s.R[-1][0,0]) * 1000)
    w.append((s.R[1+i][2,0] - s.R[-1][2,0]) * 1000)

'''
for i in range(inp1.numTimeStep2):
    t.append(i*0.5/(inp1.numTimeStep2))
    R = np.matmul(base.theta2Q(s.iniAxiVector),s.lambd[1+i])
    u2 = R[1,2] * 0.0025
    u.append((s.R[1+i][1,0] - s.R[-1][1,0] + u2) * 1000)
    w.append((s.R[1+i][1,0] - s.R[-1][1,0]) * 1000)
'''


'''
for i in range(inp1.numTimeStep1):
    t.append(i*10.0/(inp1.numTimeStep1-1))
    R = np.matmul(base.theta2Q(s.iniAxiVector),s.lambd[1+i])
    u2 = R[1,2] * 0.0025
    u.append((s.R[1+i][1,0] - s.R[1][1,0] + u2) * 1000)
    w.append(0)
'''
  


'''
for i in range(inp1.numTimeStep1):
    t.append(i*10.0/(inp1.numTimeStep1-1))
    u.append((s.R[1+i][0,0] - s.R[1][0,0]) * 1000)
    w.append((s.R[1+i][2,0] - s.R[1][2,0]) * 1000)

'''

plt.plot(t,u)
plt.plot(t,w) 
plt.show()


