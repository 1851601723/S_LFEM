# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 15:59:09 2023

@author: XinLi
"""

import numpy as np
import math as m
from os import getcwd
import sys
path = getcwd()
path = list(path)
path = path + ['\\','B','a','s','e']
path = ''.join(path)
sys.path.append(path)

from CrossSection.Rod3D import Rod3D
from Constraint.Constraint import Constraint 
from Element.HM3DClassicalRodL2 import HM3DClassicalRodL2
from Material.IsotropicLE import IsotropicLE
from Node.HM3DCurvedRod import HM3DCurvedRod

Name = '''ClassicalCircularRingBend'''


'''calculation settings'''
numTimeStep = 50
dt = []
T = []
for i in range(numTimeStep):
    dt.append(1000.0 / numTimeStep)
    T.append(1000.0 / numTimeStep * (i+1))
    
'''material setting'''
# steel
materials = [] 
materials.append(IsotropicLE(1e6,0.0,1.0)) 

'''cross-section setting'''
#joist steel GB 706-88 10
h = 0.1
d = 0.005
b = 0.07
tc = 0.01
'''
A = 2 * tc * b + d * (h - 2 * tc)
I1 = (2 * b**3 * tc + d**3 * (h - 2*tc)) / 12
I2 = d * (h - 2*tc)**3 / 12 + b * tc ** 3 / 6 + 2 * b * tc * (h / 2 - tc / 2) ** 2
'''
A = h * b
I1 = h * b**3 / 12
I2 = b * h**3 / 12
crossSections = []
crossSections.append(Rod3D(A,I1,I2,np.array([[5.0/6.0],[5.0/6.0]])))

'''applied setting'''
AMList = []
gList = []
fList = []
for i in range(numTimeStep):
    AMList.append([np.zeros([3,1]),np.zeros([3,3])])
    gList.append(np.zeros([3,1]))
    fList.append(np.zeros([3,1]))

'''node and element setting'''
numElement = 12
nodes = []                                                                   
elements = []
r = 0.1
for i in range(numElement+1):
    theta = i * m.pi / numElement
    
    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0],[0],[0]]),
                               [np.array([[-r * m.cos(theta)],[0.0],[r * m.sin(theta)]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[theta],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList))


for i in range(numElement):
    elements.append(HM3DClassicalRodL2([nodes[i],nodes[i+1]],1))
                  
   
'''boundary condition setting'''
Dconstraints = []                                                              # Displacement boundary condition list
Fconstraints = []                                                              # Force boundary condition list

DLoading = []  
FLoading = []
M = materials[0].E * I2 / (m.pi * r)
F = -150.0
for i in range(numTimeStep):
    DLoading.append(0.0)
    FLoading.append(F * (i + 1) / numTimeStep)
    
for i in range(6):
    Dconstraints.append(Constraint(i,
                                   'displacement boundary',
                                   DLoading))

Fconstraints.append(Constraint(nodes[-1].dofID[0],
                               'force boundary',
                               FLoading))
                        

'''loading matrix and stiffness matrix'''    
numFreedom = nodes[-1].dofID[-1] + 1
P = np.zeros([numFreedom,1]) 
K = np.zeros([numFreedom,numFreedom])

'''coefficient matrix'''
charLength = (crossSections[0].J / crossSections[0].A)**0.5
coefficient = [1.0,1.0,1.0,charLength,charLength,charLength] * len(nodes)
coefficient = np.array([coefficient]).T 
