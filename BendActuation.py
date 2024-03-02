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
from Element.HM3DCurvedRodL3 import HM3DCurvedRodL3
from Material.IsotropicLE import IsotropicLE
from Node.HM3DCurvedRod import HM3DCurvedRod

Name = '''BendActuation'''

'''calculation settings'''
numTimeStep = 28
dt = []
T = []
for i in range(numTimeStep):
    dt.append(100000.0)
    T.append(100000.0 * (i+1))
    
'''material setting'''
# steel
materials = [] 
materials.append(IsotropicLE(114210,0.499,1763.27)) 

'''cross-section setting'''
#joist steel GB 706-88 10
L = 0.03
D = 0.007 

A = D ** 2
I1 = D**4 / 12
I2 = D**4 / 12
crossSections = []
crossSections.append(Rod3D(A,I1,I2,np.array([[5.0/6.0],[5.0/6.0]])))

'''applied setting'''
AMList = []
gList = []
fList = []
for i in range(numTimeStep):
    AMList.append([np.array([[0.014 * (i+1) / numTimeStep],[0.0],[0.0]]),np.zeros([3,3])])
    gList.append(np.array([[0],[0],[9.8]]))
    fList.append(np.zeros([3,1]))

'''node and element setting'''
numElement = 5
nodes = []                                                                   
elements = []
for i in range(2*numElement+1):
    
    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0],[0],[0.0785]]),
                               [np.array([[0.0],[0.0],[i * L / (2 * numElement)]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList))


for i in range(numElement):
    elements.append(HM3DCurvedRodL3([nodes[2*i],nodes[2*i+1],nodes[2*i+2]],1))
                  
   
'''boundary condition setting'''
Dconstraints = []                                                              # Displacement boundary condition list
Fconstraints = []                                                              # Force boundary condition list

DLoading = []  

for i in range(numTimeStep):
    DLoading.append(0.0)
    
    
for i in range(6):
    Dconstraints.append(Constraint(i,
                                   'displacement boundary',
                                   DLoading))


                        

'''loading matrix and stiffness matrix'''    
numFreedom = nodes[-1].dofID[-1] + 1
P = np.zeros([numFreedom,1]) 
K = np.zeros([numFreedom,numFreedom])

'''coefficient matrix'''
charLength = (crossSections[0].J / crossSections[0].A)**0.5
coefficient = [1.0,1.0,1.0,charLength,charLength,charLength] * len(nodes)
coefficient = np.array([coefficient]).T 

