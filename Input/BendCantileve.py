# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 20:33:08 2023

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
from Element.HM3DCurvedRodL2 import HM3DCurvedRodL2
from Material.IsotropicLE import IsotropicLE
from Node.HM3DCurvedRod import HM3DCurvedRod


'''BendCantileve'''
'''calculation settings'''
numTimeStep = 8
dt = []
T = []
for i in range(numTimeStep):
    dt.append(100000.0)
    T.append(100000.0 * (i+1))
    
'''material setting'''
materials = [] 
materials.append(IsotropicLE(1e7,0.0,1e4)) 

'''cross-section setting'''
crossSections = []
crossSections.append(Rod3D(1.0,1.0/12.0,1.0/12.0,np.array([[1.0],[1.0]])))

'''applied setting'''
AMList = []
gList = []
fList = []
for i in range(numTimeStep):
    AMList.append([np.zeros([3,1]),np.zeros([3,3])])
    gList.append(np.zeros([3,1]))
    fList.append(np.zeros([3,1]))

'''node and element setting'''
numElement = 10
nodes = []                                                                   
elements = []
R = 100

for i in range(numElement+1):
    theta = i * m.pi / (4 * numElement)
    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0],[0],[0]]),
                               [np.array([[R - R * m.cos(theta)],[0.0],[R * m.sin(theta)]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[theta],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList))


for i in range(numElement):
    elements.append(HM3DCurvedRodL2([nodes[i],nodes[i+1]],1))

'''boundary condition setting'''
Dconstraints = []                                                              # Displacement boundary condition list
Fconstraints = []   
FLoading = [300,450,600,1000,1500,1000,600,450]  
DLoading = []

for i in range(numTimeStep):
    DLoading.append(0.0)
    
    
for i in range(6):
    Dconstraints.append(Constraint(i,
                                   'displacement boundary',
                                   DLoading))
Fconstraints.append(Constraint(6 * numElement + 1,
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
    