# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 17:44:39 2023

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

Name = '''BucklingBend'''


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
materials.append(IsotropicLE(1.0,0.0,0.0001)) 

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
L = 10.0
for i in range(numElement+1):
    
    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0],[0],[0]]),
                               [np.array([[0.0],[0.0],[L * i / numElement]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[0.0],[0.0]]),
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
FLoading1 = []
FLoading2 = []
M = -2 * materials[0].E * crossSections[0].I1 * m.pi**2 / (L**2)
F = 1e-3 * M
for i in range(numTimeStep):
    DLoading.append(0.0)
    FLoading1.append(M * i / numTimeStep)
    FLoading2.append(F)
    
    
for i in range(3):
    Dconstraints.append(Constraint(i,
                                   'displacement boundary',
                                   DLoading))
Dconstraints.append(Constraint(nodes[-1].dofID[0],
                               'displacement boundary',
                                DLoading))
Dconstraints.append(Constraint(nodes[-1].dofID[1],
                               'displacement boundary',
                                DLoading))
ii = int(numElement / 2)
Fconstraints.append(Constraint(nodes[-1].dofID[2],
                               'force boundary',
                               FLoading1))
Fconstraints.append(Constraint(nodes[ii].dofID[0],
                               'force boundary',
                               FLoading2))
                        

'''loading matrix and stiffness matrix'''    
numFreedom = nodes[-1].dofID[-1] + 1
P = np.zeros([numFreedom,1]) 
K = np.zeros([numFreedom,numFreedom])

'''coefficient matrix'''
charLength = (crossSections[0].J / crossSections[0].A)**0.5
coefficient = [1.0,1.0,1.0,charLength,charLength,charLength] * len(nodes)
coefficient = np.array([coefficient]).T 

