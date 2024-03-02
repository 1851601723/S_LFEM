# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 16:54:12 2023

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


'''PureBend'''


'''calculation settings'''
numTimeStep = 8
dt = []
T = []
for i in range(numTimeStep):
    dt.append(100000.0)
    T.append(100000.0 * (i+1))
    
'''material setting'''
materials = [] 
materials.append(IsotropicLE(1.0,0.0,1e3)) 

'''cross-section setting'''
crossSections = []
crossSections.append(Rod3D(10.0,1.0,1.0,np.array([[2.0],[2.0]])))

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

for i in range(numElement+1):
    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0],[0],[0]]),
                               [np.array([[0.0],[0.0],[1.0 * i / numElement]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList))


for i in range(numElement):
    elements.append(HM3DCurvedRodL2([nodes[i],nodes[i+1]],1))
'''
Jrho = np.diag([10,10,20])
for i in elements:
    for j in i.guassPoints:
        j.Jrho = np.matmul(j.Q0,np.matmul(Jrho,j.Q0.T))
'''
'''boundary condition setting'''
Dconstraints = []                                                              # Displacement boundary condition list
Fconstraints = []                                                              # Force boundary condition list

FLoading = [1.0*m.pi,2.0*m.pi,3.0*m.pi,4.0*m.pi,3.0*m.pi,2.0*m.pi,1.0*m.pi]  
DLoading = []

for i in range(numTimeStep):
    DLoading.append(0.0)
    FLoading.append(0.0)
    
for i in range(6):
    Dconstraints.append(Constraint(i,
                                   'displacement boundary',
                                   DLoading))
Fconstraints.append(Constraint(6 * numElement + 4,
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
    
