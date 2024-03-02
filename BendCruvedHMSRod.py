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

Name = '''BendCruvedHMSRod'''

'''calculation settings'''
numTimeStep1 = 101
numTimeStep2 = 0
numTimeStep = numTimeStep1 + numTimeStep2
dt = []
T = []
for i in range(numTimeStep1):
    dt.append(100000.0)
    T.append(100000.0 * (i+1))
T0 = T[-1]
for i in range(numTimeStep2):
    dt.append(1.0)
    T.append(T0 + 1.0 * (i+1))
    
'''material setting'''

materials = [] 

materials.append(IsotropicLE(171400,0.499,1741.23)) 
'''
materials.append(IsotropicLE(190000,0.499,1871.72)) 
'''
'''cross-section setting'''
R = 0.015
a = 0.005


A = a**2
I2 = a**4 / 12
I1 = a**4 / 12
crossSections = []
crossSections.append(Rod3D(A,I1,I2,np.array([[5.0/6.0],[5.0/6.0]])))

'''applied setting'''
AMList = []
gList = []
fList = []
for i in range(numTimeStep1):
    AMList.append([np.array([[i*0.01/(numTimeStep1-1)],[0.0],[0.0]]),np.zeros([3,3])])
    gList.append(np.array([[0],[0],[9.8]]))
    fList.append(np.zeros([6,1]))
for i in range(numTimeStep2):
    AMList.append([np.array([[0.0],[0.0],[0.0]]),np.zeros([3,3])])
    gList.append(np.array([[0],[0],[9.8]]))
    fList.append(np.zeros([6,1]))

'''node and element setting'''
numElement = 10
nodes = []                                                                   
elements = []
for i in range(2*numElement+1):
    theta = m.pi * 0.5 * i / (2*numElement)
    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0.0],[0.0],[0.07710]]),
                               [np.array([[R - R*m.cos(theta)],[0.0],[R*m.sin(theta)]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[theta],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList,
                               KI = 0.5))



for i in range(numElement):
    elements.append(HM3DCurvedRodL3([nodes[2*i],nodes[2*i+1],nodes[2*i+2]],1))
                  
   
'''boundary condition setting'''
Dconstraints = []                                                              # Displacement boundary condition list
Fconstraints = []                                                              # Force boundary condition list

DLoading1 = [] 
FLoading1 = [] 

for i in range(numTimeStep1):
    FLoading1.append([0.0,0.0,0.0])
    DLoading1.append(0.0)

   
for i in range(numTimeStep2):
    FLoading1.append(np.array([0.000000,0.0,0.0]))
    DLoading1.append(0.0)
    
    
    
    
    
for i in range(6):
    Dconstraints.append(Constraint([0,i],
                                   'displacement boundary',
                                   DLoading1))

for i in range(len(nodes)):
    Fconstraints.append(Constraint([i,0],
                                   'viscous boundary',
                                   FLoading1,)) 
    Fconstraints.append(Constraint([i,1],
                                   'viscous boundary',
                                   FLoading1,)) 
    Fconstraints.append(Constraint([i,2],
                                   'viscous boundary',
                                   FLoading1,)) 

                        

'''loading matrix and stiffness matrix'''    
numFreedom = nodes[-1].dofID[-1] + 1
P = np.zeros([numFreedom,1]) 
K = np.zeros([numFreedom,numFreedom])

'''coefficient matrix'''
charLength = (crossSections[0].J / crossSections[0].A)**0.5
coefficient = [1.0,1.0,1.0,charLength,charLength,charLength] * len(nodes)
coefficient = np.array([coefficient]).T 






