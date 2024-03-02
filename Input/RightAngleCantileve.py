# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 20:28:44 2023

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


Name = '''RightAngleCantileve'''


'''calculation settings'''
numTimeStep = 120
dt = []
T = []
for i in range(numTimeStep):
    dt.append(30.0 / numTimeStep)
    T.append(30.0 / numTimeStep * (i+1))
    
'''material setting'''
materials = [] 
materials.append(IsotropicLE(1e6,0.0,1.0)) 

'''cross-section setting'''
crossSections = []
crossSections.append(Rod3D(1.0,1e-3,1e-3,np.array([[2.0],[2.0]])))

'''applied setting'''
AMList = []
gList = []
fList = []
for i in range(numTimeStep):
    AMList.append([np.zeros([3,1]),np.zeros([3,3])])
    gList.append(np.zeros([3,1]))
    fList.append(np.zeros([3,1]))

'''node and element setting'''
numElement = 5
nodes = []                                                                   
elements = []

for i in range(numElement * 2):
    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0],[0],[0]]),
                               [np.array([[0.0],[0.0],[i*5.0/numElement]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList))

nodes.append(HM3DCurvedRod(list(range(6*2*numElement,6*2*numElement+6)),
                           materials[0],
                           crossSections[0],
                           np.array([[0],[0],[0]]),
                           [np.array([[0.0],[0.0],[10.0]]),
                            np.array([[0.0],[0.0],[0.0]]),
                            np.array([[0.0],[0.0],[0.0]])],
                           [np.array([[0.0],[m.pi*0.0],[0.0]]),
                            np.array([[0.0],[0.0],[0.0]]),
                            np.array([[0.0],[0.0],[0.0]])],
                           AMList,
                           gList,
                           fList))

for i in range(numElement*2):
    nodes.append(HM3DCurvedRod(list(range(6 *(2*numElement+1+i) ,6*(2*numElement+2+i))),
                               materials[0],
                               crossSections[0],
                               np.array([[0],[0],[0]]),
                               [np.array([[(i+1)*5.0/numElement],[0.0],[10.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[m.pi*0.5],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList))

for i in range(numElement):
    elements.append(HM3DCurvedRodL3([nodes[2*i],nodes[2*i+1],nodes[2*i+2]],1))
nodes[numElement * 2].iniAxiVector = np.array([[0.0],[m.pi*0.5],[0.0]])
for i in range(numElement):
    elements.append(HM3DCurvedRodL3([nodes[2*(i+numElement)],nodes[2*(i+numElement)+1],nodes[2*(i+numElement)+2]],1))


D = np.diag([1e6,1e6,1e6,1e3,1e3,1e3])
R = np.zeros([6,6])
R[0:3,0:3] = R[3:6,3:6] = elements[-1].guassPoints[0].Q0
for i in elements[numElement].guassPoints:
    i.D = np.matmul(R,np.matmul(D,R.T))
    i.Q0 = R[0:3,0:3]
Jrho = np.diag([10,10,20])
for i in elements:
    for j in i.guassPoints:
        j.Jrho = np.matmul(j.Q0,np.matmul(Jrho,j.Q0.T))


   
'''boundary condition setting'''
Dconstraints = []                                                              # Displacement boundary condition list
Fconstraints = []                                                              # Force boundary condition list

FLoading = []  
DLoading = []
ii = int(1.0 / dt[0])
for i in range(ii):
    FLoading.append((i + 1) * dt[0] * 50.0)
for i in range(ii):
    FLoading.append(50.0 - (i + 1) * dt[0] * 50.0)
for i in range(numTimeStep):
    DLoading.append(0.0)
    FLoading.append(0.0)
    
for i in range(6):
    Dconstraints.append(Constraint(i,
                                   'displacement boundary',
                                   DLoading))
Fconstraints.append(Constraint(6 * 2 * numElement + 1,
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
    