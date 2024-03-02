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

Name = '''SnapThrough1'''

'''calculation settings'''
numTimeStep1 = 10
numTimeStep2 = 1000
numTimeStep = numTimeStep1 + numTimeStep2
dt = []
T = []
for i in range(numTimeStep1):
    dt.append(100000.0)
    T.append(100000.0 * (i+1))
dt.append(0.069)
T.append(T[-1]+dt[-1])
T0 = T[-1]
for i in range(numTimeStep2-1):
    dt.append(0.001)
    T.append(T0 + 0.001 * (i+1))
    
'''material setting'''

materials = [] 
materials.append(IsotropicLE(4200000,0.3,1733.5)) 

'''cross-section setting'''
#joist steel GB 706-88 10
L0 = 0.120486
#L0 = 0.120925
L = 0.12

k = 0.004 * m.pi / L

H = 0.0025
W = 0.020
'''
L1 = (3 * L**3 * (1 + 0.25*k**2)*0.75) / (3*0.75*L**2-m.pi**3* H**2)
me = 3 * 0.004**2 *2000 / 8
e = (H**2 * k**2 / (3 * 0.75 * 0.004**2) + L1 * (1 + 0.75 * k**2)/L0 - 1) * 4200000 * k**2 / 2
d = m.log(2.07)
c = 2 * d * (me * e)**0.5 / (4*m.pi**2+d**2)**0.5
'''


A = H * W
I2 = W * H**3 / 12
I1 = H * W**3 / 12
crossSections = []
crossSections.append(Rod3D(A,I1,I2,np.array([[5.0/6.0],[5.0/6.0]])))

'''applied setting'''
AMList = []
gList = []
fList = []
for i in range(numTimeStep1):
    AMList.append([np.array([[0.0],[0.0],[0.0]]),np.zeros([3,3])])
    gList.append(np.array([[0],[0],[0]]))
    fList.append(np.zeros([3,1]))
for i in range(numTimeStep2):
    AMList.append([np.array([[0.0],[0.0],[0.0]]),np.zeros([3,3])])
    gList.append(np.array([[0],[0],[0]]))
    fList.append(np.zeros([3,1]))

'''node and element setting'''
numElement = 11
nodes = []                                                                   
elements = []
for i in range(numElement):

    nodes.append(HM3DCurvedRod(list(range(6*i,6*i+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0.0],[0.0],[-0.0723]]),
                               [np.array([[0.0],[0.0],[L0 * i / (2*numElement)]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               [np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]]),
                                np.array([[0.0],[0.0],[0.0]])],
                               AMList,
                               gList,
                               fList))

nodes.append(HM3DCurvedRod(list(range(6*numElement,6*numElement+6)),
                           materials[0],
                           crossSections[0],
                           np.array([[0.0],[0.0],[0.0]]),
                          [np.array([[0.0],[0.0],[L0 / 2.0]]),
                           np.array([[0.0],[0.0],[0.0]]),
                           np.array([[0.0],[0.0],[0.0]])],
                          [np.array([[0.0],[0.0],[0.0]]),
                           np.array([[0.0],[0.0],[0.0]]),
                           np.array([[0.0],[0.0],[0.0]])],
                          AMList,
                          gList,
                          fList))

for i in range(numElement):
     
    nodes.append(HM3DCurvedRod(list(range(6*(numElement+1+i),6*(numElement+1+i)+6)),
                               materials[0],
                               crossSections[0],
                               np.array([[0.0],[0.0],[0.0723]]),
                               [np.array([[0.0],[0.0],[L0 * (i+1) / (2*numElement) + L0 / 2]]),
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

DLoading1 = [] 
DLoading2 = []
DLoading3 = []

FLoading1 = [] 
FLoading2 = [] 
X = L - L0

for i in range(numTimeStep1):
    FLoading1.append(0.001)
    FLoading2.append('pass')
    DLoading1.append(0.0)
    DLoading2.append(X*(i+1)/numTimeStep1)
    DLoading3.append('pass')
    
FLoading1[-1] = 0.0
FLoading1[-2] = 0.0
FLoading1[-3] = 0.0
FLoading2.append('pass')
u = 0.004001727444909364 -0.00032

DLoading3.append(u)

   
for i in range(numTimeStep2):
    FLoading1.append(0.0)
    FLoading2.append(np.array([0.011,0.0,0.0]))
    DLoading1.append(0.0)
    DLoading2.append(X)
    DLoading3.append('pass')
    
    
    
    
for i in range(6):
    Dconstraints.append(Constraint([0,i],
                                   'displacement boundary',
                                   DLoading1))
    Dconstraints.append(Constraint([-1,i],
                                   'displacement boundary',
                                   DLoading1))
Dconstraints[5].loadList = DLoading2
Dconstraints.append(Constraint([numElement,0],
                                   'displacement boundary',
                                   DLoading3))
Fconstraints.append(Constraint([numElement-1,0],
                                   'force boundary',
                                   FLoading1))
for i in range(len(nodes)):
    Fconstraints.append(Constraint([i,0],
                                   'viscous boundary',
                                   FLoading2,)) 

'''
Fconstraints.append(Constraint(nodes[-1].dofID[1],
                                   'viscous boundary',
                                   FLoading,dof2=[-1,1])) 
Fconstraints.append(Constraint(nodes[-1].dofID[2],
                                   'viscous boundary',
                                   FLoading,dof2=[-1,2]))  
'''
                        

'''loading matrix and stiffness matrix'''    
numFreedom = nodes[-1].dofID[-1] + 1
P = np.zeros([numFreedom,1]) 
K = np.zeros([numFreedom,numFreedom])

'''coefficient matrix'''
charLength = (crossSections[0].J / crossSections[0].A)**0.5
coefficient = [1.0,1.0,1.0,charLength,charLength,charLength] * len(nodes)
coefficient = np.array([coefficient]).T 







