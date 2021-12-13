
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 15:21:08 2021

@author: XinLi
"""

import numpy as np
import math as m
import geometry
import material
import node
import element
import constraint
import base


Name = '45degreeBend'

'''Calculation settings'''

Ttotal = 1                                                                     # total time
numIncrement = 3                                                               # Loading times  
Loading = np.array(range(1,numIncrement+1)) / numIncrement                     # Loading path  
TLoading1 = np.array(Loading) 
TLoading2 = np.array(Loading)
TLoading3 = np.array(Loading)
dt = Ttotal / numIncrement                                                     # time element
T = [0.0]                                                                      # time list

'''Geometric features'''  
R = 100.0                                                                      # the radiu of the bend (unit: m)   
geometrys = []                                                                 # Geometry list
geometrys.append(geometry.g3DBeam(1,1/12,1/12,1))                              # (unit of area : m^2) (unit of Moment of inertia : m^4)

'''material setting'''
materials = []                                                                 # material list
materials.append(material.mLEIsotropy(10**7,0.0))                              # (unit of modulus : N*m^-2)                                                      

'''node and element setting'''       
numElement = 8                                                                 # the number of the elements
nodes = []                                                                     # node list
elements = []                                                                  # element list

for i in range(numElement+1):
    nodes.append(node.n3DRoboticArm(np.array([[R * m.cos(i * m.pi/(4 * numElement)) - R],[R * m.sin(i * m.pi/(4 * numElement))],[0]]),
                                        list(range(6*i,6*i+6)),
                                        base.Q2theta(np.array([[0,m.cos(i * m.pi/(4 * numElement)),-m.sin(i * m.pi/(4 * numElement))],
                                                               [0,m.sin(i * m.pi/(4 * numElement)), m.cos(i * m.pi/(4 * numElement))],
                                                               [1,0,0]])),
                                        [[0.00,0.00,0.00],[np.array([0.5,0.0]),np.array([ -0.25,0.25*3**0.5]),np.array([-0.25,-0.25*3**0.5])]],
                                        np.array([[0],[0],[0],[0],[0],[0]])
                                        )
                 )

for i in range(numElement):
    elements.append(element.eLE3DRoboticArmL2([nodes[i],nodes[i+1]],
                                                  materials[0],
                                                  geometrys[0],
                                                  Loading,
                                                  [TLoading1,TLoading2,TLoading3],
                                                  1
                                                  )
                    )


'''boundary condition setting'''
Dconstraints = []                                                              # Displacement boundary condition list
Oconstraints = []                                                              # Force boundary condition list

for i in range(6):
    Dconstraints.append(constraint.constraint(i,
                                              'displacement boundary',
                                              [0.0],
                                              Loading
                                              )
                        )


Oconstraints.append(constraint.constraint(6 * (numElement + 1) - 2,
                                              'elastic boundary',
                                              [0.0,0.0],
                                              Loading
                                              )
                    )
Oconstraints.append(constraint.constraint(6 * (numElement + 1) - 4,
                                              'elastic boundary',
                                              [0.0,600.0],
                                              Loading
                                              )
                    )    


'''loading matrix and stiffness matrix'''    
numFreedom = nodes[-1].dofId[-1] + 1
P = np.zeros([numFreedom,1]) 
K = np.zeros([numFreedom,numFreedom])


'''coefficient matrix'''
charLength = (geometrys[0].J / geometrys[0].A)**0.5
coefficient = [1.0,1.0,1.0,charLength,charLength,charLength] * len(nodes)
coefficient = np.array([coefficient]).T 
    


















