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


Name = 'RoboticArm'

'''Calculation settings'''

Ttotal = 1                                                                     # total time
numIncrement = 4                                                               # Loading times  
Loading = np.array([0.1,0.2,0.5,1.0])                                      # Loading path  
TLoading = np.array(Loading) 

dt = Ttotal / numIncrement                                                     # time element
T = [0.0]                                                                      # time list

'''Geometric features'''  
L = 1.0                                                                        # the length of the Robotic arm (unit: mm)   
geometrys = []                                                                 # Geometry list
geometrys.append(geometry.g3DBeam(1,0.005,0.005,1))                            # (unit of area : mm^2) (unit of Moment of inertia : mm^4)

'''material setting'''
materials = []                                                                 # material list
materials.append(material.mLEIsotropy(1,-0.5))                                 # (unit of modulus : N*mm^-2)                                                      

'''node and element setting'''       
numElement = 3                                                                 # the number of the elements
nodes = []                                                                     # node list
elements = []                                                                  # element list

for i in range(2*numElement+1):
    nodes.append(node.n3DRoboticArm(np.array([[0.0],[0.0],[i * L / (2 * numElement)]]),
                                        list(range(6*i,6*i+6)),
                                        np.array([[0],[0],[0]]),
                                        [[0.05*m.pi],[np.array([0.1,0.0])]],
                                        np.array([[0],[0],[0],[0],[0],[0]])
                                        )
                 )

for i in range(numElement):
    elements.append(element.eLE3DRoboticArmL3([nodes[2*i],nodes[2*i+1],nodes[2*i+2]],
                                                  materials[0],
                                                  geometrys[0],
                                                  Loading,
                                                  [TLoading],
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


Oconstraints.append(constraint.constraint(6 * (2*numElement + 1) - 5,
                                              'elastic boundary',
                                              [0.0,-0.020],
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
    



















    