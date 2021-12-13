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


Name = 'NeoHookeanRoboticArm'

'''Calculation settings'''

Ttotal = 1                                                                     # total time
numIncrement = 10                                                             # Loading times 
Loading = np.ones(numIncrement)
TLoading = np.array(range(numIncrement)) / float(numIncrement-1)  


dt = Ttotal / numIncrement                                                     # time element
T = [0.0]                                                                      # time list

'''Geometric features'''  
L = 90                                                                        # the length of the Robotic arm (unit: mm)   
geometrys = []                                                                 # Geometry list
geometrys.append(geometry.g3DBeam(407.96,13649.73,26092.79,1))                            # (unit of area : mm^2) (unit of Moment of inertia : mm^4)

'''material setting'''
materials = []                                                                 # material list
materials.append(material.mLEIsotropy(0.1,0.5))                                 # (unit of modulus : N*mm^-2)                                                      

'''node and element setting''' 
knotVectors = [0,0,0,1,2,3,4,5,5,5]                                            # the knot vectors
[numElement,numNode,nodeList,C] = base.NUBSextraction(knotVectors,2)

      
nodes = []                                                                     # node list
elements = []                                                                  # element list

for i in range(numNode):
    nodes.append(node.n3DRoboticArm(np.array([[0.0],[0.0],[i * L / (numNode-1)]]),
                                        list(range(6*i,6*i+6)),
                                        np.array([[0],[0],[0.0]]),
                                        [[10.0],[np.array([9.0,0.0])]],
                                        np.array([[0],[0],[ 1.0*0.004508],[0],[0],[0]])
                                        )
                 )

for i in range(numElement):
    elements.append(element.kNH3DRoboticArmL3([nodes[nodeList[i][0]],nodes[nodeList[i][1]],nodes[nodeList[i][2]]],
                                                  materials[0],
                                                  geometrys[0],
                                                  Loading,
                                                  [TLoading],
                                                  1,
                                                  C[i]
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


Oconstraints.append(constraint.constraint(6 * numNode - 4,
                                              'elastic boundary',
                                              [0.0,2.0],
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
    


