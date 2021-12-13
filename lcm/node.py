#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 10:58:55 2021

@author: lixin
"""

import numpy as np
import base 

class n2DSolid():
    '''2D solid node'''
    dof = 2                                 #Degrees of freedom of the node


    def __init__(self, xyz, dofId, f=0.0):
        self.f = f                          #Body force
        self.historyXYZ = [xyz]             #Historical location of the node
        self.dofId = dofId                  #Degrees of freedom number
        
        self.U = np.array([[0.0],[0.0]])    #Node displacement
        self.dU = np.array([[0.0],[0.0]])   #Node displacement increment
     
                                                
    def updateU(self,incrementU):
        self.dU = np.array([[incrementU[self.dofId[0],0]],[incrementU[self.dofId[1],0]]])                
        self.U += self.dU
        
    def updateXYZ(self):
        self.historyXYZ.append(self.historyXYZ[0]+self.U)
        
        
class n3DSolid():
    '''3D solid node'''
    dof = 3                                 #Degrees of freedom of the node
    
    
    def __init__(self, xyz, dofId, f=0.0):
        self.f = f                          #Body force
        self.historyXYZ = [xyz]             #Historical location of the node
        self.dofId = dofId                  #Degrees of freedom number
        
        self.U = np.array([[0.0],[0.0],[0.0]])
        self.dU = np.array([[0.0],[0.0],[0.0]])
        

    def updateU(self,incrementU):
        self.dU = np.array([[incrementU[self.dofId[0],0]],[incrementU[self.dofId[1],0]],[incrementU[self.dofId[2],0]]])                
        self.U += self.dU
        
    def updateXYZ(self):
        self.historyXYZ.append(self.historyXYZ[0]+self.U)

class n3DRoboticArm():
    '''3D robotic arm node'''
    dof = 6                                 #Degrees of freedom of the node
    
    
    def __init__(self, xyz, dofId, theta, cable, f): 
        self.f = f                          #Body force
        self.historyXYZ = [xyz]             #Historical location of the node
        self.dofId = dofId                  #Degrees of freedom number
        
        self.U = np.array([[0.0],[0.0],[0.0]])
        self.dU = np.array([[0.0],[0.0],[0.0]])
        
        self.theta = np.array(theta)        #Axial vector

        
        self.dTheta = np.array([[0.0],[0.0],[0.0]])       
        
        self.cableT = cable[0]             #Magnitude of the tendon tension
        self.cableXia = cable[1]           #Coordinates of cable off center(use the list)
        

            
        
    def updateU(self,incrementU):
        self.dU = np.array([[incrementU[self.dofId[0],0]],[incrementU[self.dofId[1],0]],[incrementU[self.dofId[2],0]]])                
        self.U += self.dU
        
        self.dTheta = np.array([[incrementU[self.dofId[3],0]],[incrementU[self.dofId[4],0]],[incrementU[self.dofId[5],0]]])

    def updateXYZ(self):
        self.historyXYZ.append(self.historyXYZ[0]+self.U)
        
        
        
        
     
        
     
        
     
        
     
        
     
        
     
        
     
        
     
        
     
        
     
        
     
        