# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:22:55 2023

@author: XinLi
"""

from os import getcwd, path
import sys

paths = getcwd()
paths = path.dirname(paths)
paths = list(paths)
paths = paths + ['\\','B','a','s','e']
paths = ''.join(paths)
sys.path.append(paths)


import numpy as np
import base

class HM3DCurvedRod():
    
    def __init__(self,dofID,material,crossSection,RM,translation,rotation,AMList,gList,fList,KI = 1):
        
        '''the quantities never change in the calculation process'''
        self.KI = KI
        self.dof = 6                                # degrees of freedom
        self.dofID = dofID                          # degrees of freedom number
        self.material = material                    # material parameters
        self.crossSection = crossSection            # geometric parameters of the cross-section
        self.RM = RM                                # residual magnetic flux density
        self.iniAxiVector = rotation[0]             # initial axial vector
        self.R0 = base.theta2Q(rotation[0])
        
        '''the quantities change in the calculation process'''
        self.R = [translation[0]]                   # radius vector list
        self.dR = [translation[1]]                  # linear velocity list
        self.ddR = [translation[2]]                 # linear acceleration list
        
        
        self.axiVector = [np.zeros([3,1])]          # axial vector list
        self.lambd = [np.eye(3)]                    # effective rotation matrix list 
        self.W = [rotation[1]]                      # material angular velocity list
        self.dW = [rotation[2]]                     # material angular acceleration list
        self.w = [np.matmul(self.R0,self.W[0])]
        self.dw = [np.matmul(self.R0,self.dW[0])]
        
        self.U = np.zeros([3,1])                    # linear displacement vector in time step
        self.dU = np.zeros([3,1])                   # linear displacement incremental vector in iteration
        
        self.Theta = np.zeros([3,1])                # material angular displacement vector in time step
        self.theta = np.zeros([3,1])                # angular displacement vector in time step
        self.dtheta = np.zeros([3,1])               # angular displacement incremental vector in iteration
        
        self.AMList = AMList                        # applied magnetic flux density List
        self.gList = gList                          # body force List
        self.fList = fList                          # distuributed force List
    
    def timeStepUpdate(self,i,dt,beta,gamma):
        self.AM = self.AMList[i]
        self.g = self.gList[i]
        self.f = self.fList[i]
        
        
        self.U = np.zeros([3,1]) 
        self.R.append(self.R[-1] + self.U)
        ddR = -self.dR[-1] / (dt * beta) - (0.5 - beta) * self.ddR[-1] / beta
        dR = self.dR[-1] + dt * ((1.0-gamma) * self.ddR[-1] + gamma * ddR)
        self.dR.append(np.array(dR))
        self.ddR.append(np.array(ddR))
        
        self.Theta = np.zeros([3,1])
        dW = -self.W[-1] / (dt * beta) - (0.5 - beta) * self.dW[-1] / beta
        W = self.W[-1] + dt * ((1.0-gamma) * self.dW[-1] + gamma * dW)
        self.dW.append(np.array(dW))
        self.W.append(np.array(W))
        self.theta = np.matmul(self.lambd[-1],self.Theta)
        self.lambd.append(np.matmul(base.theta2Q(self.theta),self.lambd[-1]))
        self.axiVector.append(base.Q2theta(self.lambd[-1]))
        Rt = np.matmul(self.R0,self.lambd[-1])
        self.w.append(np.matmul(Rt,W))
        self.dw.append(np.matmul(Rt,dW))
        
    def iterationUpdate(self,DU,dt,beta,gamma):
        self.dU = np.array([[DU[self.dofID[0],0]],
                            [DU[self.dofID[1],0]],
                            [DU[self.dofID[2],0]]])  
        self.dtheta = np.array([[DU[self.dofID[3],0]],
                                [DU[self.dofID[4],0]],
                                [DU[self.dofID[5],0]]])
        
        self.U = self.U + self.dU
        self.R[-1] = self.R[-1] + self.dU
        self.dR[-1] = self.dR[-1] + self.dU * gamma / (dt * beta)
        self.ddR[-1] = self.ddR[-1] + self.dU  / (dt**2 * beta)
        
        self.lambd[-1] = np.matmul(base.theta2Q(self.dtheta),self.lambd[-1])
        self.axiVector[-1] = base.Q2theta(self.lambd[-1])
        self.theta = base.Q2theta(np.matmul(base.theta2Q(self.dtheta),base.theta2Q(self.theta)))
        Theta = np.matmul(self.lambd[-2].T,self.theta)
        DTheta = Theta - self.Theta
        self.W[-1] = self.W[-1] + DTheta * gamma / (dt * beta)
        self.dW[-1] = self.dW[-1] + DTheta / (dt**2 * beta)
        self.Theta = self.Theta + DTheta
        Rt = np.matmul(self.R0,self.lambd[-1])
        self.w[-1] = np.matmul(Rt,self.W[-1])
        self.dw[-1] = np.matmul(Rt,self.dW[-1])
        
    
 
                               
                