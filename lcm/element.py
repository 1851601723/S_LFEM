#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 21:10:13 2021

@author: lixin
"""

import math as m 
import numpy as np
import base 

class eLE3DRoboticArmL2():
    
    '''linear elasticity 3D robotic arm L2 element'''
    
    def __init__(self, nodes, material, geometry, FLoading, TLoading, Reduce):
        
        self.eps = 10 ** -7
        self.numNode = 2                                                       #Number of the nodes in the element
        self.numDof = 6                                                        #Number of the degrees of freedom of the node
        
        '''input the information of the element'''
        self.nodes = nodes                                                     #Nodes of the element
        self.dofId = []                                                        #Degrees of freedom of the element
        for i in self.nodes:
            self.dofId.extend(i.dofId)
            
        self.material = material                                               #Material of the element
        self.geometry = geometry                                               #Geometry of the element
        
        self.FLoading = FLoading                                               #Load path of the body force
        self.TLoading = TLoading                                               #Load path of the tendon tension
        
        if(Reduce == 1):
            self.numGuass = [1]                                                #reduced integration 
        elif(Reduce == 0):
            self.numGuass = [2]                                                #classical integration
        else:
            self.numGuass = [3]                                                #Highest order integration
            
        '''Define material stiffness matrix'''
        self.CN0 = np.diag([material.G * geometry.A / geometry.k, 
                            material.G * geometry.A / geometry.k,
                            material.E * geometry.A]) 
        
        self.CM0 = np.diag([material.E * geometry.Ix,
                            material.E * geometry.Iy,
                            material.G * geometry.J]) 
        
        '''Define guass integration points'''
        self.guassPoints = []
        for i in range(self.numGuass[0]):
            self.guassPoints.append(self.guassPoint(self, 
                                               base.guassLocation[self.numGuass[0]-1][i],
                                               base.guassWeight[self.numGuass[0]-1][i]))
        '''Body force matrix'''
        self.BF = np.zeros([self.numNode * self.numDof, 1])
        for i in range(self.numNode):
            for j in range(len(self.guassPoints)):
                for k in range(self.numDof):
                    self.BF[i*self.numDof+k,0] += (self.guassPoints[j].guassWeight * self.guassPoints[j].f[k] * 
                                                self.guassPoints[j].N[i,0] * self.guassPoints[j].detJ)
                    
    
    '''guass intergration point class'''
    class guassPoint():
            
        def __init__(self, element, guassLocation, guassWeight):
                
            self.guassLocation = guassLocation                                 #Location of the guass intergation point                            
            self.guassWeight = guassWeight                                     #Weight of the guass intergation point
                
            '''shape function and its derivatives'''
            self.N = np.array([[(1 - self.guassLocation) / 2],                 #node 1 (-1)
                               [(1 + self.guassLocation) / 2]])                #node 2 ( 1)     
                
            self.dN = np.array([[-0.5],                                    
                                [ 0.5]])
                
            '''Jacobian matrix and its determinant'''
            self.J = np.zeros([1,1])
            dL = np.zeros([3,1])                                           
            for i in range(1):
                for j in range(1):
                    for k in range(element.numNode):
                        dL += self.dN[k,j] * element.nodes[k].historyXYZ[0]
            self.J[0,0] += float(sum(dL ** 2)) ** 0.5
            self.detJ = np.linalg.det(self.J)
            assert self.detJ>=0, 'Negative volume in element'
                
            self.dN_ds = np.matmul(self.dN,np.linalg.inv(self.J))              #Derivatives of the shape function with respect to the arc length
                
                
            self.f = np.zeros([element.numDof,1])                              #Body force in the guass intergation point
            for i in range(element.numNode):
                self.f += self.N[i,0] * element.nodes[i].f
                    

                
            
                
            '''variables with respect to configuration'''
                
            self.lambd = np.eye(3)                                             #configuration matrix
                
            self.dTheta = np.zeros([3,1])                                      #axial vector increment
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].theta
                
            self.Q = np.eye(3)                                                 #roration matrix
                
            self.dr_ds = np.zeros([3,1])                                       #Derivative of vector diameter with respect to arc length
            self.dTheta_ds = np.zeros([3,1])                                   #Derivative of axial vector increment with respect to arc length
            for i in range(element.numNode):
                self.dr_ds += self.dN_ds[i,0] * element.nodes[i].historyXYZ[0]
                self.dTheta_ds += self.dN_ds[i,0] * element.nodes[i].theta     #float(dN_ds[i][0])
                
            '''variables with respect to strain'''
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5               
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                self.Q = np.matmul(dlambd,self.Q) 
                omega0 = (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                          + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                          + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
            else:
                omega0 = self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T                #strain with respect to the displacement in current configuration 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                    #strain with respect to the displacement in reference configuration 
            self.omega = np.zeros([3,1])                                       #strain with respect to the roration in current configuration
            self.kappa = np.matmul(self.lambd.T,self.omega)                    #strain with respect to the roration in reference configuration                      
            self.Omega0 = base.V2skew(omega0)                                  #Initial curvature
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
            self.E1 = np.array([self.Q[:,0]]).T                                #Vector E1
            self.E2 = np.array([self.Q[:,1]]).T                                #Vector E2
            self.E = [self.E1,self.E2]           
                 
                    
            '''variables with respect to cable should be think about'''
            self.cableXia = []                                                 #Coordinates of cable off center
            self.dXia_ds = []                                                  #Derivatives of the Coordinates of the cable off center with respect to the arc length
            self.cableT = []                                                   #Magnitude of the tendon tension    
            for i in range(len(element.nodes[0].cableXia)):
                self.cableXia.append(np.zeros(2))
                self.dXia_ds.append(np.zeros(2))
                self.cableT.append(0)
                for j in range(element.numNode):
                    self.cableXia[i] += self.N[j,0] * element.nodes[j].cableXia[i]
                    self.dXia_ds[i] += self.dN_ds[j,0] * element.nodes[j].cableXia[i]
                    self.cableT[i] += self.N[j,0] * element.nodes[j].cableT[i]
                        
            self.cableXi0 = []
                
            for i in self.cableXia:
                self.cableXi0.append(i[0] * self.E1 + i[1] * self.E2)          #Vector of cable off center before deformation
            self.cableXi = list(self.cableXi0)   
                          #Vector of cable off center after deformation
            self.cable_dr_ds = []
            self.cable_dr_ds_bar = []
            for i in range(len(self.cableXia)):
                self.cable_dr_ds.append(np.array(self.dr_ds))
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar.append(np.zeros([3,1]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                
              
              
            '''variables with respect to interal force'''
            self.CN = np.matmul(self.Q,np.matmul(element.CN0,self.Q.T))        #the material elastic tensor in material configuration
            self.CM = np.matmul(self.Q,np.matmul(element.CM0,self.Q.T))        #the material elastic tensor in material configuration
                
            self.cableF = []                                                   #the F caused by cable
            self.cableM = []                                                   #the M caused by cable
            for i in range(len(self.cable_dr_ds_bar)):
                self.cableF.append(np.zeros([3,1]))
                self.cableM.append(np.zeros([3,1]))
                    
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
                
            self.NuB = np.matmul(self.CN,self.Gamma)                        
            self.MuB = np.matmul(self.CM,self.kappa)
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)                

                
            self.nu = np.zeros([3,1])
            self.mu = np.zeros([3,1])
                
                
        def updateGuassPoints(self,element,numIncrement):
                
            self.dTheta = np.zeros([3,1]) 
            self.dr_ds = np.zeros([3,1])
            self.dTheta_ds = np.zeros([3,1])
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].dTheta
                self.dr_ds += self.dN_ds[i,0] * (element.nodes[i].historyXYZ[0]+element.nodes[i].U)
                self.dTheta_ds += self.dN_ds[i,0] * (element.nodes[i].dTheta)
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta 
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                    
                self.lambd = np.matmul(dlambd,self.lambd) 
                self.Q = np.matmul(dlambd,self.Q) 
                self.omega = np.matmul(dlambd,self.omega)
                self.omega = self.omega + (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                                           + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                                           + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
                    
            else:
                self.omega = self.omega + self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)    
                
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                                                                                                                                     #由于转动不均匀导致的当前构型下内力矩应变
            self.kappa = np.matmul(self.lambd.T,self.omega)
            
            for i in range(len(self.cableT)):
                self.cableXi[i] = np.matmul(self.lambd,self.cableXi0[i])
                self.cable_dr_ds[i] = np.array(self.dr_ds)
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
            
            self.NuB = np.matmul(self.CN,self.Gamma)                        
            self.MuB = np.matmul(self.CM,self.kappa)
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)
            
            for i in range(len(self.cableT)):
                self.Nu += self.cableF[i] * element.TLoading[i][numIncrement]
                self.Mu += self.cableM[i] * element.TLoading[i][numIncrement]
            
            self.nu = np.matmul(self.lambd,self.Nu)                                             
            self.mu = np.matmul(self.lambd,self.Mu)  
                
           
    
    def get_K_P(self,numIncrement):
        '''calculate the Nu, Mu, nu, and mu'''
        for i in self.guassPoints:
            i.Nu = np.array(i.NuB)
            i.Mu = np.array(i.MuB)    
            
            for j in range(len(i.cableF)):
                i.Nu += i.cableF[j] * self.TLoading[j][numIncrement]
                i.Mu += i.cableM[j] * self.TLoading[j][numIncrement]
        
            i.nu = np.matmul(i.lambd,i.Nu)
            i.mu = np.matmul(i.lambd,i.Mu)
        
        
        P = np.zeros([self.numDof * self.numNode, 1])                          #unbalanced nodal force
        S = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element stiffness matrix
        T = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element geometric stiffness matrix
        Kcable = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])                        #element cable stiffness matrix
        
        self.K = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])
        self.P = np.zeros([self.numDof * self.numNode, 1]) 
        
        
        
        for i in self.guassPoints:
            c = np.zeros([self.numDof,self.numDof])
            Y = []
            c[0:3,0:3] = np.matmul(i.lambd,np.matmul(i.CN,i.lambd.T))
            c[3:6,3:6] = np.matmul(i.lambd,np.matmul(i.CM,i.lambd.T))
            for j in range(len(i.cable_dr_ds_bar)):
                Y.append((np.eye(3) - np.matmul(i.cable_dr_ds_bar[j],i.cable_dr_ds_bar[j].T)) / float(sum(i.cable_dr_ds[j]**2)**0.5))
            
            for j in range(self.numNode):
                Ej = np.zeros([6,6])
                Ej[0:3,0:3] = Ej[3:6,3:6] = float(i.dN_ds[j,0]) * np.eye(3)
                Ej[3:6,0:3] = -float(i.N[j,0]) * base.V2skew(i.dr_ds)
                #print(Ej)
                
                P[j*self.numDof:(j+1)*self.numDof] -= (i.guassWeight * np.matmul(Ej,np.r_[i.nu,i.mu]) * i.detJ)
              
                for k in range(self.numNode):
                    Ek = np.zeros([6,6])
                    Ek[0:3,0:3] = Ek[3:6,3:6] = float(i.dN_ds[k,0]) * np.eye(3)
                    Ek[3:6,0:3] = -float(i.N[k,0]) * base.V2skew(i.dr_ds)
                
                    
                    
                    for l in range(len(i.cable_dr_ds)):
                        
                        Rk = np.zeros([6,6])
                        Rk[0:3,0:3] = Y[l] * float(i.dN_ds[k])
                        Rk[0:3,3:6] = (base.V2skew(i.cable_dr_ds_bar[l] * float(i.N[k])) -
                                       np.matmul(Y[l],base.V2skew(i.cableXi[l])) * float(i.dN_ds[k]) -
                                       np.matmul(Y[l],base.V2skew(np.matmul(i.OmegaT,i.cableXi[l]) + i.dXia_ds[l][0] * i.E1 + i.dXia_ds[l][1] * i.E2)) * float(i.N[k]))
                        Rk[3:6,0:3] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,0:3])
                        Rk[3:6,3:6] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,3:6])
                        Rk = self.TLoading[l][numIncrement] * i.cableT[l] * Rk  
                        
                        Kcable[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,Rk) * i.guassWeight * i.detJ
                    
                    S[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,np.matmul(c,Ek.T)) * i.guassWeight * i.detJ
                    
                    
                    T[j*self.numDof:j*self.numDof+3,k*self.numDof+3:k*self.numDof+6] +=  -base.V2skew(i.nu) * float(i.dN_ds[j]) * float(i.N[k]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof:k*self.numDof+3] +=   base.V2skew(i.nu) * float(i.dN_ds[k]) * float(i.N[j]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof+3:k*self.numDof+6] += ((-base.V2skew(i.mu) * float(i.dN_ds[j]) * float(i.N[k]) +
                                                                                            (np.matmul(i.nu,i.dr_ds.T) - np.matmul(i.nu.T,i.dr_ds) * np.eye(3)) * float(i.N[j]) * float(i.N[k]))
                                                                                           * i.guassWeight * i.detJ) 

        self.K += S + T + Kcable
        #print(Kcable)
        P += self.BF * self.FLoading[numIncrement]
        self.P += P           
                    
                
                
    def update(self,numIncrement):
        
        for i in self.guassPoints:
            i.updateGuassPoints(self,numIncrement)
                            
class eLE3DRoboticArmL3():
    
    '''linear elasticity 3D robotic arm L3 element'''
    
    def __init__(self, nodes, material, geometry, FLoading, TLoading, Reduce):
        
        self.eps = 10 ** -7
        self.numNode = 3                                                       #Number of the nodes in the element
        self.numDof = 6                                                        #Number of the degrees of freedom of the node
        
        '''input the information of the element'''
        self.nodes = nodes                                                     #Nodes of the element
        self.dofId = []                                                        #Degrees of freedom of the element
        for i in self.nodes:
            self.dofId.extend(i.dofId)
            
        self.material = material                                               #Material of the element
        self.geometry = geometry                                               #Geometry of the element
        
        self.FLoading = FLoading                                               #Load path of the body force
        self.TLoading = TLoading                                               #Load path of the tendon tension
        
        if(Reduce == 1):
            self.numGuass = [2]                                                #reduced integration 
        elif(Reduce == 0):
            self.numGuass = [3]                                                #classical integration
        else:
            self.numGuass = [4]                                                #Highest order integration
            
        '''Define material stiffness matrix'''
        self.CN0 = np.diag([material.G * geometry.A / geometry.k, 
                            material.G * geometry.A / geometry.k,
                            material.E * geometry.A]) 
        
        self.CM0 = np.diag([material.E * geometry.Ix,
                            material.E * geometry.Iy,
                            material.G * geometry.J]) 
        
        '''Define guass integration points'''
        self.guassPoints = []
        for i in range(self.numGuass[0]):
            self.guassPoints.append(self.guassPoint(self, 
                                               base.guassLocation[self.numGuass[0]-1][i],
                                               base.guassWeight[self.numGuass[0]-1][i]))
        '''Body force matrix'''
        self.BF = np.zeros([self.numNode * self.numDof, 1])
        for i in range(self.numNode):
            for j in range(len(self.guassPoints)):
                for k in range(self.numDof):
                    self.BF[i*self.numDof+k,0] += (self.guassPoints[j].guassWeight * self.guassPoints[j].f[k] * 
                                                self.guassPoints[j].N[i,0] * self.guassPoints[j].detJ)
                    
    
    '''guass intergration point class'''
    class guassPoint():
            
        def __init__(self, element, guassLocation, guassWeight):
                
            self.guassLocation = guassLocation                                 #Location of the guass intergation point                            
            self.guassWeight = guassWeight                                     #Weight of the guass intergation point
                
            '''shape function and its derivatives'''
            self.N = np.array([[(-1 + self.guassLocation) * self.guassLocation / 2],   #node 1 (-1)
                               [1 - self.guassLocation ** 2],                         #node 2 ( 0)
                               [(1 + self.guassLocation) * self.guassLocation / 2]])  #node 3 ( 1)     
                
            self.dN = np.array([[self.guassLocation - 0.5],
                                [-2.0 * self.guassLocation],                                        
                                [self.guassLocation + 0.5]])
                
            '''Jacobian matrix and its determinant'''
            self.J = np.zeros([1,1])
            dL = np.zeros([3,1])                                           
            for i in range(1):
                for j in range(1):
                    for k in range(element.numNode):
                        dL += self.dN[k,j] * element.nodes[k].historyXYZ[0]
            self.J[0,0] += float(sum(dL ** 2)) ** 0.5
            self.detJ = np.linalg.det(self.J)
            assert self.detJ>=0, 'Negative volume in element'
                
            self.dN_ds = np.matmul(self.dN,np.linalg.inv(self.J))              #Derivatives of the shape function with respect to the arc length
                
                
            self.f = np.zeros([element.numDof,1])                              #Body force in the guass intergation point
            for i in range(element.numNode):
                self.f += self.N[i,0] * element.nodes[i].f
                    

                
            
                
            '''variables with respect to configuration'''
                
            self.lambd = np.eye(3)                                             #configuration matrix
                
            self.dTheta = np.zeros([3,1])                                      #axial vector increment
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].theta
                
            self.Q = np.eye(3)                                                 #roration matrix
                
            self.dr_ds = np.zeros([3,1])                                       #Derivative of vector diameter with respect to arc length
            self.dTheta_ds = np.zeros([3,1])                                   #Derivative of axial vector increment with respect to arc length
            for i in range(element.numNode):
                self.dr_ds += self.dN_ds[i,0] * element.nodes[i].historyXYZ[0]
                self.dTheta_ds += self.dN_ds[i,0] * element.nodes[i].theta     #float(dN_ds[i][0])
                
            '''variables with respect to strain'''
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5               
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                self.Q = np.matmul(dlambd,self.Q) 
                omega0 = (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                          + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                          + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
            else:
                omega0 = self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T                #strain with respect to the displacement in current configuration 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                    #strain with respect to the displacement in reference configuration 
            self.omega = np.zeros([3,1])                                       #strain with respect to the roration in current configuration
            self.kappa = np.matmul(self.lambd.T,self.omega)                    #strain with respect to the roration in reference configuration                      
            self.Omega0 = base.V2skew(omega0)                                  #Initial curvature
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
            self.E1 = np.array([self.Q[:,0]]).T                                #Vector E1
            self.E2 = np.array([self.Q[:,1]]).T                                #Vector E2
            self.E = [self.E1,self.E2]           
                 
                    
            '''variables with respect to cable should be think about'''
            self.cableXia = []                                                 #Coordinates of cable off center
            self.dXia_ds = []                                                  #Derivatives of the Coordinates of the cable off center with respect to the arc length
            self.cableT = []                                                   #Magnitude of the tendon tension    
            for i in range(len(element.nodes[0].cableXia)):
                self.cableXia.append(np.zeros(2))
                self.dXia_ds.append(np.zeros(2))
                self.cableT.append(0)
                for j in range(element.numNode):
                    self.cableXia[i] += self.N[j,0] * element.nodes[j].cableXia[i]
                    self.dXia_ds[i] += self.dN_ds[j,0] * element.nodes[j].cableXia[i]
                    self.cableT[i] += self.N[j,0] * element.nodes[j].cableT[i]
                        
            self.cableXi0 = []
                
            for i in self.cableXia:
                self.cableXi0.append(i[0] * self.E1 + i[1] * self.E2)          #Vector of cable off center before deformation
            self.cableXi = list(self.cableXi0)   
                          #Vector of cable off center after deformation
            self.cable_dr_ds = []
            self.cable_dr_ds_bar = []
            for i in range(len(self.cableXia)):
                self.cable_dr_ds.append(np.array(self.dr_ds))
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar.append(np.zeros([3,1]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                
              
              
            '''variables with respect to interal force'''
            self.CN = np.matmul(self.Q,np.matmul(element.CN0,self.Q.T))        #the material elastic tensor in material configuration
            self.CM = np.matmul(self.Q,np.matmul(element.CM0,self.Q.T))        #the material elastic tensor in material configuration
                
            self.cableF = []                                                   #the F caused by cable
            self.cableM = []                                                   #the M caused by cable
            for i in range(len(self.cable_dr_ds_bar)):
                self.cableF.append(np.zeros([3,1]))
                self.cableM.append(np.zeros([3,1]))
                    
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
                
            self.NuB = np.matmul(self.CN,self.Gamma)                        
            self.MuB = np.matmul(self.CM,self.kappa)
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)                

                
            self.nu = np.zeros([3,1])
            self.mu = np.zeros([3,1])
                
                
        def updateGuassPoints(self,element,numIncrement):
                
            self.dTheta = np.zeros([3,1]) 
            self.dr_ds = np.zeros([3,1])
            self.dTheta_ds = np.zeros([3,1])
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].dTheta
                self.dr_ds += self.dN_ds[i,0] * (element.nodes[i].historyXYZ[0]+element.nodes[i].U)
                self.dTheta_ds += self.dN_ds[i,0] * (element.nodes[i].dTheta)
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta 
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                    
                self.lambd = np.matmul(dlambd,self.lambd) 
                self.Q = np.matmul(dlambd,self.Q) 
                self.omega = np.matmul(dlambd,self.omega)
                self.omega = self.omega + (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                                           + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                                           + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
                    
            else:
                self.omega = self.omega + self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)    
                
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                                                                                                                                     #由于转动不均匀导致的当前构型下内力矩应变
            self.kappa = np.matmul(self.lambd.T,self.omega)
            
            for i in range(len(self.cableT)):
                self.cableXi[i] = np.matmul(self.lambd,self.cableXi0[i])
                self.cable_dr_ds[i] = np.array(self.dr_ds)
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
            
            self.NuB = np.matmul(self.CN,self.Gamma)                        
            self.MuB = np.matmul(self.CM,self.kappa)
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)
            
            for i in range(len(self.cableT)):
                self.Nu += self.cableF[i] * element.TLoading[i][numIncrement]
                self.Mu += self.cableM[i] * element.TLoading[i][numIncrement]
            
            self.nu = np.matmul(self.lambd,self.Nu)                                             
            self.mu = np.matmul(self.lambd,self.Mu)  
                
           
    
    def get_K_P(self,numIncrement):
        '''calculate the Nu, Mu, nu, and mu'''
        for i in self.guassPoints:
            i.Nu = np.array(i.NuB)
            i.Mu = np.array(i.MuB)    
            
            for j in range(len(i.cableF)):
                i.Nu += i.cableF[j] * self.TLoading[j][numIncrement]
                i.Mu += i.cableM[j] * self.TLoading[j][numIncrement]
        
            i.nu = np.matmul(i.lambd,i.Nu)
            i.mu = np.matmul(i.lambd,i.Mu)
        
        
        P = np.zeros([self.numDof * self.numNode, 1])                          #unbalanced nodal force
        S = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element stiffness matrix
        T = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element geometric stiffness matrix
        Kcable = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])                        #element cable stiffness matrix
        
        self.K = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])
        self.P = np.zeros([self.numDof * self.numNode, 1]) 
        
        
        
        for i in self.guassPoints:
            c = np.zeros([self.numDof,self.numDof])
            Y = []
            c[0:3,0:3] = np.matmul(i.lambd,np.matmul(i.CN,i.lambd.T))
            c[3:6,3:6] = np.matmul(i.lambd,np.matmul(i.CM,i.lambd.T))
            for j in range(len(i.cable_dr_ds_bar)):
                Y.append((np.eye(3) - np.matmul(i.cable_dr_ds_bar[j],i.cable_dr_ds_bar[j].T)) / float(sum(i.cable_dr_ds[j]**2)**0.5))
            
            for j in range(self.numNode):
                Ej = np.zeros([6,6])
                Ej[0:3,0:3] = Ej[3:6,3:6] = float(i.dN_ds[j,0]) * np.eye(3)
                Ej[3:6,0:3] = -float(i.N[j,0]) * base.V2skew(i.dr_ds)
                #print(Ej)
                
                P[j*self.numDof:(j+1)*self.numDof] -= (i.guassWeight * np.matmul(Ej,np.r_[i.nu,i.mu]) * i.detJ)
              
                for k in range(self.numNode):
                    Ek = np.zeros([6,6])
                    Ek[0:3,0:3] = Ek[3:6,3:6] = float(i.dN_ds[k,0]) * np.eye(3)
                    Ek[3:6,0:3] = -float(i.N[k,0]) * base.V2skew(i.dr_ds)
                
                    
                    
                    for l in range(len(i.cable_dr_ds)):
                        
                        Rk = np.zeros([6,6])
                        Rk[0:3,0:3] = Y[l] * float(i.dN_ds[k])
                        Rk[0:3,3:6] = (base.V2skew(i.cable_dr_ds_bar[l] * float(i.N[k])) -
                                       np.matmul(Y[l],base.V2skew(i.cableXi[l])) * float(i.dN_ds[k]) -
                                       np.matmul(Y[l],base.V2skew(np.matmul(i.OmegaT,i.cableXi[l]) + i.dXia_ds[l][0] * i.E1 + i.dXia_ds[l][1] * i.E2)) * float(i.N[k]))
                        Rk[3:6,0:3] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,0:3])
                        Rk[3:6,3:6] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,3:6])
                        Rk = self.TLoading[l][numIncrement] * i.cableT[l] * Rk  
                        
                        Kcable[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,Rk) * i.guassWeight * i.detJ
                    
                    S[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,np.matmul(c,Ek.T)) * i.guassWeight * i.detJ
                    
                    
                    T[j*self.numDof:j*self.numDof+3,k*self.numDof+3:k*self.numDof+6] +=  -base.V2skew(i.nu) * float(i.dN_ds[j]) * float(i.N[k]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof:k*self.numDof+3] +=   base.V2skew(i.nu) * float(i.dN_ds[k]) * float(i.N[j]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof+3:k*self.numDof+6] += ((-base.V2skew(i.mu) * float(i.dN_ds[j]) * float(i.N[k]) +
                                                                                            (np.matmul(i.nu,i.dr_ds.T) - np.matmul(i.nu.T,i.dr_ds) * np.eye(3)) * float(i.N[j]) * float(i.N[k]))
                                                                                           * i.guassWeight * i.detJ) 

        self.K += S + T + Kcable
        #print(Kcable)
        P += self.BF * self.FLoading[numIncrement]
        self.P += P           
                    
                
                
    def update(self,numIncrement):
        
        for i in self.guassPoints:
            i.updateGuassPoints(self,numIncrement)                
        
            
            
class kLE3DRoboticArmL3():
    
    '''linear elasticity 3D robotic arm L3 knot'''
    
    def __init__(self, nodes, material, geometry, FLoading, TLoading, Reduce, C):
        
        self.eps = 10 ** -7
        self.numNode = 3                                                       #Number of the nodes in the element
        self.numDof = 6                                                        #Number of the degrees of freedom of the node
        
        '''input the information of the element'''
        self.nodes = nodes                                                     #Nodes of the element
        self.dofId = []                                                        #Degrees of freedom of the element
        for i in self.nodes:
            self.dofId.extend(i.dofId)
            
        self.material = material                                               #Material of the element
        self.geometry = geometry                                               #Geometry of the element
        
        self.FLoading = FLoading                                               #Load path of the body force
        self.TLoading = TLoading                                               #Load path of the tendon tension
        self.C = C                                                             #The Bezier extraction operator
        
        if(Reduce == 1):
            self.numGuass = [2]                                                #reduced integration 
        elif(Reduce == 0):
            self.numGuass = [3]                                                #classical integration
        else:
            self.numGuass = [4]                                                #Highest order integration
            
        '''Define material stiffness matrix'''
        self.CN0 = np.diag([material.G * geometry.A / geometry.k, 
                            material.G * geometry.A / geometry.k,
                            material.E * geometry.A]) 
        
        self.CM0 = np.diag([material.E * geometry.Ix,
                            material.E * geometry.Iy,
                            material.G * geometry.J]) 
        
        '''Define guass integration points'''
        self.guassPoints = []
        for i in range(self.numGuass[0]):
            self.guassPoints.append(self.guassPoint(self, 
                                               base.guassLocation[self.numGuass[0]-1][i],
                                               base.guassWeight[self.numGuass[0]-1][i]))
        '''Body force matrix'''
        self.BF = np.zeros([self.numNode * self.numDof, 1])
        for i in range(self.numNode):
            for j in range(len(self.guassPoints)):
                for k in range(self.numDof):
                    self.BF[i*self.numDof+k,0] += (self.guassPoints[j].guassWeight * self.guassPoints[j].f[k] * 
                                                self.guassPoints[j].N[i,0] * self.guassPoints[j].detJ)
                    
    
    '''guass intergration point class'''
    class guassPoint():
            
        def __init__(self, element, guassLocation, guassWeight):
                
            self.guassLocation = guassLocation                                 #Location of the guass intergation point                            
            self.guassWeight = guassWeight                                     #Weight of the guass intergation point
                
            '''shape function and its derivatives'''
            self.N = np.array([[(1 - self.guassLocation) ** 2 / 4],            #node 1 
                               [(1 - self.guassLocation ** 2) / 2],            #node 2 
                               [(1 + self.guassLocation) ** 2 / 4]])           #node 3  
            self.N = np.matmul(element.C,self.N)
                
            self.dN = np.array([[0.5*self.guassLocation - 0.5],
                                [-self.guassLocation],                                        
                                [0.5*self.guassLocation + 0.5]])
            self.dN = np.matmul(element.C,self.dN)
                
            '''Jacobian matrix and its determinant'''
            self.J = np.zeros([1,1])
            dL = np.zeros([3,1])                                           
            for i in range(1):
                for j in range(1):
                    for k in range(element.numNode):
                        dL += self.dN[k,j] * element.nodes[k].historyXYZ[0]
            self.J[0,0] += float(sum(dL ** 2)) ** 0.5
            self.detJ = np.linalg.det(self.J)
            assert self.detJ>=0, 'Negative volume in element'
                
            self.dN_ds = np.matmul(self.dN,np.linalg.inv(self.J))              #Derivatives of the shape function with respect to the arc length
                
                
            self.f = np.zeros([element.numDof,1])                              #Body force in the guass intergation point
            for i in range(element.numNode):
                self.f += self.N[i,0] * element.nodes[i].f
                    

                
            
                
            '''variables with respect to configuration'''
                
            self.lambd = np.eye(3)                                             #configuration matrix
                
            self.dTheta = np.zeros([3,1])                                      #axial vector increment
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].theta
                
            self.Q = np.eye(3)                                                 #roration matrix
                
            self.dr_ds = np.zeros([3,1])                                       #Derivative of vector diameter with respect to arc length
            self.dTheta_ds = np.zeros([3,1])                                   #Derivative of axial vector increment with respect to arc length
            for i in range(element.numNode):
                self.dr_ds += self.dN_ds[i,0] * element.nodes[i].historyXYZ[0]
                self.dTheta_ds += self.dN_ds[i,0] * element.nodes[i].theta     #float(dN_ds[i][0])
                
            '''variables with respect to strain'''
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5               
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                self.Q = np.matmul(dlambd,self.Q) 
                omega0 = (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                          + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                          + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
            else:
                omega0 = self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T                #strain with respect to the displacement in current configuration 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                    #strain with respect to the displacement in reference configuration 
            self.omega = np.zeros([3,1])                                       #strain with respect to the roration in current configuration
            self.kappa = np.matmul(self.lambd.T,self.omega)                    #strain with respect to the roration in reference configuration                      
            self.Omega0 = base.V2skew(omega0)                                  #Initial curvature
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
            self.E1 = np.array([self.Q[:,0]]).T                                #Vector E1
            self.E2 = np.array([self.Q[:,1]]).T                                #Vector E2
            self.E = [self.E1,self.E2]           
                 
                    
            '''variables with respect to cable should be think about'''
            self.cableXia = []                                                 #Coordinates of cable off center
            self.dXia_ds = []                                                  #Derivatives of the Coordinates of the cable off center with respect to the arc length
            self.cableT = []                                                   #Magnitude of the tendon tension    
            for i in range(len(element.nodes[0].cableXia)):
                self.cableXia.append(np.zeros(2))
                self.dXia_ds.append(np.zeros(2))
                self.cableT.append(0)
                for j in range(element.numNode):
                    self.cableXia[i] += self.N[j,0] * element.nodes[j].cableXia[i]
                    self.dXia_ds[i] += self.dN_ds[j,0] * element.nodes[j].cableXia[i]
                    self.cableT[i] += self.N[j,0] * element.nodes[j].cableT[i]
                        
            self.cableXi0 = []
                
            for i in self.cableXia:
                self.cableXi0.append(i[0] * self.E1 + i[1] * self.E2)          #Vector of cable off center before deformation
            self.cableXi = list(self.cableXi0)   
                          #Vector of cable off center after deformation
            self.cable_dr_ds = []
            self.cable_dr_ds_bar = []
            for i in range(len(self.cableXia)):
                self.cable_dr_ds.append(np.array(self.dr_ds))
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar.append(np.zeros([3,1]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                
              
              
            '''variables with respect to interal force'''
            self.CN = np.matmul(self.Q,np.matmul(element.CN0,self.Q.T))        #the material elastic tensor in material configuration
            self.CM = np.matmul(self.Q,np.matmul(element.CM0,self.Q.T))        #the material elastic tensor in material configuration
                
            self.cableF = []                                                   #the F caused by cable
            self.cableM = []                                                   #the M caused by cable
            for i in range(len(self.cable_dr_ds_bar)):
                self.cableF.append(np.zeros([3,1]))
                self.cableM.append(np.zeros([3,1]))
                    
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
                
            self.NuB = np.matmul(self.CN,self.Gamma)                        
            self.MuB = np.matmul(self.CM,self.kappa)
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)                

                
            self.nu = np.zeros([3,1])
            self.mu = np.zeros([3,1])
                
                
        def updateGuassPoints(self,element,numIncrement):
                
            self.dTheta = np.zeros([3,1]) 
            self.dr_ds = np.zeros([3,1])
            self.dTheta_ds = np.zeros([3,1])
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].dTheta
                self.dr_ds += self.dN_ds[i,0] * (element.nodes[i].historyXYZ[0]+element.nodes[i].U)
                self.dTheta_ds += self.dN_ds[i,0] * (element.nodes[i].dTheta)
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta 
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                    
                self.lambd = np.matmul(dlambd,self.lambd) 
                self.Q = np.matmul(dlambd,self.Q) 
                self.omega = np.matmul(dlambd,self.omega)
                self.omega = self.omega + (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                                           + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                                           + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
                    
            else:
                self.omega = self.omega + self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)    
                
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                                                                                                                                     #由于转动不均匀导致的当前构型下内力矩应变
            self.kappa = np.matmul(self.lambd.T,self.omega)
            
            for i in range(len(self.cableT)):
                self.cableXi[i] = np.matmul(self.lambd,self.cableXi0[i])
                self.cable_dr_ds[i] = np.array(self.dr_ds)
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
            
            self.NuB = np.matmul(self.CN,self.Gamma)                        
            self.MuB = np.matmul(self.CM,self.kappa)
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)
            
            for i in range(len(self.cableT)):
                self.Nu += self.cableF[i] * element.TLoading[i][numIncrement]
                self.Mu += self.cableM[i] * element.TLoading[i][numIncrement]
            
            self.nu = np.matmul(self.lambd,self.Nu)                                             
            self.mu = np.matmul(self.lambd,self.Mu)  
                
           
    
    def get_K_P(self,numIncrement):
        '''calculate the Nu, Mu, nu, and mu'''
        for i in self.guassPoints:
            i.Nu = np.array(i.NuB)
            i.Mu = np.array(i.MuB)    
            
            for j in range(len(i.cableF)):
                i.Nu += i.cableF[j] * self.TLoading[j][numIncrement]
                i.Mu += i.cableM[j] * self.TLoading[j][numIncrement]
        
            i.nu = np.matmul(i.lambd,i.Nu)
            i.mu = np.matmul(i.lambd,i.Mu)
        
        
        P = np.zeros([self.numDof * self.numNode, 1])                          #unbalanced nodal force
        S = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element stiffness matrix
        T = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element geometric stiffness matrix
        Kcable = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])                        #element cable stiffness matrix
        
        self.K = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])
        self.P = np.zeros([self.numDof * self.numNode, 1]) 
        
        
        
        for i in self.guassPoints:
            c = np.zeros([self.numDof,self.numDof])
            Y = []
            c[0:3,0:3] = np.matmul(i.lambd,np.matmul(i.CN,i.lambd.T))
            c[3:6,3:6] = np.matmul(i.lambd,np.matmul(i.CM,i.lambd.T))
            for j in range(len(i.cable_dr_ds_bar)):
                Y.append((np.eye(3) - np.matmul(i.cable_dr_ds_bar[j],i.cable_dr_ds_bar[j].T)) / float(sum(i.cable_dr_ds[j]**2)**0.5))
            
            for j in range(self.numNode):
                Ej = np.zeros([6,6])
                Ej[0:3,0:3] = Ej[3:6,3:6] = float(i.dN_ds[j,0]) * np.eye(3)
                Ej[3:6,0:3] = -float(i.N[j,0]) * base.V2skew(i.dr_ds)
                #print(Ej)
                
                P[j*self.numDof:(j+1)*self.numDof] -= (i.guassWeight * np.matmul(Ej,np.r_[i.nu,i.mu]) * i.detJ)
              
                for k in range(self.numNode):
                    Ek = np.zeros([6,6])
                    Ek[0:3,0:3] = Ek[3:6,3:6] = float(i.dN_ds[k,0]) * np.eye(3)
                    Ek[3:6,0:3] = -float(i.N[k,0]) * base.V2skew(i.dr_ds)
                
                    
                    
                    for l in range(len(i.cable_dr_ds)):
                        
                        Rk = np.zeros([6,6])
                        Rk[0:3,0:3] = Y[l] * float(i.dN_ds[k])
                        Rk[0:3,3:6] = (base.V2skew(i.cable_dr_ds_bar[l] * float(i.N[k])) -
                                       np.matmul(Y[l],base.V2skew(i.cableXi[l])) * float(i.dN_ds[k]) -
                                       np.matmul(Y[l],base.V2skew(np.matmul(i.OmegaT,i.cableXi[l]) + i.dXia_ds[l][0] * i.E1 + i.dXia_ds[l][1] * i.E2)) * float(i.N[k]))
                        Rk[3:6,0:3] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,0:3])
                        Rk[3:6,3:6] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,3:6])
                        Rk = self.TLoading[l][numIncrement] * i.cableT[l] * Rk  
                        
                        Kcable[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,Rk) * i.guassWeight * i.detJ
                    
                    S[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,np.matmul(c,Ek.T)) * i.guassWeight * i.detJ
                    
                    
                    T[j*self.numDof:j*self.numDof+3,k*self.numDof+3:k*self.numDof+6] +=  -base.V2skew(i.nu) * float(i.dN_ds[j]) * float(i.N[k]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof:k*self.numDof+3] +=   base.V2skew(i.nu) * float(i.dN_ds[k]) * float(i.N[j]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof+3:k*self.numDof+6] += ((-base.V2skew(i.mu) * float(i.dN_ds[j]) * float(i.N[k]) +
                                                                                            (np.matmul(i.nu,i.dr_ds.T) - np.matmul(i.nu.T,i.dr_ds) * np.eye(3)) * float(i.N[j]) * float(i.N[k]))
                                                                                           * i.guassWeight * i.detJ) 

        self.K += S + T + Kcable
        P += self.BF * self.FLoading[numIncrement]
        #print(Kcable)
        self.P += P           
                    
                
                
    def update(self,numIncrement):
        
        for i in self.guassPoints:
            i.updateGuassPoints(self,numIncrement)                
        
class kNH3DRoboticArmL3():
    
    '''Incompressible neo-hookean 3D robotic arm L3 knot'''
    
    def __init__(self, nodes, material, geometry, FLoading, TLoading, Reduce, C):
        
        self.eps = 10 ** -7
        self.numNode = 3                                                       #Number of the nodes in the element
        self.numDof = 6                                                        #Number of the degrees of freedom of the node
        
        '''input the information of the element'''
        self.nodes = nodes                                                     #Nodes of the element
        self.dofId = []                                                        #Degrees of freedom of the element
        for i in self.nodes:
            self.dofId.extend(i.dofId)
            
        self.material = material                                               #Material of the element
        self.geometry = geometry                                               #Geometry of the element
        
        self.FLoading = FLoading                                               #Load path of the body force
        self.TLoading = TLoading                                               #Load path of the tendon tension
        self.C = C                                                             #The Bezier extraction operator
        
        if(Reduce == 1):
            self.numGuass = [2]                                                #reduced integration 
        elif(Reduce == 0):
            self.numGuass = [3]                                                #classical integration
        else:
            self.numGuass = [4]                                                #Highest order integration
            
        '''Define material and geometric parameters'''
        self.GA = material.G * geometry.A
        self.GI1 = material.G * geometry.Ix
        self.GI2 = material.G * geometry.Iy
        self.GJ = material.G * geometry.J 
        
        '''Define guass integration points'''
        self.guassPoints = []
        for i in range(self.numGuass[0]):
            self.guassPoints.append(self.guassPoint(self, 
                                               base.guassLocation[self.numGuass[0]-1][i],
                                               base.guassWeight[self.numGuass[0]-1][i]))
        '''Body force matrix'''
        self.BF = np.zeros([self.numNode * self.numDof, 1])
        for i in range(self.numNode):
            for j in range(len(self.guassPoints)):
                for k in range(self.numDof):
                    self.BF[i*self.numDof+k,0] += (self.guassPoints[j].guassWeight * self.guassPoints[j].f[k] * 
                                                self.guassPoints[j].N[i,0] * self.guassPoints[j].detJ)
                    
    
    '''guass intergration point class'''
    class guassPoint():
            
        def __init__(self, element, guassLocation, guassWeight):
                
            self.guassLocation = guassLocation                                 #Location of the guass intergation point                            
            self.guassWeight = guassWeight                                     #Weight of the guass intergation point
                
            '''shape function and its derivatives'''
            self.N = np.array([[(1 - self.guassLocation) ** 2 / 4],            #node 1 
                               [(1 - self.guassLocation ** 2) / 2],            #node 2 
                               [(1 + self.guassLocation) ** 2 / 4]])           #node 3  
            self.N = np.matmul(element.C,self.N)
                
            self.dN = np.array([[0.5*self.guassLocation - 0.5],
                                [-self.guassLocation],                                        
                                [0.5*self.guassLocation + 0.5]])
            self.dN = np.matmul(element.C,self.dN)
                
            '''Jacobian matrix and its determinant'''
            self.J = np.zeros([1,1])
            dL = np.zeros([3,1])                                           
            for i in range(1):
                for j in range(1):
                    for k in range(element.numNode):
                        dL += self.dN[k,j] * element.nodes[k].historyXYZ[0]
            self.J[0,0] += float(sum(dL ** 2)) ** 0.5
            self.detJ = np.linalg.det(self.J)
            assert self.detJ>=0, 'Negative volume in element'
                
            self.dN_ds = np.matmul(self.dN,np.linalg.inv(self.J))              #Derivatives of the shape function with respect to the arc length
                
                
            self.f = np.zeros([element.numDof,1])                              #Body force in the guass intergation point
            for i in range(element.numNode):
                self.f += self.N[i,0] * element.nodes[i].f
                    

                
            
                
            '''variables with respect to configuration'''
            self.GA = element.GA
            self.GI1 = element.GI1
            self.GI2 = element.GI2
            self.GJ = element.GJ                
            self.lambd = np.eye(3)                                             #configuration matrix
                
            self.dTheta = np.zeros([3,1])                                      #axial vector increment
            
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].theta
                
            self.Q = np.eye(3)                                                 #roration matrix
                
            self.dr_ds = np.zeros([3,1])                                       #Derivative of vector diameter with respect to arc length
            self.dTheta_ds = np.zeros([3,1])                                   #Derivative of axial vector increment with respect to arc length
            for i in range(element.numNode):
                self.dr_ds += self.dN_ds[i,0] * element.nodes[i].historyXYZ[0]
                self.dTheta_ds += self.dN_ds[i,0] * element.nodes[i].theta     #float(dN_ds[i][0])
                
            '''variables with respect to strain'''
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5               
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                self.Q = np.matmul(dlambd,self.Q) 
                omega0 = (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                          + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                          + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
            else:
                omega0 = self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T                #strain with respect to the displacement in current configuration 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                    #strain with respect to the displacement in reference configuration 
            self.omega = np.zeros([3,1])                                       #strain with respect to the roration in current configuration
            self.kappa = np.matmul(self.lambd.T,self.omega)                    #strain with respect to the roration in reference configuration                      
            self.Omega0 = base.V2skew(omega0)                                  #Initial curvature
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
            self.Q0 = np.array(self.Q)                                         #Initial roration matrix
            self.E1 = np.array([self.Q[:,0]]).T                                #Vector E1
            self.E2 = np.array([self.Q[:,1]]).T                                #Vector E2
            self.E3 = np.array([self.Q[:,2]]).T                                #Vector E3
            self.E = [self.E1,self.E2]  
            self.w1 = float(np.matmul(self.E1.T,self.kappa))
            self.w2 = float(np.matmul(self.E2.T,self.kappa))
            self.w3 = float(np.matmul(self.E3.T,self.kappa))
            self.T1 = float(np.matmul(self.E1.T,self.Gamma))
            self.T2 = float(np.matmul(self.E2.T,self.Gamma))
            self.T3 = float(np.matmul(self.E3.T,self.Gamma)) + 1.0
                 
                    
            '''variables with respect to cable should be think about'''
            self.cableXia = []                                                 #Coordinates of cable off center
            self.dXia_ds = []                                                  #Derivatives of the Coordinates of the cable off center with respect to the arc length
            self.cableT = []                                                   #Magnitude of the tendon tension    
            for i in range(len(element.nodes[0].cableXia)):
                self.cableXia.append(np.zeros(2))
                self.dXia_ds.append(np.zeros(2))
                self.cableT.append(0)
                for j in range(element.numNode):
                    self.cableXia[i] += self.N[j,0] * element.nodes[j].cableXia[i]
                    self.dXia_ds[i] += self.dN_ds[j,0] * element.nodes[j].cableXia[i]
                    self.cableT[i] += self.N[j,0] * element.nodes[j].cableT[i]
                        
            self.cableXi0 = []
                
            for i in self.cableXia:
                self.cableXi0.append(i[0] * self.E1 + i[1] * self.E2)          #Vector of cable off center before deformation
            self.cableXi = list(self.cableXi0)   
                          #Vector of cable off center after deformation
            self.cable_dr_ds = []
            self.cable_dr_ds_bar = []
            for i in range(len(self.cableXia)):
                self.cable_dr_ds.append(np.array(self.dr_ds))
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar.append(np.zeros([3,1]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                
              
              
            '''variables with respect to interal force'''
            self.c = np.zeros([6,6])
            self.c[0,0] = float(self.GA)
            self.c[0,2] = 2 * self.GI1 * self.w1 * self.w3 / (self.T3 ** 3)
            self.c[0,3] = - self.GI1 * self.w3 / (self.T3 ** 2)
            self.c[0,5] = - self.GI1 * self.w1 / (self.T3 ** 2)
            self.c[1,1] = float(self.GA)
            self.c[1,2] = 2 * self.GI2 * self.w2 * self.w3 / (self.T3 ** 3)
            self.c[1,4] = - self.GI2 * self.w3 / (self.T3 ** 2)
            self.c[1,5] = - self.GI2 * self.w2 / (self.T3 ** 2)
            self.c[2,2] = (1 + 2 * (self.T3 ** -3)) * self.GA - 2 * (self.GI1 * (self.w1 ** 2) + self.GI2 * (self.w2 ** 2)) / (self.T3 ** 3)
            self.c[2,3] = 2 * self.GI1 * self.w1 / (self.T3 ** 2)
            self.c[2,4] = 2 * self.GI2 * self.w2 / (self.T3 ** 2)
            self.c[3,2] = -(2 * (self.T3 ** -2) + 4 * (self.T3 ** -5)) * self.GI1 * self.w1
            self.c[3,3] = self.GI1 * (2 * self.T3 ** (-1) + self.T3 ** (-4))
            self.c[4,2] = -(2 * (self.T3 ** -2) + 4 * (self.T3 ** -5)) * self.GI2 * self.w2
            self.c[4,4] = self.GI2 * (2 * (self.T3 ** -1) + (self.T3 ** -4))
            self.c[5,0] = -(self.GI1 * self.w1) / (self.T3 ** 2)
            self.c[5,1] = -(self.GI2 * self.w2) / (self.T3 ** 2)
            self.c[5,2] = -(self.GJ * self.w3) / (self.T3 ** 2)
            self.c[5,5] = self.GJ / self.T3
            
            R = np.zeros([6,6])
            R[0:3,0:3] = R[3:6,3:6] = np.array(self.Q0)
            self.c = np.matmul(R,np.matmul(self.c,R.T))
            
            
            
            
                
            self.cableF = []                                                   #the F caused by cable
            self.cableM = []                                                   #the M caused by cable
            for i in range(len(self.cable_dr_ds_bar)):
                self.cableF.append(np.zeros([3,1]))
                self.cableM.append(np.zeros([3,1]))
                    
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
            
            self.F = np.zeros([6,1])
            self.F[0,0] = self.GA * self.T1 - self.GI1 * self.w1 * self.w3 / (self.T3 ** 2)
            self.F[1,0] = self.GA * self.T2 - self.GI2 * self.w2 * self.w3 / (self.T3 ** 2)
            self.F[2,0] = self.GA * (self.T3 - self.T3 ** (-2)) + (self.GI1 * (self.w1 ** 2) + self.GI2 * (self.w2 ** 2)) / (self.T3 ** 2)
            self.F[3,0] = self.GI1 * (2 * self.T3 ** (-1) + self.T3 ** (-4)) * self.w1
            self.F[4,0] = self.GI2 * (2 * self.T3 ** (-1) + self.T3 ** (-4)) * self.w2
            self.F[5,0] = (self.GJ * self.w3 * self.T3 - self.GI2 * self.w2 * self.T2 - self.GI1 * self.w1 * self.T1) / (self.T3 ** 2)
            
            
            
            self.NuB = np.matmul(self.Q0,self.F[0:3])                        
            self.MuB = np.matmul(self.Q0,self.F[3:6])
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)                

                
            self.nu = np.zeros([3,1])
            self.mu = np.zeros([3,1])
                
                
        def updateGuassPoints(self,element,numIncrement):
                
            self.dTheta = np.zeros([3,1]) 
            self.dr_ds = np.zeros([3,1])
            self.dTheta_ds = np.zeros([3,1])
            for i in range(element.numNode):
                self.dTheta += self.N[i,0] * element.nodes[i].dTheta
                self.dr_ds += self.dN_ds[i,0] * (element.nodes[i].historyXYZ[0]+element.nodes[i].U)
                self.dTheta_ds += self.dN_ds[i,0] * (element.nodes[i].dTheta)
                
            abs_dTheta = float(sum(self.dTheta ** 2)) ** 0.5
            if abs_dTheta >= element.eps:
                bar_dTheta = m.tan(abs_dTheta/2) * self.dTheta / abs_dTheta 
                dlambd = (np.eye(3) + 
                          2  * (base.V2skew(bar_dTheta) + np.matmul(base.V2skew(bar_dTheta),base.V2skew(bar_dTheta))) 
                          / (1 + m.tan(abs_dTheta/2) ** 2))
                    
                self.lambd = np.matmul(dlambd,self.lambd) 
                self.Q = np.matmul(dlambd,self.Q) 
                self.omega = np.matmul(dlambd,self.omega)
                self.omega = self.omega + (m.sin(abs_dTheta) * self.dTheta_ds / abs_dTheta
                                           + (1 - m.sin(abs_dTheta) / abs_dTheta) * (np.matmul(self.dTheta.T,self.dTheta_ds) / abs_dTheta) * self.dTheta / abs_dTheta
                                           + 0.5 * (2 * m.sin(0.5 * abs_dTheta) / abs_dTheta) ** 2 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds))
                    
            else:
                self.omega = self.omega + self.dTheta_ds + 0.5 * np.matmul(base.V2skew(self.dTheta),self.dTheta_ds)    
                
            self.OmegaT = base.V2skew(self.omega) + np.matmul(self.lambd,np.matmul(self.Omega0,self.lambd.T))
                
            self.gamma = self.dr_ds - np.array([self.Q[:,2]]).T 
            self.Gamma = np.matmul(self.lambd.T,self.gamma)                                                                                                                                     #由于转动不均匀导致的当前构型下内力矩应变
            self.kappa = np.matmul(self.lambd.T,self.omega)
            self.w1 = float(np.matmul(self.E1.T,self.kappa))
            self.w2 = float(np.matmul(self.E2.T,self.kappa))
            self.w3 = float(np.matmul(self.E3.T,self.kappa))
            self.T1 = float(np.matmul(self.E1.T,self.Gamma))
            self.T2 = float(np.matmul(self.E2.T,self.Gamma))
            self.T3 = float(np.matmul(self.E3.T,self.Gamma)) + 1.0
            
            for i in range(len(self.cableT)):
                self.cableXi[i] = np.matmul(self.lambd,self.cableXi0[i])
                self.cable_dr_ds[i] = np.array(self.dr_ds)
                for j in range(2):
                    self.cable_dr_ds[i] += np.matmul((self.cableXia[i][j]*self.OmegaT+self.dXia_ds[i][j]*np.eye(3)),np.matmul(self.lambd,self.E[j]))
                self.cable_dr_ds_bar[i] = self.cable_dr_ds[i] / float(sum(self.cable_dr_ds[i]**2)**0.5)
                self.cableF[i] = np.matmul(self.lambd.T,self.cableT[i] * self.cable_dr_ds_bar[i])
                self.cableM[i] = np.matmul(base.V2skew(self.cableXi0[i]),self.cableF[i])
            
            self.c = np.zeros([6,6])
            self.c[0,0] = float(self.GA)
            self.c[0,2] = 2 * self.GI1 * self.w1 * self.w3 / (self.T3 ** 3)
            self.c[0,3] = - self.GI1 * self.w3 / (self.T3 ** 2)
            self.c[0,5] = - self.GI1 * self.w1 / (self.T3 ** 2)
            self.c[1,1] = float(self.GA)
            self.c[1,2] = 2 * self.GI2 * self.w2 * self.w3 / (self.T3 ** 3)
            self.c[1,4] = - self.GI2 * self.w3 / (self.T3 ** 2)
            self.c[1,5] = - self.GI2 * self.w2 / (self.T3 ** 2)
            self.c[2,2] = (1 + 2 * (self.T3 ** -3)) * self.GA - 2 * (self.GI1 * (self.w1 ** 2) + self.GI2 * (self.w2 ** 2)) / (self.T3 ** 3)
            self.c[2,3] = 2 * self.GI1 * self.w1 / (self.T3 ** 2)
            self.c[2,4] = 2 * self.GI2 * self.w2 / (self.T3 ** 2)
            self.c[3,2] = -(2 * (self.T3 ** -2) + 4 * (self.T3 ** -5)) * self.GI1 * self.w1
            self.c[3,3] = self.GI1 * (2 * self.T3 ** (-1) + self.T3 ** (-4))
            self.c[4,2] = -(2 * (self.T3 ** -2) + 4 * (self.T3 ** -5)) * self.GI2 * self.w2
            self.c[4,4] = self.GI2 * (2 * (self.T3 ** -1) + (self.T3 ** -4))
            self.c[5,0] = -(self.GI1 * self.w1) / (self.T3 ** 2)
            self.c[5,1] = -(self.GI2 * self.w2) / (self.T3 ** 2)
            self.c[5,2] = -(self.GJ * self.w3) / (self.T3 ** 2)
            self.c[5,5] = self.GJ / self.T3
            R = np.zeros([6,6])
            R[0:3,0:3] = R[3:6,3:6] = np.array(self.Q0)
            self.c = np.matmul(R,np.matmul(self.c,R.T))
            
            self.F[0,0] = self.GA * self.T1 - self.GI1 * self.w1 * self.w3 / (self.T3 ** 2)
            self.F[1,0] = self.GA * self.T2 - self.GI2 * self.w2 * self.w3 / (self.T3 ** 2)
            self.F[2,0] = self.GA * (self.T3 - self.T3 ** (-2)) + (self.GI1 * (self.w1 ** 2) + self.GI2 * (self.w2 ** 2)) / (self.T3 ** 2)
            self.F[3,0] = self.GI1 * (2 * self.T3 ** (-1) + self.T3 ** (-4)) * self.w1
            self.F[4,0] = self.GI2 * (2 * self.T3 ** (-1) + self.T3 ** (-4)) * self.w2
            self.F[5,0] = (self.GJ * self.w3 * self.T3 - self.GI2 * self.w2 * self.T2 - self.GI1 * self.w1 * self.T1) / (self.T3 ** 2)
            
            self.NuB = np.matmul(self.Q0,self.F[0:3])                        
            self.MuB = np.matmul(self.Q0,self.F[3:6])
            self.Nu = np.array(self.NuB)
            self.Mu = np.array(self.MuB)
            
            for i in range(len(self.cableT)):
                self.Nu += self.cableF[i] * element.TLoading[i][numIncrement]
                self.Mu += self.cableM[i] * element.TLoading[i][numIncrement]
            
            self.nu = np.matmul(self.lambd,self.Nu)                                             
            self.mu = np.matmul(self.lambd,self.Mu)  
                
           
    
    def get_K_P(self,numIncrement):
        '''calculate the Nu, Mu, nu, and mu'''
        for i in self.guassPoints:
            i.Nu = np.array(i.NuB)
            i.Mu = np.array(i.MuB)    
            
            for j in range(len(i.cableF)):
                i.Nu += i.cableF[j] * self.TLoading[j][numIncrement]
                i.Mu += i.cableM[j] * self.TLoading[j][numIncrement]
        
            i.nu = np.matmul(i.lambd,i.Nu)
            i.mu = np.matmul(i.lambd,i.Mu)
        
        
        P = np.zeros([self.numDof * self.numNode, 1])                          #unbalanced nodal force
        S = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element stiffness matrix
        T = np.zeros([self.numDof * self.numNode, self.numDof * self.numNode]) #element geometric stiffness matrix
        Kcable = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])                        #element cable stiffness matrix
        
        self.K = np.zeros([self.numDof * self.numNode, 
                           self.numDof * self.numNode])
        self.P = np.zeros([self.numDof * self.numNode, 1]) 
        
        
        
        for i in self.guassPoints:
            c = np.zeros([self.numDof,self.numDof])
            Y = []
            PI = np.zeros([self.numDof,self.numDof])
            PI[0:3,0:3] = np.array(i.lambd)
            PI[3:6,3:6] = np.array(i.lambd)
            c = np.matmul(PI,np.matmul(i.c,PI.T))
    
            for j in range(len(i.cable_dr_ds_bar)):
                Y.append((np.eye(3) - np.matmul(i.cable_dr_ds_bar[j],i.cable_dr_ds_bar[j].T)) / float(sum(i.cable_dr_ds[j]**2)**0.5))
            
            for j in range(self.numNode):
                Ej = np.zeros([6,6])
                Ej[0:3,0:3] = Ej[3:6,3:6] = float(i.dN_ds[j,0]) * np.eye(3)
                Ej[3:6,0:3] = -float(i.N[j,0]) * base.V2skew(i.dr_ds)
                #print(Ej)
                
                P[j*self.numDof:(j+1)*self.numDof] -= (i.guassWeight * np.matmul(Ej,np.r_[i.nu,i.mu]) * i.detJ)
              
                for k in range(self.numNode):
                    Ek = np.zeros([6,6])
                    Ek[0:3,0:3] = Ek[3:6,3:6] = float(i.dN_ds[k,0]) * np.eye(3)
                    Ek[3:6,0:3] = -float(i.N[k,0]) * base.V2skew(i.dr_ds)
                
                    
                    
                    for l in range(len(i.cable_dr_ds)):
                        
                        Rk = np.zeros([6,6])
                        Rk[0:3,0:3] = Y[l] * float(i.dN_ds[k])
                        Rk[0:3,3:6] = (base.V2skew(i.cable_dr_ds_bar[l] * float(i.N[k])) -
                                       np.matmul(Y[l],base.V2skew(i.cableXi[l])) * float(i.dN_ds[k]) -
                                       np.matmul(Y[l],base.V2skew(np.matmul(i.OmegaT,i.cableXi[l]) + i.dXia_ds[l][0] * i.E1 + i.dXia_ds[l][1] * i.E2)) * float(i.N[k]))
                        Rk[3:6,0:3] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,0:3])
                        Rk[3:6,3:6] = np.matmul(base.V2skew(i.cableXi[l]),Rk[0:3,3:6])
                        Rk = self.TLoading[l][numIncrement] * i.cableT[l] * Rk  
                        
                        Kcable[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,Rk) * i.guassWeight * i.detJ
                    
                    S[j*self.numDof:(j+1)*self.numDof,k*self.numDof:(k+1)*self.numDof] += np.matmul(Ej,np.matmul(c,Ek.T)) * i.guassWeight * i.detJ
                    
                    
                    T[j*self.numDof:j*self.numDof+3,k*self.numDof+3:k*self.numDof+6] +=  -base.V2skew(i.nu) * float(i.dN_ds[j]) * float(i.N[k]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof:k*self.numDof+3] +=   base.V2skew(i.nu) * float(i.dN_ds[k]) * float(i.N[j]) * i.guassWeight * i.detJ
                    T[j*self.numDof+3:j*self.numDof+6,k*self.numDof+3:k*self.numDof+6] += ((-base.V2skew(i.mu) * float(i.dN_ds[j]) * float(i.N[k]) +
                                                                                            (np.matmul(i.nu,i.dr_ds.T) - np.matmul(i.nu.T,i.dr_ds) * np.eye(3)) * float(i.N[j]) * float(i.N[k]))
                                                                                           * i.guassWeight * i.detJ) 

        self.K += S + T + Kcable
        #print(Kcable)
        P += self.BF * self.FLoading[numIncrement]
        self.P += P           
                    
                
                
    def update(self,numIncrement):
        
        for i in self.guassPoints:
            i.updateGuassPoints(self,numIncrement)                            
                            
                

                                    
                    
                
                    
                
                
                
            
        