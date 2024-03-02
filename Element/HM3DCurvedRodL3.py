# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 16:44:37 2023

@author: XinLi
"""



import numpy as np
import base
import math as m

class HM3DCurvedRodL3():
    
    '''3D hard-magnetic curved linear elastic rod L3 element'''
    
    def __init__(self,nodes,reduce):
        
        self.eps = 1e-7
        self.numNodes = 3                                                      # number of the nodes in the element
        self.numDof = 6                                                        # number of the degree of freedom of the node
        
        '''input the information of the elements'''
        
        self.nodes = nodes                                                     # nodes list
        self.dofID = []                                                        # degree of freedom of the element
        for i in nodes:
            self.dofID.extend(i.dofID)
            
        if(reduce == 1):
            self.numGuass = [2]                                                # reduced integration
        elif(reduce == 0):
            self.numGuass = [3]                                                # classical integration
        else:
            self.numGuass = [4]                                                # Highest order integration
            
        '''Define guass integration points'''
        
        self.guassPoints = []
        for i in range(self.numGuass[0]):
            self.guassPoints.append(self.guassPoint(self, 
                                                    base.guassLocation[self.numGuass[0]-1][i],
                                                    base.guassWeight[self.numGuass[0]-1][i]))
    
    def get_K_P(self,dt,beta,gamma):
        
        self.K = np.zeros([self.numDof * self.numNodes, 
                           self.numDof * self.numNodes])                       # element stiffiness matrix
        self.P = np.zeros([self.numDof * self.numNodes, 1])                    # element out-of-balance matrix
        
        self.K1 = np.zeros([self.numDof * self.numNodes, 
                           self.numDof * self.numNodes])                       # classical element stiffiness matrix
        self.Ka = np.zeros([self.numDof * self.numNodes, 
                           self.numDof * self.numNodes])                       # applied force stiffiness matrix
        self.Kd = np.zeros([self.numDof * self.numNodes, 
                           self.numDof * self.numNodes])                       # dynamical stiffiness matrix
        
        for g in self.guassPoints:
            
            PI = np.zeros([6,6])
            PI[0:3,0:3] = np.array(g.lambd1)
            PI[3:6,3:6] = np.array(g.lambd1)
            d = np.matmul(PI,np.matmul(g.D,PI.T))
            Y = np.zeros([9,9])
            Y[0:3,6:9] = -base.V2skew(g.nu)
            Y[3:6,6:9] = -base.V2skew(g.mu)
            Y[6:9,0:3] =  base.V2skew(g.nu)
            Y[6:9,6:9] =  np.matmul(g.nu,g.dR_ds.T) - np.matmul(g.nu.T,g.dR_ds) * np.eye(3)
            
            qg = np.zeros([6,6])
            qg[3:6,3:6] = g.density * g.A * np.matmul(base.V2skew(g.g),base.V2skew(g.Rc))
            qm = np.zeros([6,6])
            RMt = np.matmul(g.lambd1,g.RM)
            qm[0:3,3:6] = -np.matmul(g.AM[1].T,base.V2skew(RMt))
            qm[3:6,0:3] =  np.matmul(base.V2skew(RMt),g.AM[1])
            qm[3:6,3:6] =  np.matmul(base.V2skew(g.AM0),base.V2skew(RMt))
            qm = qm * g.A / g.Gm
            
            H1 = np.zeros([6,6])
            H2 = np.zeros([6,6])
            
            abstheta = float(sum(g.theta ** 2)) ** 0.5
            e2 = abstheta / 2
            if abstheta >= self.eps:
                e = g.theta / abstheta
                Z = np.matmul(e,e.T) + (e2 / m.tan(e2)) * (np.eye(3) - np.matmul(e,e.T)) - base.V2skew(g.theta) / 2
            else:
                Z = np.eye(3) - base.V2skew(g.theta) / 2
            
            
            H1[0:3,0:3] = g.Arho * np.eye(3) / (beta * dt**2)
            a = g.Jrho / (beta * dt**2) + gamma * (np.matmul(base.V2skew(g.W),g.Jrho) - base.V2skew(np.matmul(g.Jrho,g.W))) / (beta * dt) 
            b = np.matmul(g.Jrho,g.dW) + np.matmul(base.V2skew(g.W),np.matmul(g.Jrho,g.W))
            H1[3:6,3:6] =  np.matmul(np.matmul(g.lambd1,a),np.matmul(g.lambd0.T,Z)) - base.V2skew(np.matmul(g.lambd1,b))
            
            H2[3:6,0:3] = g.Arho * base.V2skew(g.Rc) / (beta * dt**2)
            H2[3:6,0:3] = g.Arho * np.matmul(base.V2skew(g.ddR),base.V2skew(g.Rc))
            a = base.V2skew(g.Rc0) / (beta * dt**2) + gamma * (base.V2skew(np.matmul(base.V2skew(g.W),g.Rc0)) + np.matmul(base.V2skew(g.W),base.V2skew(g.Rc0))) / (beta * dt)
            b = np.matmul(base.V2skew(g.dW),g.Rc0) + np.matmul(base.V2skew(g.W),np.matmul(base.V2skew(g.W),g.Rc0))
            H2[0:3,3:6] = -(np.matmul(np.matmul(g.lambd1,a),np.matmul(g.lambd0.T,Z)) + base.V2skew(np.matmul(g.lambd1,b))) * g.Arho
            
            for i in range(self.numNodes):
                EI = np.zeros([6,6])
                EI[0:3,0:3] = EI[3:6,3:6] = float(g.dN_ds[i,0]) * np.eye(3)
                EI[3:6,0:3] = -float(g.N[i,0]) * base.V2skew(g.dR_ds)
                
                YI = np.zeros([6,9])
                YI[0:3,0:3] = YI[3:6,3:6] = float(g.dN_ds[i,0]) * np.eye(3)
                YI[3:6,6:9] = float(g.N[i,0]) * np.eye(3)
                
                self.P[i*self.numDof:(i+1)*self.numDof] += (g.guassWeight * float(g.N[i,0]) * (g.f + g.fm + g.fg - g.fd) * g.detJ) 
                self.P[i*self.numDof:(i+1)*self.numDof] -= (g.guassWeight * np.matmul(EI,np.r_[g.nu,g.mu]) * g.detJ)
                
                for j in range(self.numNodes):
                    EJ = np.zeros([6,6])
                    EJ[0:3,0:3] = EJ[3:6,3:6] = float(g.dN_ds[j,0]) * np.eye(3)
                    EJ[3:6,0:3] = -float(g.N[j,0]) * base.V2skew(g.dR_ds)
                    
                    YJ = np.zeros([6,9])
                    YJ[0:3,0:3] = YJ[3:6,3:6] = float(g.dN_ds[j,0]) * np.eye(3)
                    YJ[3:6,6:9] = float(g.N[j,0]) * np.eye(3)
                    
                    self.K1[i*self.numDof:(i+1)*self.numDof,j*self.numDof:(j+1)*self.numDof] += ((np.matmul(EI,np.matmul(d,EJ.T)) 
                                                                                                 + np.matmul(YI,np.matmul(Y,YJ.T)))
                                                                                                 * g.guassWeight * g.detJ)
                    
                    self.Ka[i*self.numDof:(i+1)*self.numDof,j*self.numDof:(j+1)*self.numDof] -= ((qg + qm) * g.guassWeight * g.detJ * float(g.N[i,0]) * float(g.N[j,0]))
                    self.Kd[i*self.numDof:(i+1)*self.numDof,j*self.numDof:(j+1)*self.numDof] += ((H1 + H2) * g.guassWeight * g.detJ * float(g.N[i,0]) * float(g.N[j,0]))
            
        self.K = self.K1 + self.Ka + self.Kd
                    
    
    
    
    
    
    def timeStepUpdate(self,dt,beta,gamma):
        for i in self.guassPoints:
            i.timeStepUpdate(self,dt,beta,gamma)
    
    def iterationUpdate(self,dt,beta,gamma):
        for i in self.guassPoints:
            i.iterationUpdate(self,dt,beta,gamma)

        
    
    '''guass intergration point class'''
    class guassPoint():
        
        def __init__(self,element,guassLocation,guassWeight):
            
            self.guassLocation = guassLocation                                 # location of the guass intergation point                            
            self.guassWeight = guassWeight                                     # weight of the guass intergation point
            
            '''shape function and its derivatives'''
            self.N = np.array([[(1 - self.guassLocation) ** 2 / 4],            # node 1 
                               [(1 - self.guassLocation ** 2) / 2],            # node 2 
                               [(1 + self.guassLocation) ** 2 / 4]])           # node 3
            
            self.dN = np.array([[0.5*self.guassLocation - 0.5],
                                [-self.guassLocation],                                        
                                [0.5*self.guassLocation + 0.5]])
            
            '''Jacobian matrix and its determinant'''
            self.Ja = np.zeros([1,1])
            dL = np.zeros([3,1])
            for i in range(1):
                for j in range(1):
                    for k in range(element.numNodes):
                        dL += self.dN[k] * element.nodes[k].R[0]
            self.Ja[0,0] += float(sum(dL ** 2)) ** 0.5
            self.detJ = np.linalg.det(self.Ja)
            assert self.detJ>0, 'Negative volume in element'
            
            self.dN_ds = np.matmul(self.dN,np.linalg.inv(self.Ja))             # derivatives of the shape function with respect to the arc length
            
            '''the quantities never change in the calculation process'''
            
            self.E = 0                                                         # Young's modulus
            self.G = 0                                                         # shear modulus
            self.Gm = 4*m.pi * 1e-7                                            # magnetic permeability 
            self.density = 0                                                   # density
            self.A = 0                                                         # area of the cross-section
            self.k = np.array([[0.0],[0.0]])                                   # reasonable factor
            self.I1 = 0                                                        # moment of inertia of the cross-section at 1-dirction
            self.I2 = 0                                                        # moment of inertia of the cross-section at 2-dirction
            self.RM = 0                                                        # residual magnetic flux density
            self.iniAxiVector = 0                                              # initial axial vector
            self.KI = 0
            
            for i in range(element.numNodes):
                self.E += self.N[i,0] * element.nodes[i].material.E
                self.G += self.N[i,0] * element.nodes[i].material.G
                self.density += self.N[i,0] * element.nodes[i].material.density
                self.A += self.N[i,0] * element.nodes[i].crossSection.A
                self.k += self.N[i,0] * element.nodes[i].crossSection.k
                self.I1 += self.N[i,0] * element.nodes[i].crossSection.I1
                self.I2 += self.N[i,0] * element.nodes[i].crossSection.I2
                self.RM += self.N[i,0] * element.nodes[i].RM
                self.iniAxiVector += self.N[i,0] * element.nodes[i].iniAxiVector
                self.KI += self.N[i,0] * element.nodes[i].KI
            self.J = self.I1 + self.I2
            self.Arho = self.A * self.density * self.KI
            
            self.Q0 = base.theta2Q(self.iniAxiVector)                          # initial rotation matrix
            self.E1 = np.array([self.Q0[:,0]]).T                               # vector E1
            self.E2 = np.array([self.Q0[:,1]]).T                               # vector E2
            self.E3 = np.array([self.Q0[:,2]]).T                               # vector E3
            
            self.dR0_ds = np.zeros([3,1])                                      # initial derivative of vector radius  
            dAxi_ds = np.zeros([3,1])
            
            for i in range(element.numNodes):
                self.dR0_ds += self.dN_ds[i,0] * element.nodes[i].R[0]
                dAxi_ds += self.dN_ds[i,0] * element.nodes[i].iniAxiVector
            
            self.dR_ds = np.array(self.dR0_ds)
                
            abs_iniAxiVector = float(sum(self.iniAxiVector ** 2)) ** 0.5
            if abs_iniAxiVector >= element.eps:
                omega0 = (m.sin(abs_iniAxiVector) * dAxi_ds / abs_iniAxiVector
                          + (1 - m.sin(abs_iniAxiVector) / abs_iniAxiVector) * (np.matmul(self.iniAxiVector.T,dAxi_ds) / abs_iniAxiVector) * self.iniAxiVector / abs_iniAxiVector
                          + 0.5 * (2 * m.sin(0.5 * abs_iniAxiVector) / abs_iniAxiVector) ** 2 * np.matmul(base.V2skew(self.iniAxiVector),dAxi_ds))
            else:
                omega0 = dAxi_ds + 0.5 * np.matmul(base.V2skew(self.iniAxiVector),dAxi_ds)
            
            self.omega01 = float(np.matmul(omega0.T,self.E1))                  # initial curvature in E1 dirction
            self.omega02 = float(np.matmul(omega0.T,self.E2))                  # initial curvature in E2 dirction
            self.Rc0 = (-self.omega02 * self.I2 * self.E1 / self.A
                        +self.omega01 * self.I1 * self.E2 / self.A)            # initial relative coordinates of the mass of centers 
            self.Rc = np.array(self.Rc0)
            
            Jrho = np.diag([self.density * self.I1,
                            self.density * self.I2,
                            self.density * self.J])
            Jrho = Jrho * self.KI
            self.Jrho = np.matmul(self.Q0,np.matmul(Jrho,self.Q0.T))
            
            
            R = np.zeros([6,6])
            R[0:3,0:3] = R[3:6,3:6] = np.array(self.Q0)
            D = np.diag([self.G * self.A * self.k[0,0],
                         self.G * self.A * self.k[1,0],
                         self.E * self.A, 
                         self.E * self.I1,
                         self.E * self.I2,
                         self.G * self.J])
            D[0,5] = D[5,0] = self.omega01 * self.G * self.I1
            D[1,5] = D[5,1] = self.omega02 * self.G * self.I2
            D[2,3] = D[3,2] =-self.omega01 * self.E * self.I1
            D[2,4] = D[4,2] =-self.omega02 * self.E * self.I2
            self.D = np.matmul(R,np.matmul(D,R.T))
            
            
            '''configuration quantities'''
            self.R = np.zeros([3,1])
            self.dR = np.zeros([3,1])
            self.ddR = np.zeros([3,1])
            self.W = np.zeros([3,1])
            self.dW = np.zeros([3,1])
            for i in range(element.numNodes):
                self.R += self.N[i,0] * element.nodes[i].R[-1]
                self.dR += self.N[i,0] * element.nodes[i].dR[-1]
                self.ddR += self.N[i,0] * element.nodes[i].ddR[-1]
                self.W += self.N[i,0] * element.nodes[i].W[-1] 
                self.dW += self.N[i,0] * element.nodes[i].dW[-1] 
            
            self.lambd0 = np.eye(3)                                            # rotation matrix at the time step start
            self.lambd1 = np.eye(3)                                            # rotation matrix in the time step
            
            '''strain and stress'''
            self.omega = np.zeros([3,1])                                       # rotational strain in the spatial description
            self.Strain = np.zeros([6,1])                                      # strain in the material description
            self.Stress = np.zeros([6,1])                                      # stress in the material description
            self.nu = np.zeros([3,1])                                          # translational stress in the spatial description
            self.mu = np.zeros([3,1])                                          # rotational stress in the spatial description   
            
            

        def timeStepUpdate(self,element,dt,beta,gamma):
            
            '''update the configuration quantities'''
            ddR = -self.dR / (dt * beta) - (0.5 - beta) * self.ddR / beta
            dR = self.dR + dt * ((1.0-gamma) * self.ddR + gamma * ddR)
            self.dR = np.array(dR)                                             # linear velocity
            self.ddR = np.array(ddR)                                           # linear acceleration
            
            dW = -self.W / (dt * beta) - (0.5 - beta) * self.dW / beta
            W = self.W + dt * ((1.0-gamma) * self.dW + gamma * dW)
            self.dW = np.array(dW)                                             # angular velocity
            self.W = np.array(W)                                               # angular acceleration
            self.Theta = np.zeros([3,1])                                       # material angular displacement vector in time step
            self.theta = np.zeros([3,1])                                       # angular displacement vector in time step
            
            
            self.axiVector = base.Q2theta(self.lambd1)                         # axial vector
            self.lambd0 = np.array(self.lambd1)
            
            '''update the applied quantities'''
            self.AM = [np.zeros([3,1]),np.zeros([3,3])]                        # applied magnetic flux density
            self.f = np.zeros([6,1])                                           # distributed force
            self.g = np.zeros([3,1])                                           # body force density
            
            for i in range(element.numNodes):
                self.AM[0] += self.N[i,0] * element.nodes[i].AM[0]
                self.AM[1] += self.N[i,0] * element.nodes[i].AM[1]
                self.f += self.N[i,0] * element.nodes[i].f
                self.g += self.N[i,0] * element.nodes[i].g
            
            self.AM0 = self.AM[0] + np.matmul(self.AM[1],self.R)
            self.fm = np.zeros([6,1])                                          # magnetic force
            self.fg = np.zeros([6,1])                                          # body force
            
            self.fm[0:3] = self.A * np.matmul(self.AM[1].T,np.matmul(self.lambd1,self.RM)) / self.Gm
            self.fm[3:6] = self.A * np.matmul(base.V2skew(np.matmul(self.lambd1,self.RM)),self.AM0) / self.Gm
            self.fg[0:3] = self.density * self.A * self.g
            self.fg[3:6] = self.density * self.A * np.matmul(base.V2skew(self.Rc),self.g)
            
            '''update the dynamic force'''
            self.fd = np.zeros([6,1])                                          # dynamic force
            self.fd[0:3] += self.Arho * self.ddR
            self.fd[0:3] += np.matmul(self.lambd1,np.matmul(base.V2skew(self.dW),self.Rc0)) * self.Arho
            self.fd[0:3] += np.matmul(self.lambd1,np.matmul(base.V2skew(self.W),np.matmul(base.V2skew(self.W),self.Rc0))) * self.Arho
            
            self.fd[3:6] += self.Arho * np.matmul(base.V2skew(self.Rc),self.ddR)
            self.fd[3:6] += np.matmul(self.lambd1,np.matmul(self.Jrho,self.dW))
            self.fd[3:6] += np.matmul(self.lambd1,np.matmul(base.V2skew(self.W),np.matmul(self.Jrho,self.W)))
                            
        
        def iterationUpdate(self,element,dt,beta,gamma):
            
            '''update the configuration quantities'''
            self.dtheta = np.zeros([3,1])                                      # angular displacement incremental vector in iteration
            self.dU = np.zeros([3,1])                                          # linear displacement incremental vector in iteration
            
            for i in range(element.numNodes):
                self.dtheta += self.N[i,0] * element.nodes[i].dtheta
                self.dU += self.N[i,0] * element.nodes[i].dU
            
            self.R = self.R + self.dU
            self.dR = self.dR + self.dU * gamma / (dt * beta)
            self.ddR = self.ddR + self.dU  / (dt**2 * beta)
            
            self.lambd1 = np.matmul(base.theta2Q(self.dtheta),self.lambd1)
            self.theta = base.Q2theta(np.matmul(base.theta2Q(self.dtheta),base.theta2Q(self.theta)))
            Theta = np.matmul(self.lambd0.T,self.theta)
            DTheta = Theta - self.Theta
            self.W = self.W + DTheta * gamma / (dt * beta)
            self.dW = self.dW + DTheta / (dt**2 * beta)
            self.Theta = self.Theta + DTheta
            
            self.Rc = np.matmul(self.lambd1,self.Rc0)
            
            '''update the strain and stress'''
            self.dR_ds = np.zeros([3,1])
            self.dtheta_ds = np.zeros([3,1])
            
            for i in range(element.numNodes):
                self.dR_ds += self.dN_ds[i,0] * element.nodes[i].R[-1]
                self.dtheta_ds += self.dN_ds[i,0] * element.nodes[i].dtheta
            
            self.Strain[0:3] = np.matmul(self.lambd1.T,self.dR_ds) - self.dR0_ds
            
            abs_dtheta = float(sum(self.dtheta ** 2)) ** 0.5
            if abs_dtheta >= element.eps:
                
                self.omega = np.matmul(base.theta2Q(self.dtheta),self.omega)
                self.omega += ((m.sin(abs_dtheta) / abs_dtheta) * self.dtheta_ds
                               + (1 - m.sin(abs_dtheta) / abs_dtheta) * (np.matmul(self.dtheta.T,self.dtheta_ds) / abs_dtheta) * self.dtheta / abs_dtheta
                               + 2.0 * (m.sin(0.5 * abs_dtheta) / abs_dtheta)**2 * np.matmul(base.V2skew(self.dtheta),self.dtheta_ds))
            else:
                self.omega += self.dtheta_ds + 0.5 * np.matmul(base.V2skew(self.dtheta),self.dtheta_ds)
            
            self.Omega = np.matmul(self.lambd1.T,self.omega)
            self.Strain[3:6] = np.array(self.Omega)
            
            self.Stress = np.matmul(self.D,self.Strain)
            self.nu = np.matmul(self.lambd1,self.Stress[0:3])
            self.mu = np.matmul(self.lambd1,self.Stress[3:6])
            
            '''update the applied quantities'''
            self.AM0 = self.AM[0] + np.matmul(self.AM[1],self.R)
            self.fm = np.zeros([6,1])                                          # magnetic force
            self.fg = np.zeros([6,1])                                          # body force
            
            self.fm[0:3] = self.A * np.matmul(self.AM[1].T,np.matmul(self.lambd1,self.RM)) / self.Gm
            self.fm[3:6] = self.A * np.matmul(base.V2skew(np.matmul(self.lambd1,self.RM)),self.AM0) / self.Gm
            self.fg[0:3] = self.A * self.density * self.g
            self.fg[3:6] = self.A * self.density * np.matmul(base.V2skew(self.Rc),self.g)
            
            '''update the dynamic force'''
            self.fd = np.zeros([6,1])                                          
            self.fd[0:3] += self.Arho * self.ddR
            self.fd[0:3] += np.matmul(self.lambd1,np.matmul(base.V2skew(self.dW),self.Rc0)) * self.Arho
            self.fd[0:3] += np.matmul(self.lambd1,np.matmul(base.V2skew(self.W),np.matmul(base.V2skew(self.W),self.Rc0))) * self.Arho
            
            self.fd[3:6] += self.Arho * np.matmul(base.V2skew(self.Rc),self.ddR)
            self.fd[3:6] += np.matmul(self.lambd1,np.matmul(self.Jrho,self.dW))
            self.fd[3:6] += np.matmul(self.lambd1,np.matmul(base.V2skew(self.W),np.matmul(self.Jrho,self.W)))
            
            
            
                
            
            
            

            
            
        
        