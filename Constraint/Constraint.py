# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 19:32:03 2023

@author: XinLi
"""

import numpy as np

class Constraint():
    
    def __init__(self,dof,kind,loadList):
        
        self.kind = kind                         # kind of constraint
        self.dof = dof                           # constrained degrees of freedom
        self.loadList = loadList                 # loading list
        self.reaction = [np.zeros(1)]            # constraint reaction
        self.K_save = 0                          # original stiffness
        self.P_save = 0                          # original load
        
    
    def change_K_P(self,K,P,ii,j,nodes,h,beta,gamma):
        
        if self.loadList[ii] == 'pass':
            pass
        else:
            dofID = nodes[self.dof[0]].dofID[self.dof[1]]
            self.K_save = np.zeros([1,len(P)])
            self.K_save = np.array(K[dofID])
            self.P_save = float(P[dofID])
            
            if(self.kind == 'displacement boundary'):
                if self.dof[1] == 0 or self.dof[1] == 1 or self.dof[1] == 2:
                    load = self.loadList[ii] + nodes[self.dof[0]].R[0][self.dof[1],0] - nodes[self.dof[0]].R[-1][self.dof[1],0]
                elif self.dof[1] == 3 or self.dof[1] == 4 or self.dof[1] == 5:
                    load = self.loadList[ii] + nodes[self.dof[0]].axiVector[0][self.dof[1]-3,0] - nodes[self.dof[0]].axiVector[-1][self.dof[1]-3,0]
                
                for i in range(len(P)):
                    P[i] = P[i] - K[i,dofID] *  load   
                
                P[dofID] = K[dofID,dofID] * load 
                K[dofID] = np.zeros([1,len(P)])
                K[:,dofID] = np.zeros(len(P))
                K[dofID,dofID] = self.K_save[dofID]            
            
            elif(self.kind == 'force boundary'):
            
                P[dofID] += self.loadList[ii]                                   
        
            elif(self.kind == 'viscous boundary'):
                if self.dof[1] == 0 or self.dof[1] == 1 or self.dof[1] == 2:
                    V = nodes[self.dof[0]].dR[-1][self.dof[1],0]
                else:
                    V = nodes[self.dof[0]].w[-1][self.dof[1]-3,0]
                a = self.loadList[ii][0]
                b = self.loadList[ii][0]
                c = self.loadList[ii][0]
                P[dofID] -= a * V + b * V * abs(V) + c * V**3
                if self.dof[1] == 0 or self.dof[1] == 1 or self.dof[1] == 2:
                    K[dofID,dofID] += a * gamma / (beta * h)
        
                                                                     
                                                    
            
    def restore_K_P(self,K,P,U,nodes):
        dofID = nodes[self.dof[0]].dofID[self.dof[1]]
        K[dofID] = self.K_save
        if self.kind == 'displacement boundary':
            P[dofID] = self.P_save
            self.reaction.append(np.matmul(K[dofID],U) - self.P_save)

            

