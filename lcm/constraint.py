# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 20:36:07 2021

@author: XinLi
"""
'''
boundary condition 
'''

import numpy as np

class constraint():
    
    def __init__(self,dof,kind,cov,Loading):
        
        self.kind = kind                #kind of constraint
        self.Loading = Loading          #Loading step
        self.cov = cov                  #Constraint constant
        self.dof = dof                  #Constrained degrees of freedom
        self.reaction = 0               #Constraint reaction
        self.K_save = 0                 #Original stiffness
        self.P_save = 0                 #Original load
        
        
    def change_K_P(self,K,P,j,ii):
        
        self.K_save = np.zeros([1,len(P)])
        self.K_save = np.array(K[self.dof])
        self.P_save = float(P[self.dof]) 
        if ii == 0:
            Load = self.Loading[ii]
        else:
            Load = self.Loading[ii] - self.Loading[ii-1]
        
        if(self.kind == 'displacement boundary'):
            
            if j == 0:
                for i in range(len(P)):
                    P[i] = P[i] - K[i,self.dof] * self.cov[0] * Load   
                P[self.dof] = K[self.dof,self.dof] * self.cov[0] * Load 
                                                                 
            else:
                P[self.dof] = 0                                        #The second Newton iteration does not add repeated displacement load
            
            K[self.dof] = np.zeros([1,len(P)])
            K[:,self.dof] = np.zeros(len(P))
            K[self.dof,self.dof] = self.K_save[self.dof]            
            
        elif(self.kind == 'elastic boundary') :
            K[self.dof,self.dof] = K[self.dof,self.dof] + self.cov[0]  #The first elastic boundary constant is the external equivalent spring stiffness coefficient
            P[self.dof] += self.cov[1] * self.Loading[ii]              #The second elastic boundary constant is the external equivalent force
                                                                     
                                                    
            
    def restore_K_P(self,K,P,U):
        K[self.dof] = self.K_save
        if self.kind == 'displacement boundary':
            P[self.dof] = self.P_save
            self.reaction = np.matmul(K[self.dof],U) - self.P_save
