# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 20:20:46 2023

@author: XinLi
"""

from sys import exit

class Rod3D():
    
    def __init__(self,A,I1,I2,k):
        

        self.k = k           # reasonable factor
        self.A = A           # area of the cross-section
        self.I1 = I1         # moment of inertia of the cross-section at 1-dirction
        self.I2 = I2         # moment of inertia of the cross-section at 2-dirction
        self.J = I1 + I2     # polar moment of inertia of the cross-section
        
        
        if((A<0) or (I1<0) or (I2<0)):
            print('The geometric parameters of the cross-section are incorrect')
            exit(0)
