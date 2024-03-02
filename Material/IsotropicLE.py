# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 20:41:21 2023

@author: XinLi
"""

from sys import exit

class IsotropicLE():
    
    "isotropic linear elastic constant"
    def __init__(self,E,V,d):
        
        self.E = E                                  # Young's modulus
        self.V = V                                  # Poisson's ratio
        self.G = E / (2 * (1 + V))                  # shear modulus
        self.lame = V * E / ((1 + V) * (1 - 2*V))   # lame constant
        self.density = d                            # density
 
        
        if((E<0.0) or (V>0.5) or (V<-1.0) or (d<0.0)):
            print('isotropic linear elastic constant is incorrect')
            exit(0)
