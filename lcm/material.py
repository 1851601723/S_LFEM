#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 10:40:07 2021

@author: lixin
"""

class mLEIsotropy():
    '''Isotropic linear elastic material'''
    
    def __init__(self,E,V):
        
        self.E = E                                  #Young modulus
        self.V = V                                  #Poisson's ratio
        self.G = E / (2 * (1 + V))                  #Shear modulus 
        self.lamde = V * E / ((1 + V) * (1 - 2*V) + 0.000001)  #Lame constant
        
        assert self.E > 0,'Young modulus is less than zero'
        assert self.G > 0,'Shear modulus is less than zero'
        
        
class mNeoHookean():
    '''Neo Hookean model material'''
    
    def __init__(self,mu,D):
        
        self.mu = mu                                  #Initial shear modulus
        self.D = D                                    #Incompressible coefficient
        
        assert self.mu > 0, 'Initial shear modulus is less than zero'
        assert self.D > 0 , 'Incompressible coefficient is less than zero'


        