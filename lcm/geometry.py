#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 20:13:32 2021

@author: lixin
"""

class g3DBeam():
    
    def __init__(self,A,Ix,Iy,k):
        
        self.k = k       #Correction factor
        self.A = A       #area of cross-section
        self.Ix = Ix     #X-direction moment of inertia
        self.Iy = Iy     #Y-direction moment of inertia
        self.J = Ix + Iy #Y-direction moment of inertia
        
        assert self.A > 0,  'the area of cross-section is less than zero'
        assert self.Ix > 0, 'the X-direction moment of inertia is less than zero'
        assert self.Iy > 0, 'the Y-direction moment of inertia is less than zero'
