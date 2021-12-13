# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 15:09:07 2021

@author: XinLi
"""

import math as m
import numpy as np
import time
import base 
import sys


def mNewtonStatics(inp,maxNewton=80,Iplt = 1):
    processResidual = []                                                       # The residual in the process   
    numNewton = []                                                             # The number of a iteration
    Utotal = np.zeros([inp.numFreedom,1])                                      # The total displacement
    eps = 10e-6
    
    for i in range(inp.numIncrement):
        
        processResidual.append([])
        inp.T.append(inp.T[-1] + inp.dt) 
        j = 0    
        
        while 1:
            for k in inp.elements:
                 
                k.get_K_P(i)

            [inp.K,inp.P] = base.Group_K_P(inp.elements,inp.numFreedom)
        
            for k in inp.Dconstraints:
                k.change_K_P(inp.K,inp.P,j,i)
                 
            for k in inp.Oconstraints:
                k.change_K_P(inp.K,inp.P,j,i)
            
            dU = np.linalg.solve(inp.K,inp.P)
            Utotal += dU
             
            for k in inp.nodes:
                k.updateU(dU)
            
            for k in inp.elements:
                k.update(i)
            
            j += 1
            
            residual = float(sum((dU * inp.coefficient)**2))**0.5/float(sum((Utotal * inp.coefficient)**2))**0.5
            #residual = float(sum((dU**2))**0.5)
            processResidual[i].append(residual)
            
            if (residual <= eps):
                numNewton.append(j)
                break
            
            if(j >= maxNewton):
                print('Newton iteration does not converge',i)
                sys.exit(0)
                break
        for k in inp.nodes:
            k.updateXYZ()
    for i in inp.Dconstraints:
        i.restore_K_P(inp.K,inp.P,dU)
    
    base.Visualization(inp.elements,inp.nodes,inp.T,inp.Name,Iplt)

    return [processResidual,numNewton,Utotal]
                
            
                                                                                                                           
    
    


