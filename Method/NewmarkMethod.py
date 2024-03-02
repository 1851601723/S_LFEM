# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 17:05:26 2023

@author: XinLi
"""


import numpy as np
import base
import math as m
import time
import sys

def NewmarkMethod(inp,maxNewton=160,Iplt = 1):
    processResidual = []                                                       # The residual in the process   
    numNewton = []                                                             
    eps = 10e-6                                                                # the telorance
    beta = 0.25
    gamma = 0.5
    Utotal = np.zeros([inp.numFreedom,1])
    for i in range(inp.numTimeStep):
        
        processResidual.append([])
        j = 0
        
        for node in inp.nodes:
            node.timeStepUpdate(i,inp.dt[i],beta,gamma)
        for element in inp.elements:
            element.timeStepUpdate(inp.dt[i],beta,gamma)
        
        while 1:
             
            for element in inp.elements:
                element.get_K_P(inp.dt[i],beta,gamma)

            [inp.K,inp.P] = base.Group_K_P(inp.elements,inp.numFreedom)
            
            for k in inp.Dconstraints:
                k.change_K_P(inp.K,inp.P,i,j,inp.nodes,inp.dt[i],beta,gamma)
                 
            for k in inp.Fconstraints:
                k.change_K_P(inp.K,inp.P,i,j,inp.nodes,inp.dt[i],beta,gamma)
            
            dU = np.linalg.solve(inp.K,inp.P)
            
            Utotal += dU
            
            for node in inp.nodes:
                node.iterationUpdate(dU,inp.dt[i],beta,gamma)
            for element in inp.elements:
                element.iterationUpdate(inp.dt[i],beta,gamma)
            
            j += 1
            residual = float(sum((dU * inp.coefficient)**2))**0.5/float(sum((Utotal * inp.coefficient)**2))**0.5
            #residual = float(sum((dU**2))**0.5)
            processResidual[i].append(residual)
            
            if (residual <= eps):
                numNewton.append(j)
                break
            
            if(j >= maxNewton):
                print('Newton iteration does not converge',i,j)
                sys.exit(0)
                break
        print(i,j)
        for j in inp.Dconstraints:
            j.restore_K_P(inp.K,inp.P,dU,inp.nodes)
    
    base.Visualization(inp.elements,inp.nodes,inp.Dconstraints,inp.T,inp.Name,Iplt)
    return [processResidual,numNewton,Utotal]
            
        
                 