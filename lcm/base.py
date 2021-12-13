#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 10:20:38 2021

@author: lixin
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import shutil

'''the location and weight of Gaussian integral point'''

guassLocation =   [[0.000000000000000],
                  [-0.577350269189626,  0.577350269189626],
                  [-0.774596669241483,  0.000000000000000,  0.774596669241483],
                  [-0.861136311594053, -0.339981043584856,  0.339981043584856,  0.861136311594053],
                  [-0.906179845938664, -0.538469310105683,  0.000000000000000,  0.538469310105683,  0.906179845938664],
                  [-0.932469514203152, -0.661209386466265, -0.238619186083197,  0.238619186083197,  0.661209386466265,  0.932469514203152]]

guassWeight =   [[2.000000000000000],
                [1.000000000000000, 1.000000000000000],
                [0.555555555555554, 0.888888888888889, 0.555555555555554],
                [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454],
                [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189],
                [0.171324492379171, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379171]]







def quat2Q(q):
    
    '''Quaternion to rotation matrix'''   
    Q = 2 * np.array([[float(q[0]**2 + q[1]**2 - 0.5), float(q[1]*q[2] - q[3]*q[0]), float(q[1]*q[3] + q[2]*q[0])],
                      [float(q[1]*q[2] + q[3]*q[0]), float(q[0]**2 + q[2]**2 - 0.5), float(q[3]*q[2] - q[1]*q[0])],
                      [float(q[1]*q[3] - q[2]*q[0]), float(q[2]*q[3] + q[1]*q[0]), float(q[0]**2 + q[3]**2 - 0.5)]])
    
    return Q



def theta2quat(theta,eps=10**-7):

    '''Axis vector to Quaternion'''
    q = np.zeros((4,1))
    t = float((sum(theta ** 2)) ** 0.5)
    if t <= eps:
        q = np.array([[1],[0],[0],[0]])
    else:
        q[0][0] = float(math.cos(t / 2))
        q[1:4] = math.sin(t / 2) * theta / t
    
    return q

def quat2theta(quat,eps=10**-7):

    '''Quaternion to Axis vector'''  
    theta = np.zeros((3,1))
    a = 2 * math.acos(float(quat[0,0]))
    
    if a <= eps:
        theta = np.array([[0],[0],[0]])
    else:
        theta = quat[1:4] * a / math.sin(a/2)
    return theta

def Q2quat(Q):

    '''rotation matrix to Quaternion'''    
    q = np.zeros([4,1])
    M = max([Q.trace(),Q[0][0],Q[1][1],Q[2][2]])
    
    if M == Q.trace():
        
        q[0][0] = 0.5 * (1 + M) ** 0.5
        for i in range(1,4):
            q[i][0] = 0.25 * (Q[(i+1)%3][i%3] - Q[i%3][(i+1)%3]) / q[0][0]
    
    else: 
        
        for i in range(1,4):
            if M == Q[i-1,i-1]:
                q[i,0] = (0.5 * M + 0.25 * (1 - Q.trace())) ** 0.5
                q[0,0] = 0.25 * (Q[(i+1)%3,i%3] - Q[i%3,(i+1)%3]) / q[i,0]
                for j in range(1,4):
                    if i != j:
                        q[j,0] = 0.25 * (Q[i - 1,j - 1] + Q[j - 1,i - 1]) / q[i,0]
                break
            else:
                
                continue
    t = float(sum(q ** 2) ** 0.5)
    q = q / t
    return q
        
def theta2Q(theta):
    
    Q = quat2Q(theta2quat(theta))
    return Q
    
def Q2theta(Q):
    
    theta = quat2theta(Q2quat(Q))
    return theta
    
def V2skew(V):
    
    '''Axis vector to Antisymmetric matrix''' 
    skew = np.array([[0.0,-float(V[2]),float(V[1])],
                     [float(V[2]),0.0,-float(V[0])],
                     [-float(V[1]),float(V[0]),0.0]])
    
    return skew

def skew2V(skew):
    
    '''Antisymmetric matrix to Axis vector'''   
    V = np.array([[skew[2,1]],
                  [skew[0,2]],
                  [skew[1,0]]])
    return V



    
    


def Group_K_P(elements,numFreedom):
    K = np.zeros([numFreedom,numFreedom])
    P = np.zeros([numFreedom,1])
    for i in elements:
        for j in range(len(i.dofId)):
            P[i.dofId[j]] = P[i.dofId[j]] + i.P[j]
            for k in range(len(i.dofId)):
                K[i.dofId[j],i.dofId[k]] = K[i.dofId[j],i.dofId[k]] + i.K[j,k]
    return [K,P]



def Visualization(elements,nodes,T,Name,Iplt):
    if Iplt == 1:
        resultPath = os.getcwd()
    
        resultPath = resultPath + '\\oup\\' + Name
        if (os.path.exists(resultPath)):
            shutil.rmtree(resultPath)
        os.mkdir(resultPath)
        os.mkdir(resultPath+'\\3D')
        os.mkdir(resultPath+'\\XY')
        os.mkdir(resultPath+'\\XZ')
        os.mkdir(resultPath+'\\YZ')
        
        n = len(nodes)
        m = len(T)
        xyz = []
        for i in range(m):
            xyz.append(np.zeros([3,n]))
            for j in range(n):
                for k in range(3):
                    xyz[i][k,j] = nodes[j].historyXYZ[i][k,0]
        xyz = np.array(xyz)
        Lmax = np.max(xyz)
        Lmin = np.min(xyz)
        D = max(abs(1.1*Lmin - 0.1*Lmax),1.1*Lmax - 0.1*Lmin)
        for i in range(m):
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.set_xlim([-D,D])
            ax.set_ylim([-D,D])
            ax.set_zlim([-D,D])
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')
        
            for j in range(n):
                ax.scatter(xyz[i][0,j],xyz[i][1,j],xyz[i][2,j],c ='r')
                
            ax.plot(xyz[i][0,:],xyz[i][1,:],xyz[i][2,:])
            plt.title('t='+'%.2f'%T[i]+'s')
            plt.savefig(resultPath+'\\3D\\' + str(i) + '.jpg')
            plt.cla()
        
        for i in range(m):
            fig = plt.figure() 
            plt.ylim([-D,D])
            plt.axis('equal')
            plt.xlabel('x')
            plt.ylabel('z')
            for j in range(n):
                plt.scatter(xyz[i][0,j],xyz[i][2,j],c ='r')                        
            plt.plot(xyz[i][0,:],xyz[i][2,:])                
            plt.title('t='+'%.2f'%T[i]+'s')        
            plt.savefig(resultPath+'\\XZ\\' + str(i) + '.jpg')
            plt.cla() 
            
        for i in range(m):
            fig = plt.figure()
            plt.ylim([-D,D])
            plt.axis('equal')
            plt.xlabel('x')
            plt.ylabel('y')
            for j in range(n):
                plt.scatter(xyz[i][0,j],xyz[i][1,j],c ='r')                        
            plt.plot(xyz[i][0,:],xyz[i][1,:])                
            plt.title('t='+'%.2f'%T[i]+'s')        
            plt.savefig(resultPath+'\\XY\\' + str(i) + '.jpg')
            plt.cla() 
            
        for i in range(m):
            fig = plt.figure()
            plt.ylim([-D,D])
            plt.axis('equal')
            plt.xlabel('y')
            plt.ylabel('z')
            for j in range(n):
                plt.scatter(xyz[i][1,j],xyz[i][2,j],c ='r')                        
            plt.plot(xyz[i][1,:],xyz[i][2,:])                
            plt.title('t='+'%.2f'%T[i]+'s')        
            plt.savefig(resultPath+'\\YZ\\' + str(i) + '.jpg')
            plt.cla() 
            
            
def NUBSextraction(U,p):
    m = len(U)
    a = p+1
    b = a+1
    
    PL = []
    ne = 0
    nn = m - a
    nodeList = []
    C = []
    elementList = [U[0]]
    count = 0
    while count < m-1:
        
        if U[count] == U[count+1]:
            count = count + 1
        else:
            elementList.append(U[count+1])
            count = count + 1
            ne = ne + 1
            nodeList.append([])
            
            
     
    for i in range(m-a):
        PL.append([U[i],U[i+a]])
    
    for i in range(ne):
        for j in range(m-a): 
            if (elementList[i]>=PL[j][0]) and (elementList[i+1]<=PL[j][1]):
                nodeList[i].append(j)
                
                
                
    for i in range(ne):
        C.append(np.zeros([a,a]))
    C[0] = np.eye(p+1)
    nb = 1
    while b<m:
        if nb < ne:
            
            C[nb] = np.eye(p+1)
            
        i = b
        while (b<m) and (U[b]==U[b-1]):
            b = b + 1

        
        multi = b - i + 1
        
        if multi <= p:

            numer = U[b-1] - U[a-1]
            alphas = np.zeros([p-multi,1])
            for j in range(p,multi,-1):
                alphas[j-multi-1,0] = numer / (U[a+j-1]-U[a-1])
            r = p - multi
            
            for j in range(1,r+1):
                save = r - j + 1
                s = multi + j
               
                for k in range(p+1,s,-1):
                    
                    alpha = alphas[k-s-1,0]
                    C[nb-1][:,k-1] = alpha * C[nb-1][:,k-1] + (1.0 - alpha) * C[nb-1][:,k-2]
                
                if b<m:
                    C[nb][save-1:save+j,save-1] = np.array(C[nb-1][p-j:p+1,p])
            
            nb = nb + 1
            if b<m:
                a = b
                b = b + 1
                


    return [ne,nn,nodeList,C]
            
                                       
            