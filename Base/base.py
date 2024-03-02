# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 19:04:01 2023

@author: XinLi
"""

import numpy as np
import math
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
import shutil
import xlsxwriter as xw


'''iterative parameters'''
eps = 10 ** -7            # error tolerance
alpha = 10 ** 10          # penalty function factor
Newton_eps = 10 ** -5     # error tolerance on Newton Laplace iteration
Newton_num = 80           # maximum iterative times of the Newton Laplace iteration



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
        for j in range(len(i.dofID)):
            P[i.dofID[j]] = P[i.dofID[j]] + i.P[j]
            for k in range(len(i.dofID)):
                K[i.dofID[j],i.dofID[k]] = K[i.dofID[j],i.dofID[k]] + i.K[j,k]
    return [K,P]

def Visualization(elements,nodes,Dconstraints,T,Name,Iplt):
    T.insert(0,0.00)
    n = len(nodes)
    m = len(T)
    e = len(elements)
    xyz = []
    track = []
    for i in range(m):
        track.append([])
        xyz.append(np.zeros([3,n]))
        for j in range(n):
            for k in range(3):
                xyz[i][k,j] = nodes[j].R[i][k,0]
        for j in range(e):
            track[i].append(np.zeros([3,len(elements[j].nodes)]))
            for k in range(len(elements[j].nodes)):
                for l in range(3):
                    track[i][j][l,k] = elements[j].nodes[k].R[i][l,0]
    xyz = np.array(xyz)
    resultPath = os.getcwd()
    resultPath = resultPath + '\\Output\\' + Name
    if (os.path.exists(resultPath)):
        shutil.rmtree(resultPath)
    os.mkdir(resultPath)
    result = resultPath + '\\' + 'result.xlsx'
    workBook = xw.Workbook(result)
  
    workSheet1 = workBook.add_worksheet('node radius vector')
    workSheet1.write(0, 0, "Time")
    for i in range(n):
        nodeName = "node"+str(i+1)
        workSheet1.write(0, 3*i+1, nodeName)
    for i in range(m):
        workSheet1.write(i+1, 0, T[i])
        for j in range(n):
            workSheet1.write(i+1, 3*j+1, xyz[i][0,j])
            workSheet1.write(i+1, 3*j+2, xyz[i][1,j])
            workSheet1.write(i+1, 3*j+3, xyz[i][2,j])
    workSheet2 = workBook.add_worksheet('node linear velocity')
    workSheet2.write(0, 0, "Time")
    for i in range(n):
        nodeName = "node"+str(i+1)
        workSheet2.write(0, 3*i+1, nodeName)
    for i in range(m):
        workSheet2.write(i+1, 0, T[i])
        for j in range(n):
            workSheet2.write(i+1, 3*j+1, nodes[j].dR[i][0,0])
            workSheet2.write(i+1, 3*j+2, nodes[j].dR[i][1,0])
            workSheet2.write(i+1, 3*j+3, nodes[j].dR[i][2,0])
    workSheet3 = workBook.add_worksheet('node linear acceleration')
    workSheet3.write(0, 0, "Time")
    for i in range(n):
        nodeName = "node"+str(i+1)
        workSheet3.write(0, 3*i+1, nodeName)
    for i in range(m):
        workSheet3.write(i+1, 0, T[i])
        for j in range(n):
            workSheet3.write(i+1, 3*j+1, nodes[j].ddR[i][0,0])
            workSheet3.write(i+1, 3*j+2, nodes[j].ddR[i][1,0])
            workSheet3.write(i+1, 3*j+3, nodes[j].ddR[i][2,0])
    workSheet4 = workBook.add_worksheet('node axial vector')
    workSheet4.write(0, 0, "Time")
    for i in range(n):
        nodeName = "node"+str(i+1)
        workSheet4.write(0, 3*i+1, nodeName)
    for i in range(m):
        workSheet4.write(i+1, 0, T[i])
        for j in range(n):
            workSheet4.write(i+1, 3*j+1, nodes[j].axiVector[i][0,0])
            workSheet4.write(i+1, 3*j+2, nodes[j].axiVector[i][1,0])
            workSheet4.write(i+1, 3*j+3, nodes[j].axiVector[i][2,0])
    workSheet5 = workBook.add_worksheet('node angular velocity')
    workSheet5.write(0, 0, "Time")
    for i in range(n):
        nodeName = "node"+str(i+1)
        workSheet5.write(0, 3*i+1, nodeName)
    for i in range(m):
        workSheet5.write(i+1, 0, T[i])
        for j in range(n):
            workSheet5.write(i+1, 3*j+1, nodes[j].W[i][0,0])
            workSheet5.write(i+1, 3*j+2, nodes[j].W[i][1,0])
            workSheet5.write(i+1, 3*j+3, nodes[j].W[i][2,0])
    workSheet6 = workBook.add_worksheet('node angular acceleration')
    workSheet6.write(0, 0, "Time")
    for i in range(n):
        nodeName = "node"+str(i+1)
        workSheet6.write(0, 3*i+1, nodeName)
    for i in range(m):
        workSheet6.write(i+1, 0, T[i])
        for j in range(n):
            workSheet6.write(i+1, 3*j+1, nodes[j].dW[i][0,0])
            workSheet6.write(i+1, 3*j+2, nodes[j].dW[i][1,0])
            workSheet6.write(i+1, 3*j+3, nodes[j].dW[i][2,0])
    workSheet7 = workBook.add_worksheet('reactions')
    workSheet7.write(0, 0, "Time")
    for i in range(len(Dconstraints)):
        constraintName = "constraint"+str(i+1)
        workSheet7.write(0, i+1, constraintName)
    for i in range(m):
        workSheet7.write(i+1, 0, T[i])
        for j in range(len(Dconstraints)):
            workSheet7.write(i+1, j+1, Dconstraints[j].reaction[i][0])
            
    
    workBook.close()
    
    
    
    if Iplt == 1:
        os.mkdir(resultPath+'\\3D')
        os.mkdir(resultPath+'\\XY')
        os.mkdir(resultPath+'\\XZ')
        os.mkdir(resultPath+'\\YZ')
        
        
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
                
            for j in range(e):
                ax.plot(track[i][j][0,:],track[i][j][1,:],track[i][j][2,:],c='r')
                
            plt.title('t='+'%.2f'%T[i]+'s',fontsize=20)
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
            for j in range(e):
                plt.plot(track[i][j][0,:],track[i][j][2,:],c='r')                        
                           
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
            for j in range(e):
                plt.plot(track[i][j][0,:],track[i][j][1,:],c='r')                        
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
            for j in range(e):
                plt.plot(track[i][j][1,:],track[i][j][2,:],c='r')                        
            plt.title('t='+'%.2f'%T[i]+'s')        
            plt.savefig(resultPath+'\\YZ\\' + str(i) + '.jpg')
            plt.cla() 
    
            