# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 18:49:22 2021

@author: XinLi
"""

from Material.Material_LE_Isotropy import material_LE_Isotropy
from Element.Element_3DS_LEIL_L8 import element_3DS_LEIL_L8
from Node.Node_3D_solid import node_3D_solid
from Constraint.Constraint import constraint

import numpy as np
import math
import os
import shutil

'''
算例名称
'''

Name = 'Job/3D_L8_CompressBeam'
if (os.path.exists(Name)):
    shutil.rmtree(Name)
os.mkdir(Name)
'''
进行计算设置
'''
T_total = 1.0                 #总时间
increment_N =1              #加载次数
Loading = np.array(range(1,increment_N+1)) / increment_N
dt = T_total / increment_N    #时间步
T = [0.0]                     #当前时间
'''
首先进行网格划分与初始数据输入
'''
L = 10                        #试件长
H = 1                         #试件宽
B = 2                         #试件高
n = 100                       #长方向单元个数
m = 5                         #宽方向单元个数
l = 4                         #高方向单元个数
P = -1.0e6                    #载荷大小

material = []
material.append(material_LE_Isotropy(1.2e11,0.3))



node = []
boundary = [3,[[1,2]]]
ID = np.zeros([n+1,m+1,l+1])
for i in range(n+1):
    for j in range(m+1):        
        for k in range(l+1):
            node.append(node_3D_solid(np.array([[L * i / n],[H * j / m],[B * k / l]]),0))
            node[-1].assign_dof_id((len(node) - 1) * 3)
            ID[i,j,k] = int(len(node) - 1)
        boundary[1].append(node[int(ID[i,j,0])])
        
n_freedom = node[-1].dof_num[-1]+1






element = []
for i in range(n):
    for j in range(m):
        for k in range(l):
            element.append(element_3DS_LEIL_L8([node[int(ID[i,j,k])],
                                                node[int(ID[i+1,j,k])],
                                                node[int(ID[i+1,j+1,k])],
                                                node[int(ID[i,j+1,k])],
                                                node[int(ID[i,j,k+1])],
                                                node[int(ID[i+1,j,k+1])],
                                                node[int(ID[i+1,j+1,k+1])],
                                                node[int(ID[i,j+1,k+1])]], material[0], 0))
Other_constraints = []    
D_constraints = []
for i in range(m+1):
    for j in range(l+1):
        D_constraints.append(constraint(int(ID[0,i,j]) * 3,'displacement boundary',[0],Loading))
        D_constraints.append(constraint(int(ID[0,i,j]) * 3 + 1,'displacement boundary',[0],Loading))
        D_constraints.append(constraint(int(ID[0,i,j]) * 3 + 2,'displacement boundary',[0],Loading))
        Other_constraints.append(constraint(int(ID[-1,i,j]) * 3,'elastic boundary',[0,P/((l+1)*(m+1))],Loading))
Other_constraints.append(constraint(int(ID[-1,-1,-1]) * 3+1,'elastic boundary',[0,0],Loading))


    
P = np.zeros([n_freedom,1])         #建立载荷矩阵，储存原始载荷
for i in Other_constraints:
    
    P[i.dof] += i.cov[1]

for i in element:
    for j in range(8):
        P[i.dof[j]] += i.BF[j,0]

#####################################################################################
coefficient = [1,1,1] * len(node)
coefficient = np.array([coefficient]).T 
'''
牛顿迭代收敛残差计算的系数矩阵 不同单元需要使用不同的系数矩阵
'''
######################################################################################