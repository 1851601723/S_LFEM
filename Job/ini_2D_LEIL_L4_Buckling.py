# -*- coding: utf-8 -*-
"""
Created on Sat Jun  5 15:43:16 2021

@author: XinLi
"""

from Material.Material_LE_Isotropy import material_LE_Isotropy
from Element.Element_2DS_LEIL_L4 import element_2DS_LEIL_L4
from Node.Node_2D_solid import node_2D_solid
from Constraint.Constraint import constraint

import numpy as np
import math
import os
import shutil

'''
算例名称
'''
Name = 'Job/ini_2D_LEIL_L4_Buckling'
if (os.path.exists(Name)):
    shutil.rmtree(Name)
os.mkdir(Name)
'''
进行计算设置
'''
T_total = 1.0                 #总时间
increment_N = 1              #加载次数
Loading = np.array(range(1,increment_N+1)) / increment_N
dt = T_total / increment_N    #时间步
T = [0.0]                     #当前时间

'''
首先进行网格划分与初始数据输入
'''
L = 10                        #梁长
H = 2                         #梁宽
n = 100                       #横向单元个数
m = 5                         #纵向单元个数
P = -1.0e6

material = []
material.append(material_LE_Isotropy(1.2*(10**11),0.3))



node = []
for i in range(n+1):
    for j in range(m+1):
        node.append(node_2D_solid(np.array([[L * i / n],[H * j / m]]),0))
        node[-1].assign_dof_id((len(node) - 1) * 2)
        
n_freedom = node[-1].dof_num[-1]+1


element = []
for i in range(n):
    for j in range(m):
        element.append(element_2DS_LEIL_L4([node[i*(m+1)+j],
                                            node[(i+1)*(m+1)+j],
                                            node[(i+1)*(m+1)+j+1],
                                            node[i*(m+1)+j+1]], material[0], 0))
    
D_constraints = []
for i in range(2*m+2):
    D_constraints.append(constraint(i,'displacement boundary',[0],Loading))

    
Other_constraints = []
for i in range(m+1):
    Other_constraints.append(constraint(node[-(i+1)].dof_num[0],'elastic boundary',[0,P/(m+1)],Loading))
    
P = np.zeros([n_freedom,1])         #建立载荷矩阵，储存原始载荷
for i in Other_constraints:
    
    P[i.dof] += i.cov[1]

for i in element:
    for j in range(8):
        P[i.dof[j]] += i.BF[j,0]

#####################################################################################
coefficient = [1,1] * len(node)
coefficient = np.array([coefficient]).T 
'''
牛顿迭代收敛残差计算的系数矩阵 不同单元需要使用不同的系数矩阵
'''
######################################################################################

