# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 18:55:41 2021

@author: XinLi
"""
'''
牛顿-拉普斯迭代非线性问题静态求解
'''

'''记录开始时间'''
import time
time_start = time.time()


'''将Base包放入搜索路径'''
from os import getcwd
import sys

path = getcwd()
path = list(path)
path = path + ['\\','B','a','s','e']
path = ''.join(path)
sys.path.append(path)


'''导入Base包中的支撑数据与支撑函数'''
import BaseFunction as bf
import BaseData as bd


'''进行有限元前处理与单元、载荷等的初始化'''
from ini_3DCR_L2_RobiticArm import *


'''过程中重要参数的初始化'''
Process_residual = []                                       #每次迭代过程的残差
Num_Newton = []                                             #每次迭代过程的迭代次数
U_total = np.zeros([n_freedom,1])                           #位移全量


'''进行迭代求解'''
for i in range(increment_N):
    print(i)
    Process_residual.append([])
    T.append(T[0] + (i+1)*dt)
    K = np.zeros([n_freedom,n_freedom])                     #建立切线刚度矩阵
    PP = np.zeros([n_freedom,1])                            #建立中间载荷矩阵，用于计算的载荷矩阵
    j = 0                                                   #迭代次数
    
    while 1:
        
        for k in element:            
            k.get_K_P()
        [K,P_internal] = bf.Group_K_P(element,n_freedom)    #刚度矩阵与载荷矩阵组集
        P_total = -P_internal + P * Loading[i]              #将外部载荷加入到载荷矩阵
          
        for k in D_constraints:
            k.change_K_P(K,P_total,j,i)                       #引入位移边界条件
            
        for k in Other_constraints:
            k.change_K_P(K,P_total,j,i)                       #引入弹性边界条件
        
        U = np.linalg.solve(K,P_total)
        U_total += U
         
        for k in node:
            k.update_U(U)                                   #节点构型更新
        
        for k in element:
            k.update_element()                              #单元构型更新
        
        j = j + 1
        

        residual = float(sum((U*coefficient)**2))**0.5/float(sum((U_total*coefficient)**2))**0.5 #相对残差
#        residual = float(sum((U**2))**0.5)                                                       #绝对误差
        
        Process_residual[i].append(residual)

        if (residual <= bd.Newton_eps):
            Num_Newton.append(j)
            break
        
        if(j >= bd.Newton_num):
            print('牛顿迭代不收敛')
            sys.exit(0)
            break

    for k in node:
        k.history_xyz.append(np.array(k.spa_xyz+k.U))             #将过程位置全部记录
    for k in element:
        if (k.element_kind == 'element_3DCR_LEL_L2'):
            k.change_T(i+1)
                                  
        
for i in D_constraints:
    i.restore_K_P(K, PP, U)                                   #求位移约束处的支座反力


'''记录结束时间与花费时间'''
time_end = time.time()
time_cost = time_end - time_start