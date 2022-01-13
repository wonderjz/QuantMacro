# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 15:34:08 2022

@author: ljz
"""

import numpy as np

alpha = 0.3
beta = 0.6
delta = 0.3
epsilon = 0.0001

k_min = 0.04
k_max = 0.20
K = np.linspace (k_min, k_max, 1000)

V_last = np.zeros(len(K))
V_new = -10* np.ones(len(K))  #大于 epsilon 即可，防止迭代错误

K_optimal = K.copy()
C_optimal = K.copy()
# 一次迭代
#for i in range(len(K)): # j for each value of k
#    for j in range(len(K)): # j for each value of k'
#        if K[i] ** alpha + (1 - delta) * K[i] - K[j] >= 0:
#            RHS = np.log(K[i] ** alpha + (1 - delta) * K[i] - K[j]) + beta* V_last[j]
#        if V_new[i] < RHS:
#            V_new[i] = RHS
        
while max(abs(V_new- V_last)) > epsilon:  # L1 norm for the vector V 
    V_last = V_new.copy()
          
    for i in range(len(K)): # j for each value of k
        for j in range(len(K)): # j for each value of k'
            if K[i] ** alpha + (1 - delta) * K[i] - K[j] >= 0:
                RHS = np.log(K[i] ** alpha + (1 - delta) * K[i] - K[j]) + beta* V_last[j]
            if V_new[i] < RHS:
                V_new[i] = RHS
                K_optimal[i] = K[j] # 不同i 下 的最优 资本和最优消费
                C_optimal[i] = K[i] ** alpha + (1 - delta) * K[i] - K[j]
                


import matplotlib.pyplot as plt
plt.plot(K,V_new)

plt.plot(K,K_optimal) # 这个只是一期的， 因为上面的迭代实际上是压缩映射去找一期的
#就是在 本期资本为K 下一期为 K_optimal   实际上是policy function

plt.plot(K,C_optimal)









