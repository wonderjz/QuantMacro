# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 23:05:43 2022

@author: ljz
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

T = 10
alpha = 0.3
beta = 0.9
delta = 0.05
k0 = 1
k = np.ones(T+2)

k[0] = 0
k[1] = 0.75
# set the epsilon
epsilon = 0.001
# set the iteration times
max_times = 1000
times = 1
# learning rate
adjust = 0.0001

#k[T+1] 不能初值为1
k[T+1] = 0

while abs(k0 - k[T+1]) > epsilon and times < max_times:
    if k0 > k[T+1]:
        k[1] = k[1] * (1 + adjust)
    else:
        k[1] = k[1] * (1 - adjust)
    times += 1


    for i in range(T):
        def euler_equation(x):
            A = beta * (alpha * (k[i+1])**(alpha -1) + 1 - delta)
            B = (k[i+1])**alpha + (1 - delta) * k[i+1] - k[i]
            C = x**alpha + (1-delta)*x - k[i+1]
            F = A * C - B
            return F
        #solve the euler equation
        k[i+2] = fsolve(euler_equation, k[i+1]) #给一个初始值k[i+1]

print(k[T+1]) #也就是k0


# correct range
K = np.ones(T+2)
for i in range(len(K)):
    K[i] = k[-1-i]

plt.plot(K)

C = np.ones(T+1)
for i in range(len(k) -1):
    C[i] = K[i] ** alpha + (1-delta) * K[i] - K[i+1]

plt.plot(C)

