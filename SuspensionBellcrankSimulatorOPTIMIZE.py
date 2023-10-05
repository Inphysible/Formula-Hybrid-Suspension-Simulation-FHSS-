# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 00:30:41 2023

@author: crims
"""

import numpy as np
from SuspensionBellcrankSimulatorMAIN import MAIN as FHSS_MAIN

suspension_parameters = np.loadtxt("SuspensionGeometry.txt", delimiter=",")
suspension_parameters_front = suspension_parameters[:19]
suspension_parameters_rear = suspension_parameters[19:]

for i in range(len(suspension_parameters_front)-3):
    if i == 12 or i == 0:
        continue
    else:
        suspension_parameters_front[i] = suspension_parameters_front[i]/39.3700787
        suspension_parameters_rear[i] = suspension_parameters_rear[i]/39.3700787
p = [14,15]

P = ['camber_angles']

def line(t):
    p0 = np.array([-19.5523, 11.6392, 7.6])
    p1 = np.array([-19.7877, 11.1408, 3.1])
    a = p0 - p1
    x1 = a[0]*t + p1[0]
    x2 = a[1]*t + p1[1]
    x3 = a[2]*t + p1[2]
    
    return np.array([x1, x2, x3])

points = np.array(list(map(line, np.linspace(0, 1, 50))))        

def Optimize_REAR_G(point):
    
    suspension_parameters_rear[p[0]] = point/39.3700787
    suspension_parameters_rear[p[1]] = point/39.3700787
        
    aa, suspension_geometry_dict_rear_new, c, d \
        = FHSS_MAIN(suspension_parameters_front, suspension_parameters_rear)
        
    return suspension_geometry_dict_rear_new[P[0]][:,0]

toe_angles = np.array(list(map(Optimize_REAR_G, points)))
delta_toe = np.zeros(50)

for i, a in enumerate(toe_angles):
    delta_toe[i] = a.max() - a.min()





