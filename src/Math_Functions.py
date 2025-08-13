
import numpy as np
from scipy.interpolate import interp1d

def cross(a, b):
    return np.array([a[1]*b[2] - a[2]*b[1],-(a[0]*b[2] - a[2]*b[0]),a[0]*b[1] - a[1]*b[0]])

def norm(a):
    return ((a[0]**2) + (a[1]**2) + (a[2]**2))**(1/2)

def div(a, b):
    return a/b

def multi(a, b):
    return a*b

def matrix_construct(x1, x2, x3):
    return np.array([x1, x2, x3])

def Disk(V1, V2, P1, R, count):
    #V1 & V2 must be orthonormal
    n = V1
    n = n/(norm(V1))
    bi = V2
    bi = bi/(norm(V2))
    theta = np.linspace(0, 2*np.pi, int(count))
    x, y, z =   ((R*(np.cos(theta)))*n[0])+((R*(np.sin(theta)))*bi[0]) + P1[0], \
                ((R*(np.cos(theta)))*n[1])+((R*(np.sin(theta)))*bi[1]) + P1[1], \
                ((R*(np.cos(theta)))*n[2])+((R*(np.sin(theta)))*bi[2]) + P1[2]
    
    return x, y, z

def interpolation(a):
    
    indices = np.linspace(0,len(a)-1,len(a))
    new_indices = np.linspace(0,len(a)-1,500)
    f_a = interp1d(indices, a, kind='quadratic')
    
    return f_a(new_indices)

# def Cylinder(V1, V2, P1, P2, R, h, count):
#     #V1 & V2 must be orthogonal
#     v = P2 - P1
#     v = v/(norm(v))
#     n = V1
#     n = n/(norm(V1))
#     bi = V2
#     bi = bi/(norm(V2))
#     theta = np.linspace(0, 2*np.pi, int(count))
#     vscal = np.linspace(0, h, 2)
#     theta, vscal = np.meshgrid(theta, vscal)
#     x, y, z =   ((R*(np.cos(theta)))*n[0])+((R*(np.sin(theta)))*bi[0])+(vscal*v[0]) + P1[0], \
#                 ((R*(np.cos(theta)))*n[1])+((R*(np.sin(theta)))*bi[1])+(vscal*v[1]) + P1[1], \
#                 ((R*(np.cos(theta)))*n[2])+((R*(np.sin(theta)))*bi[2])+(vscal*v[2]) + P1[2]
    
#     return x, y, z
