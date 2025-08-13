
import json
import numpy as np
import matplotlib.pyplot as plt

#Loading the results of the MAIN portion of the simulator into a dictionary for reading purposes.
with open("SuspensionDynamicGeometryFRONT.txt") as f:
    suspension_geometry_dict_front = json.load(f)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_front[a] = np.array(suspension_geometry_dict_front[a])
    
with open("SuspensionDynamicForcesFRONT.txt") as f:
    suspension_forces_dict_front = json.load(f)
for i, a in enumerate(suspension_forces_dict_front):
    suspension_forces_dict_front[a] = np.array(suspension_forces_dict_front[a])

    
with open("SuspensionDynamicGeometryFRONTalt.txt") as f:
    suspension_geometry_dict_front_alt = json.load(f)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_front_alt[a] = np.array(suspension_geometry_dict_front_alt[a])
    
with open("SuspensionDynamicForcesFRONTalt.txt") as f:
    suspension_forces_dict_front_alt = json.load(f)
for i, a in enumerate(suspension_forces_dict_front_alt):
    suspension_forces_dict_front_alt[a] = np.array(suspension_forces_dict_front_alt[a])

    
with open("SuspensionDynamicGeometryREAR.txt") as f:
    suspension_geometry_dict_rear = json.load(f)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_rear[a] = np.array(suspension_geometry_dict_rear[a])
    
with open("SuspensionDynamicForcesREAR.txt") as f:
    suspension_forces_dict_rear = json.load(f)
for i, a in enumerate(suspension_forces_dict_rear):
    suspension_forces_dict_rear[a] = np.array(suspension_forces_dict_rear[a])

    
with open("SuspensionDynamicGeometryREARalt.txt") as f:
    suspension_geometry_dict_rear_alt = json.load(f)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_rear_alt[a] = np.array(suspension_geometry_dict_rear_alt[a])
    
with open("SuspensionDynamicForcesREARalt.txt") as f:
    suspension_forces_dict_rear_alt = json.load(f)
for i, a in enumerate(suspension_forces_dict_rear_alt):
    suspension_forces_dict_rear_alt[a] = np.array(suspension_forces_dict_rear_alt[a])



suspension_parameters = np.loadtxt("SuspensionGeometry.txt", delimiter=",")
suspension_parameters_front = suspension_parameters[:19]
suspension_parameters_rear = suspension_parameters[19:]
suspension_parameters_front_alt = np.copy(suspension_parameters_front)
suspension_parameters_rear_alt = np.copy(suspension_parameters_rear)

for i in range(len(suspension_parameters_front)-3):
    if i == 12 or i == 0:
        continue
    else:
        suspension_parameters_front[i] = suspension_parameters_front[i]/39.3700787
        suspension_parameters_rear[i] = suspension_parameters_rear[i]/39.3700787
        suspension_parameters_front_alt[i] = suspension_parameters_front_alt[i]/39.3700787
        suspension_parameters_rear_alt[i] = suspension_parameters_rear_alt[i]/39.3700787
        
        suspension_parameters_front_alt[i][1] = (-1)*suspension_parameters_front_alt[i][1]
        suspension_parameters_rear_alt[i][1] = (-1)*suspension_parameters_rear_alt[i][1]
        
suspension_parameters_front_alt[15][1] = suspension_parameters_front_alt[15][1]+(2*0.060325)

del suspension_parameters, a, f, i

# radtodeg = 180/(np.pi)
# frontcamber = (suspension_geometry_dict_front['camber_angles'] - (suspension_geometry_dict_front['camber_angles'][0,0] + 0.0070162236))*radtodeg
# fronttoe = (suspension_geometry_dict_front['toe_angles'] - (suspension_geometry_dict_front['toe_angles'][0,0] + 0.00227067336))*radtodeg
# frontsteeringangles_inner = np.zeros([len(fronttoe)])

# for j, a in enumerate(fronttoe):
#         frontsteeringangles_inner[j] = np.abs(a[-1] - a[0])


# fronttoe_alt = (suspension_geometry_dict_front_alt['toe_angles'] - (suspension_geometry_dict_front_alt['toe_angles'][0,0] - 0.00227067336))*radtodeg
# frontsteeringangles_outer = np.zeros([len(fronttoe_alt)])

# for j, a in enumerate(fronttoe_alt):
#         frontsteeringangles_outer[j] = np.abs(a[-1] - a[0])

        
# rearcamber = (suspension_geometry_dict_rear['camber_angles'] - (suspension_geometry_dict_rear['camber_angles'][0,0] + 0.0070162236))*radtodeg
# reartoe = (suspension_geometry_dict_rear['toe_angles'] - (suspension_geometry_dict_rear['toe_angles'][0,0]))*radtodeg


# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_front['total_displacements']*1000, frontcamber[:,0])
# ax.set(xlabel='Heave (mm)', ylabel='Camber Angle (Degrees)',
#        title='Change in Camber In with Heave - Front')
# ax.grid()
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_rear['total_displacements']*1000, rearcamber[:,0])
# ax.set(xlabel='Heave (mm)', ylabel='Camber Angle (Degrees)',
#        title='Change in Camber In with Heave - Rear')
# ax.grid()
# plt.show()


# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_front['total_displacements']*1000, fronttoe[:,0])
# ax.set(xlabel='Heave (mm)', ylabel='Toe Angle (Degrees)',
#        title='Change in Toe In with Heave - Front')
# ax.grid()
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_rear['total_displacements']*1000, reartoe[:,0])
# ax.set(xlabel='Heave (mm)', ylabel='Toe Angle (Degrees)',
#        title='Change in Toe In with Heave - Rear')
# ax.grid()
# plt.show()


# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_rear['total_displacements']*1000, frontsteeringangles_inner)
# ax.set(xlabel='Heave (mm)', ylabel='Steering Angle (Degrees)',
#        title='Max Steering Angle with Heave - Inner Wheel')
# ax.grid()
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_rear['total_displacements']*1000, frontsteeringangles_outer)
# ax.set(xlabel='Heave (mm)', ylabel='Steering Angle (Degrees)',
#        title='Max Steering Angle with Heave - Outer Wheel')
# ax.grid()
# plt.show()


# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_front['total_displacements']*1000, suspension_forces_dict_front['v_f_m'])
# ax.set(xlabel='Heave (mm)', ylabel='Ratio',
#        title='Heave to Shock Displacement Force Ratio - Front')
# ax.grid()
# plt.show()

# fig, ax = plt.subplots()
# ax.plot(-suspension_geometry_dict_rear['total_displacements']*1000, suspension_forces_dict_rear['v_f_m'])
# ax.set(xlabel='Heave (mm)', ylabel=' Ratio',
#        title='Heave to Shock Displacement Force Ratio - Rear')
# ax.grid()
# plt.show()

