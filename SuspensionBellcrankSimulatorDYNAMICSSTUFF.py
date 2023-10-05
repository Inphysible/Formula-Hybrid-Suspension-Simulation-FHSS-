
import numpy as np 
import pandas as pd
import json
import plotly.graph_objects as pgo
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.interpolate import interp1d
from scipy.optimize import root
from scipy.signal import lsim, lti
from Math_Functions import matrix_construct

#Loading the results of the MAIN portion of the simulator into a dictionary for graphing purposes.
with open("SuspensionDynamicGeometryFRONT.txt") as f:
    suspension_geometry_dict_front = json.load(f)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_front[str(a)] = np.array(suspension_geometry_dict_front[str(a)])
with open("SuspensionDynamicForcesFRONT.txt") as f:
    suspension_forces_dict_front = json.load(f)
for i, a in enumerate(suspension_forces_dict_front):
    suspension_forces_dict_front[str(a)] = np.array(suspension_forces_dict_front[str(a)])
    
with open("SuspensionDynamicGeometryFRONTalt.txt") as f:
    suspension_geometry_dict_front_alt = json.load(f)
for i, a in enumerate(suspension_geometry_dict_front_alt):
    suspension_geometry_dict_front_alt[str(a)] = np.array(suspension_geometry_dict_front_alt[str(a)])
with open("SuspensionDynamicForcesFRONTalt.txt") as f:
    suspension_forces_dict_front_alt = json.load(f)
for i, a in enumerate(suspension_forces_dict_front_alt):
    suspension_forces_dict_front_alt[str(a)] = np.array(suspension_forces_dict_front_alt[str(a)])
    
with open("SuspensionDynamicGeometryREAR.txt") as f:
    suspension_geometry_dict_rear = json.load(f)
for i, a in enumerate(suspension_geometry_dict_rear):
    suspension_geometry_dict_rear[str(a)] = np.array(suspension_geometry_dict_rear[str(a)])
with open("SuspensionDynamicForcesREAR.txt") as f:
    suspension_forces_dict_rear = json.load(f)
for i, a in enumerate(suspension_forces_dict_rear):
    suspension_forces_dict_rear[str(a)] = np.array(suspension_forces_dict_rear[str(a)])
    
with open("SuspensionDynamicGeometryREARalt.txt") as f:
    suspension_geometry_dict_rear_alt = json.load(f)
for i, a in enumerate(suspension_geometry_dict_rear_alt):
    suspension_geometry_dict_rear_alt[str(a)] = np.array(suspension_geometry_dict_rear_alt[str(a)])
with open("SuspensionDynamicForcesREARalt.txt") as f:
    suspension_forces_dict_rear_alt = json.load(f)
for i, a in enumerate(suspension_forces_dict_rear_alt):
    suspension_forces_dict_rear_alt[str(a)] = np.array(suspension_forces_dict_rear_alt[str(a)])
    
    
suspension_parameters = np.loadtxt("SuspensionGeometry.txt", delimiter=",")
suspension_parameters_front = suspension_parameters[:19]
suspension_parameters_rear = suspension_parameters[19:]
suspension_parameters_front_alt = np.copy(suspension_parameters_front)
suspension_parameters_rear_alt = np.copy(suspension_parameters_rear)

for i in range(len(suspension_parameters_front)-3):
    if i == 12 or i == 0:
        continue
    else:
        suspension_parameters_front_alt[i][1] = (-1)*suspension_parameters_front_alt[i][1]
        suspension_parameters_rear_alt[i][1] = (-1)*suspension_parameters_rear_alt[i][1]
        
suspension_parameters_front_alt[15][1] = suspension_parameters_front_alt[15][1]+(2*0.060325)
        
#Can simply write this as a class function that takes geometry data to do the kinetics computations for
#equilibrium (ride height), with the numpy array -> dataframe script as class initialization.
indices = np.linspace(0,int(suspension_parameters[-1][0])-1,int(suspension_parameters[-1][0]))
new_indices = np.linspace(0,int(suspension_parameters[-1][0])-1,10000)
f_bottom_aarm_arc = interp1d(indices, suspension_geometry_dict_front["bot_aarm_arc"][:,2], kind='quadratic')
f_displacements = interp1d(indices, suspension_geometry_dict_front["total_displacements"], kind='quadratic')
f_vfm = interp1d(indices, suspension_forces_dict_front["v_f_m"], kind='quadratic')
total_displacements = suspension_geometry_dict_front["total_displacements"]
total_bottom_aarm_arc = suspension_geometry_dict_front["bot_aarm_arc"][:,2] - suspension_geometry_dict_front["bot_aarm_arc"][0,2]
total_vfm = suspension_forces_dict_front["v_f_m"]
compression_ratios = total_displacements/total_bottom_aarm_arc
compression_ratios[0] = compression_ratios[1] - (compression_ratios[2] - compression_ratios[1])

g = 9.81
z_chassis_motion_arr = []
z_chassis_deltavel_arr = []
shock_absorber_motion_arr = []
shock_absorber_vel_arr = []

#inputs; these are the only changing parameters

M_wheel = 7.5
M_chassis = 62.5
k_tire = 1500000
c_tire = 250
k_shock = 61295
c_shock = 4500
preload_force = 778.45

time_evol = np.linspace(0,3.88,23000)
delta_time = np.array(time_evol[1] - time_evol[0])

    #Functions of motion for the ground (surface) vertical displacement. Currently set for zero input.
U_surface_disp = (np.zeros_like(time_evol))
U_surface_vel = (np.zeros_like(time_evol))
    #Functions of normal force the ground (surface) provides. Currently set for weight of quadrant of the car in newtons.
m_oscill = -(19*128.41)*((1/(6 + 5*np.cos(30*time_evol)))-0.091)
U_surface_normal_force_total = -(((-(M_wheel+M_chassis)*g)*(np.e**(-8*time_evol)))+((M_wheel+M_chassis)*g)) + m_oscill
U_surface_normal_force_chassis = -((-((M_chassis)*g)*(np.e**(-8*time_evol))+((M_chassis)*g))) + preload_force + m_oscill

U_surface_normal_force_chassis = np.where(U_surface_normal_force_chassis > 0, 0, U_surface_normal_force_chassis)

#outputs

u1 = np.array(list(map(matrix_construct, U_surface_disp, U_surface_vel, U_surface_normal_force_total)))

A1 = np.array([[0,1],[-k_tire/(M_wheel), -c_tire/(M_wheel)]]) 
B1 = np.array([[0,0,0],[k_tire/(M_wheel), c_tire/(M_wheel), 1/(M_wheel)]])
C1 = np.array([[1,0],[0,1]])
D1 = np.array([[0,0,0],[0,0,0]])

system1 = lti(A1, B1, C1, D1)

tout, z_wheel, xout1 = lsim(system1, U=u1, T=time_evol)

z_wheel_motion_arr = z_wheel[:,0]
z_wheel_vel = z_wheel[:,1]

u2 = np.array(list(map(matrix_construct, z_wheel_motion_arr, z_wheel_vel, U_surface_normal_force_chassis)))

for i, a in enumerate(u2):
    if i == 0:    
        A2 = np.array([[0,1],[-((np.abs(compression_ratios[0])))*(k_shock/(M_chassis)), -((np.abs(compression_ratios[0])))*(c_shock/(M_chassis))]])
        B2 = np.array([[0,0,0],[(np.abs(compression_ratios[0]))*(k_shock/(M_chassis)), ((np.abs(compression_ratios[0])))*(c_shock/(M_chassis)), (np.abs(1/total_vfm[0]))*(1/(M_chassis))]])
        C2 = np.array([[1,0]])
        D2 = np.array([[0,0,0]])
        
        system2 = lti(A2, B2, C2, D2)
        
        a0 = np.array(list(map(matrix_construct, np.linspace(0,a[0],100), np.linspace(0,a[1],100), np.linspace(0,a[2],100))))
        t0 = np.linspace(0,delta_time,100)
        
        tout, z_chassis_motion_val, xout2 = lsim(system2, U=a0, T=t0)
        
        z_chassis_motion_arr.append(z_chassis_motion_val[-1])
        z_chassis_deltavel_arr.append(xout2[-1,1] - xout2[0,1])
        shock_absorber_motion_val = -compression_ratios[1]*(z_chassis_motion_val[-1] - a[0])
        shock_absorber_vel_val = -compression_ratios[1]*(xout2[-1,1] - a[1])
        shock_absorber_motion_arr.append(shock_absorber_motion_val)
        shock_absorber_vel_arr.append(shock_absorber_vel_val)
        a1 = a
        
    elif i == 1:
        find_displacement = lambda index: (np.abs(np.round(total_displacements[int(index)],7))) - (np.abs(np.round((shock_absorber_motion_val),7)))
        j = int(root(find_displacement, 0).x)
        
        A2 = np.array([[0,1],[-((np.abs(compression_ratios[j])))*(k_shock/(M_chassis)), -((np.abs(compression_ratios[j])))*(c_shock/(M_chassis))]])
        B2 = np.array([[0,0,0],[((np.abs(compression_ratios[j])))*(k_shock/(M_chassis)), ((np.abs(compression_ratios[j])))*(c_shock/(M_chassis)), (np.abs(1/total_vfm[j]))*(1/(M_chassis))]])
        C2 = np.array([[1,0]])
        D2 = np.array([[0,0,0]])
    
        system2 = lti(A2, B2, C2, D2)
        
        a0 = np.array(list(map(matrix_construct, np.linspace(a1[0],a[0],100),np.linspace(a1[1],a[1],100),np.linspace(a1[2],a[2],100))))
        
        tout, z_chassis_motion_val, xout2 = lsim(system2, U=a0, T=t0, X0=xout2[-1])
        
        z_chassis_motion_arr.append(z_chassis_motion_val[-1])
        z_chassis_deltavel_arr.append(xout2[-1,1] - xout2[0,1])
        shock_absorber_motion_val = -compression_ratios[j]*((z_chassis_motion_val[-1]) - a[0])
        shock_absorber_vel_val = -compression_ratios[1]*(xout2[-1,1] - a[1])
        shock_absorber_motion_arr.append(shock_absorber_motion_val)
        shock_absorber_vel_arr.append(shock_absorber_vel_val)
        a1 = a
        
    elif i > 1:
        find_displacement = lambda index: (np.abs(np.round(total_displacements[int(index)],7))) - (np.abs(np.round((shock_absorber_motion_val),7)))
        j = int(root(find_displacement, 0).x)
        
        A2 = np.array([[0,1],[-((np.abs(compression_ratios[j])))*(k_shock/(M_chassis)), -((np.abs(compression_ratios[j])))*(c_shock/(M_chassis))]])
        B2 = np.array([[0,0,0],[((np.abs(compression_ratios[j])))*(k_shock/(M_chassis)), ((np.abs(compression_ratios[j])))*(c_shock/(M_chassis)), (np.abs(1/total_vfm[j]))*(1/(M_chassis))]])
        C2 = np.array([[1,0]])
        D2 = np.array([[0,0,0]])
    
        system2 = lti(A2, B2, C2, D2)
        
        a0 = np.array(list(map(matrix_construct, np.linspace(a1[0],a[0],100),np.linspace(a1[1],a[1],100),np.linspace(a1[2],a[2],100))))
        
        tout, z_chassis_motion_val, xout2 = lsim(system2, U=a0, T=t0, X0=xout2[-1])
        
        z_chassis_motion_arr.append(z_chassis_motion_val[-1])
        z_chassis_deltavel_arr.append(xout2[-1,1] - xout2[0,1])
        shock_absorber_motion_val = compression_ratios[j]*(z_chassis_motion_val[-1] - a[0])
        shock_absorber_vel_val = -compression_ratios[1]*(xout2[-1,1] - a[1])
        shock_absorber_motion_arr.append(shock_absorber_motion_val)
        shock_absorber_vel_arr.append(shock_absorber_vel_val)
        a1 = a
        print(j)

z_chassis_F_arr = M_chassis*(z_chassis_deltavel_arr/delta_time) - ((-((M_chassis)*g)*(np.e**(-8*time_evol))+((M_chassis)*g)))
z_wheel_motion = pd.DataFrame(data=z_wheel_motion_arr, columns=['z'])
z_chassis_motion = pd.DataFrame(data=z_chassis_motion_arr, columns=['z'])
shock_absorber_motion = pd.DataFrame(data=shock_absorber_motion_arr, columns=['z'])
dfU_surface_disp = pd.DataFrame(data=U_surface_disp, columns=['z'])
z_tire_motion = dfU_surface_disp - z_wheel_motion

fig = make_subplots(rows=1, cols=2, subplot_titles=("Car and Shock Absorber Displacements", "Chassis Position"))
fig.add_trace(pgo.Scatter(x=time_evol, y=z_wheel_motion['z'], name="Wheel Displacement", line=dict(color='blue')), row=1, col=1)
fig.add_trace(pgo.Scatter(x=time_evol, y=z_chassis_motion['z'], name="Chassis Displacement", line=dict(color='red')), row=1, col=1)
fig.add_trace(pgo.Scatter(x=time_evol, y=shock_absorber_motion['z'], name="Shock Displacement", line=dict(color='orange')), row=1, col=1)
fig.add_trace(pgo.Scatter(x=time_evol, y=dfU_surface_disp['z'], name="Surface Displacement", line=dict(color='gray')), row=1, col=1)
fig.add_trace(pgo.Scatter(x=time_evol, y=z_tire_motion['z'], name="Tire Displacement", line=dict(color='black')), row=1, col=1)
fig.add_trace(pgo.Scatter(x=time_evol, y=z_chassis_motion['z']+suspension_parameters_rear[10][2], name="Chassis Position", \
                          line=dict(color='red')), row=1, col=2)
fig.update_xaxes(title_text = "Time, t (s)", row=1, col=1)
fig.update_xaxes(title_text = "Time, t (s)", row=1, col=2)
fig.update_yaxes(title_text = "Displacement (m)", tickvals=[0, z_wheel_motion['z'][len(time_evol)-1], \
                                                            z_chassis_motion['z'][len(time_evol)-1], \
                                                            shock_absorber_motion['z'][len(time_evol)-1],
                                                            z_tire_motion['z'][len(time_evol)-1]], \
                                                             row=1, col=1)
fig.update_yaxes(title_text = "Position, z (m)", tickvals=[0, \
                                                           z_chassis_motion['z'][0]+suspension_parameters_front[10][2], \
                                                               z_chassis_motion['z'][len(time_evol)-1]+suspension_parameters_front[10][2]], \
                                                                 row=1, col=2)
fig.update_layout(title_text="Front Suspension Response")

pio.renderers.default='browser'

fig.write_html("TESTDYNAMICS.html", auto_open=True)