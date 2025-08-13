
import numpy as np
import pandas as pd
from scipy.optimize import root
from scipy.signal import lsim, lti
from Math_Functions import matrix_construct, interpolation

class Dynamics():
    
    def __init__(self, suspension_geometry_dict, suspension_forces_dict):
        
        self._suspension_geometry_dict = suspension_geometry_dict
        self._suspension_forces_dict = suspension_forces_dict
        return 
        
    def ExampleDynamicMotion1(self):
        #Can simply write this as a class function that takes geometry data to do the kinetics computations for
        #equilibrium (ride height), with the numpy array -> dataframe script as class initialization.
        
        total_displacements = interpolation(self._suspension_geometry_dict["total_displacements"])
        total_bottom_aarm_arc = interpolation(self._suspension_geometry_dict["bot_aarm_arc"][:,2] - self._suspension_geometry_dict["bot_aarm_arc"][0,2])
        total_vfm = interpolation(self._suspension_forces_dict["v_f_m"])
        compression_ratios = total_displacements/total_bottom_aarm_arc
        compression_ratios[0] = compression_ratios[1] - (compression_ratios[2] - compression_ratios[1])
        
        g = 9.81
        z_chassis_motion_arr = []
        z_chassis_deltavel_arr = []
        z_chassis_F_arr = []
        shock_absorber_motion_arr = []
        shock_absorber_vel_arr = []
        shock_absorber_F_arr = []
        disp_indices = []
        
        #inputs; these are the only changing parameters between dynamics functions
        
        M_wheel = 7.5
        M_chassis = 62.5
        k_tire = 1500000
        c_tire = 250
        k_shock = 61295
        c_shock = 4500
        preload_force = 500
        
        time_evol = np.linspace(0.06,0.88,2300)
        delta_time = np.array(time_evol[1] - time_evol[0])
        
            #Functions of motion for the ground (surface) vertical displacement. Currently set for zero input.
        U_surface_disp = np.zeros_like(time_evol)
        U_surface_vel = np.zeros_like(time_evol)
            #Functions of normal force the ground (surface) provides. Currently set for weight of quadrant of the car in newtons, plus an oscillating, downward force.
        m_oscill = -((19*128.41)*((1/(6 + 5*np.cos(30*time_evol)))-0.091))
        U_surface_normal_force_total = -(((M_wheel+M_chassis)*g)/0.9381)*(np.sqrt(time_evol))# + m_oscill
        U_surface_normal_force_chassis = -(((M_chassis)*g)/0.9381)*(np.sqrt(time_evol)) + preload_force# + m_oscill
                
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
                shock_absorber_motion_val = compression_ratios[1]*(z_chassis_motion_val[-1] - a[0])
                shock_absorber_vel_val = compression_ratios[1]*(xout2[-1,1] - a[1])
                shock_absorber_motion_arr.append(shock_absorber_motion_val)
                shock_absorber_vel_arr.append(shock_absorber_vel_val)
                shock_absorber_F_val = (k_shock*(shock_absorber_motion_val))+ \
                                           (c_shock*(shock_absorber_vel_val))
                z_chassis_F_val = shock_absorber_F_val+ \
                                               (U_surface_normal_force_chassis[i]*(1/total_vfm[0]))
                shock_absorber_F_arr.append(shock_absorber_F_val)                                               
                z_chassis_F_arr.append(z_chassis_F_val)
                a1 = a
                
            elif i == 1:
                find_displacement = lambda index: (np.abs(np.round(total_displacements[int(index)],3))) - (np.abs(np.round((shock_absorber_motion_val),3)))
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
                shock_absorber_motion_val = compression_ratios[j]*((z_chassis_motion_val[-1]) - a[0])
                shock_absorber_vel_val = compression_ratios[j]*(xout2[-1,1] - a[1])
                shock_absorber_motion_arr.append(shock_absorber_motion_val)
                shock_absorber_vel_arr.append(shock_absorber_vel_val)
                shock_absorber_F_val = (k_shock*(shock_absorber_motion_val))+ \
                                           (c_shock*(shock_absorber_vel_val))
                z_chassis_F_val = (k_shock*((shock_absorber_motion_val))+ \
                                           c_shock*((shock_absorber_vel_val))+ \
                                               (U_surface_normal_force_chassis[i]*(1/total_vfm[j])))
                shock_absorber_F_arr.append(shock_absorber_F_val)
                z_chassis_F_arr.append(z_chassis_F_val)
                a1 = a
                disp_indices.append(j)
                
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
                shock_absorber_vel_val = compression_ratios[j]*(xout2[-1,1] - a[1])
                shock_absorber_motion_arr.append(shock_absorber_motion_val)
                shock_absorber_vel_arr.append(shock_absorber_vel_val)
                shock_absorber_F_val = (k_shock*(shock_absorber_motion_val))+ \
                                           (c_shock*(shock_absorber_vel_val))                
                z_chassis_F_val = (k_shock*((shock_absorber_motion_val)) + \
                                           c_shock*((shock_absorber_vel_val)) + \
                                               (U_surface_normal_force_chassis[i]*(1/total_vfm[j])))
                shock_absorber_F_arr.append(shock_absorber_F_val)
                z_chassis_F_arr.append(z_chassis_F_val)
                a1 = a
                disp_indices.append(j)
        
        z_chassis_F = pd.DataFrame(data=z_chassis_F_arr, columns=['N'])
        z_wheel_motion = pd.DataFrame(data=z_wheel_motion_arr, columns=['z'])
        z_chassis_motion = pd.DataFrame(data=z_chassis_motion_arr, columns=['z'])
        shock_absorber_motion = pd.DataFrame(data=shock_absorber_motion_arr, columns=['z'])
        shock_absorber_F = pd.DataFrame(data=shock_absorber_F_arr, columns=['N'])
        dfU_surface_disp = pd.DataFrame(data=U_surface_disp, columns=['z'])
        z_tire_motion = dfU_surface_disp - z_wheel_motion
                
        suspension_dynamics_dict = {"chassis_force": z_chassis_F,
                                    "wheel_motion": z_wheel_motion,
                                    "chassis_motion": z_chassis_motion,
                                    "shock_motion": shock_absorber_motion,
                                    "shock_force": shock_absorber_F,
                                    "surface_Normal_force": U_surface_normal_force_chassis,
                                    "tire_motion": z_tire_motion,
                                    "time_evol": time_evol,
                                    "delta_time": delta_time, 
                                    "disp_indices": np.asarray(disp_indices)}
        
        return suspension_dynamics_dict
    