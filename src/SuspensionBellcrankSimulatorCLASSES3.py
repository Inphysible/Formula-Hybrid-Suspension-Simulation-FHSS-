
import numpy as np 
import pandas as pd
import json
import plotly.graph_objects as pgo
import plotly.io as pio
from plotly.subplots import make_subplots
from scipy.interpolate import interp1d
from scipy.optimize import root
from scipy.signal import lsim, lti
from scipy.constants import g
from Math_Functions import matrix_construct

class Dynamics():
    
    def __init__(self, suspension_geometry_dict, suspension_forces_dict, M_wheel, M_chassis, k_tire, c_tire, k_shock, c_shock, preload_force, t0, t1):
        
        self._suspension_geometry_dict = suspension_geometry_dict
        self._suspension_forces_dict = suspension_forces_dict   
             
        self.M_wheel = M_wheel #7.5
        self.M_chassis = M_chassis #62.5
        self.k_tire = k_tire #1500000
        self.c_tire = c_tire #250
        self.k_shock = k_shock #61295
        self.c_shock = c_shock #4500
        self.preload_force = preload_force #778.45
        self.time_evol = np.linspace(t0, t1, int((t1-t0) * 1e3)) #t0=0, t1=3.88

        return 
        
    def StateSpaceDynamicMotion(self):            
            
        suspension_parameters = np.loadtxt("../SuspensionGeometry.txt", delimiter=",")
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
        indices = np.linspace(0,int(suspension_parameters[-1][0]-1),int(suspension_parameters[-1][0]), dtype=int)
        # new_indices = np.linspace(0,int(suspension_parameters[-1][0])-1,10000)
        # f_bottom_aarm_arc = interp1d(indices, self._suspension_geometry_dict["bot_aarm_arc"][:,2], kind='quadratic')
        # f_displacements = interp1d(indices, self._suspension_geometry_dict["total_displacements"], kind='quadratic')
        # f_vfm = interp1d(indices, suspension_forces_dict_front["v_f_m"], kind='quadratic')
        total_displacements = self._suspension_geometry_dict["total_displacements"]
        total_bottom_aarm_arc = self._suspension_geometry_dict["bot_aarm_arc"][:,2] - self._suspension_geometry_dict["bot_aarm_arc"][0,2]
        total_vfm = self._suspension_forces_dict["v_f_m"]
        compression_ratios = total_displacements/total_bottom_aarm_arc
        compression_ratios[0] = compression_ratios[1] - (compression_ratios[2] - compression_ratios[1])

        z_chassis_motion_arr = []
        z_chassis_deltavel_arr = []
        shock_absorber_motion_arr = []
        shock_absorber_vel_arr = []
        shock_absorber_F_arr = []


        #inputs; these are the only changing parameters

        
        delta_time = np.array(self.time_evol[1] - self.time_evol[0])

            #Functions of motion for the ground (surface) vertical displacement. Currently set for zero input.
        U_surface_disp = (np.zeros_like(self.time_evol))
        U_surface_vel = (np.zeros_like(self.time_evol))
            #Functions of normal force the ground (surface) provides. Currently set for weight of quadrant of the car in newtons.
        m_oscill = -(19*128.41)*((1/(6 + 5*np.cos(30*self.time_evol)))-0.091)
        U_surface_normal_force_total = -(((-(self.M_wheel+self.M_chassis)*g)*(np.e**(-8*self.time_evol)))+((self.M_wheel+self.M_chassis)*g)) + m_oscill
        U_surface_normal_force_chassis = -((-((self.M_chassis)*g)*(np.e**(-8*self.time_evol))+((self.M_chassis)*g))) + self.preload_force + m_oscill

        U_surface_normal_force_chassis = np.where(U_surface_normal_force_chassis > 0, 0, U_surface_normal_force_chassis)

        #outputs

        u1 = np.array(list(map(matrix_construct, U_surface_disp, U_surface_vel, U_surface_normal_force_total)))

        A1 = np.array([[0,1],[-self.k_tire/(self.M_wheel), -self.c_tire/(self.M_wheel)]]) 
        B1 = np.array([[0,0,0],[self.k_tire/(self.M_wheel), self.c_tire/(self.M_wheel), 1/(self.M_wheel)]])
        C1 = np.array([[1,0],[0,1]])
        D1 = np.array([[0,0,0],[0,0,0]])

        system1 = lti(A1, B1, C1, D1)

        tout, z_wheel, xout1 = lsim(system1, U=u1, T=self.time_evol)

        z_wheel_motion_arr = z_wheel[:,0]
        z_wheel_vel = z_wheel[:,1]

        u2 = np.array(list(map(matrix_construct, z_wheel_motion_arr, z_wheel_vel, U_surface_normal_force_chassis)))

        for i, a in enumerate(u2):
            
            if i == 0:    
                A2 = np.array([[0,1],[-((np.abs(compression_ratios[0])))*(self.k_shock/(self.M_chassis)), -((np.abs(compression_ratios[0])))*(self.c_shock/(self.M_chassis))]])
                B2 = np.array([[0,0,0],[(np.abs(compression_ratios[0]))*(self.k_shock/(self.M_chassis)), ((np.abs(compression_ratios[0])))*(self.c_shock/(self.M_chassis)), (np.abs(1/total_vfm[0]))*(1/(self.M_chassis))]])
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
                shock_absorber_F_val = (self.k_shock*(shock_absorber_motion_val))+ \
                                           (self.c_shock*(shock_absorber_vel_val))
                shock_absorber_F_arr.append(shock_absorber_F_val)
  
                a1 = a
                
                percentDone = "Calculating dynamics... " + str(np.round((i/len(u2) * 100), 0)) + "%"
                print(percentDone)
                
            elif i == 1:
                find_displacement = lambda index: (np.abs(np.round(total_displacements[int(index)],7))) - (np.abs(np.round((shock_absorber_motion_val),7)))
                j = int(root(find_displacement, 0).x)
                
                A2 = np.array([[0,1],[-((np.abs(compression_ratios[j])))*(self.k_shock/(self.M_chassis)), -((np.abs(compression_ratios[j])))*(self.c_shock/(self.M_chassis))]])
                B2 = np.array([[0,0,0],[((np.abs(compression_ratios[j])))*(self.k_shock/(self.M_chassis)), ((np.abs(compression_ratios[j])))*(self.c_shock/(self.M_chassis)), (np.abs(1/total_vfm[j]))*(1/(self.M_chassis))]])
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
                shock_absorber_F_val = (self.k_shock*(shock_absorber_motion_val))+ \
                                           (self.c_shock*(shock_absorber_vel_val))
                shock_absorber_F_arr.append(shock_absorber_F_val)

                a1 = a
                
                percentDone = "Calculating dynamics... " + str(np.round((i/len(u2) * 100), 0)) + "%"
                if (percentDone != "Calculating dynamics... " + str(np.round(((i-1)/len(u2) * 100), 0)) + "%"):
                    print(percentDone) 
                else:
                    continue
                
            elif i > 1:
                find_displacement = lambda index: (np.abs(np.round(total_displacements[int(index)],7))) - (np.abs(np.round((shock_absorber_motion_val),7)))
                j = int(root(find_displacement, 0).x)
                
                A2 = np.array([[0,1],[-((np.abs(compression_ratios[j])))*(self.k_shock/(self.M_chassis)), -((np.abs(compression_ratios[j])))*(self.c_shock/(self.M_chassis))]])
                B2 = np.array([[0,0,0],[((np.abs(compression_ratios[j])))*(self.k_shock/(self.M_chassis)), ((np.abs(compression_ratios[j])))*(self.c_shock/(self.M_chassis)), (np.abs(1/total_vfm[j]))*(1/(self.M_chassis))]])
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
                shock_absorber_F_val = (self.k_shock*(shock_absorber_motion_val))+ \
                                           (self.c_shock*(shock_absorber_vel_val))
                shock_absorber_F_arr.append(shock_absorber_F_val)

                a1 = a
                
                percentDone = "Calculating dynamics... " + str(np.round((i/len(u2) * 100), 0)) + "%"
                if (percentDone != "Calculating dynamics... " + str(np.round(((i-1)/len(u2) * 100), 0)) + "%"):
                    print(percentDone) 
                else:
                    continue

        z_chassis_F_arr = self.M_chassis*(z_chassis_deltavel_arr/delta_time) - ((-((self.M_chassis)*g)*(np.e**(-8*self.time_evol))+((self.M_chassis)*g)))
        z_tire_motion = U_surface_disp - z_wheel_motion_arr

        # fig = make_subplots(rows=1, cols=2, subplot_titles=("Car and Shock Absorber Displacements", "Chassis Position"))
        # fig.add_trace(pgo.Scatter(x=self.time_evol, y=z_wheel_motion['z'], name="Wheel Displacement", line=dict(color='blue')), row=1, col=1)
        # fig.add_trace(pgo.Scatter(x=self.time_evol, y=z_chassis_motion['z'], name="Chassis Displacement", line=dict(color='red')), row=1, col=1)
        # fig.add_trace(pgo.Scatter(x=self.time_evol, y=shock_absorber_motion['z'], name="Shock Displacement", line=dict(color='orange')), row=1, col=1)
        # fig.add_trace(pgo.Scatter(x=self.time_evol, y=dfU_surface_disp['z'], name="Surface Displacement", line=dict(color='gray')), row=1, col=1)
        # fig.add_trace(pgo.Scatter(x=self.time_evol, y=z_tire_motion['z'], name="Tire Displacement", line=dict(color='black')), row=1, col=1)
        # fig.add_trace(pgo.Scatter(x=self.time_evol, y=z_chassis_motion['z']+suspension_parameters_rear[10][2], name="Chassis Position", \
        #                         line=dict(color='red')), row=1, col=2)
        # fig.update_xaxes(title_text = "Time, t (s)", row=1, col=1)
        # fig.update_xaxes(title_text = "Time, t (s)", row=1, col=2)
        # fig.update_yaxes(title_text = "Displacement (m)", tickvals=[0, z_wheel_motion['z'][len(self.time_evol)-1], \
        #                                                             z_chassis_motion['z'][len(self.time_evol)-1], \
        #                                                             shock_absorber_motion['z'][len(self.time_evol)-1],
        #                                                             z_tire_motion['z'][len(self.time_evol)-1]], \
        #                                                             row=1, col=1)
        # fig.update_yaxes(title_text = "Position, z (m)", tickvals=[0, \
        #                                                         z_chassis_motion['z'][0]+suspension_parameters_front[10][2], \
        #                                                             z_chassis_motion['z'][len(self.time_evol)-1]+suspension_parameters_front[10][2]], \
        #                                                                 row=1, col=2)
        # fig.update_layout(title_text="Front Suspension Response")

        # pio.renderers.default='browser'

        # fig.write_html("../outputs/TESTDYNAMICS.html", auto_open=False) 
               
        suspension_dynamics_dict = {"chassis_force": np.asarray(z_chassis_F_arr),
                                            "wheel_motion": np.asarray(z_wheel_motion_arr),
                                            "chassis_motion": np.asarray(z_chassis_motion_arr),
                                            "shock_motion": np.asarray(shock_absorber_motion_arr),
                                            "shock_force": np.asarray(shock_absorber_F_arr),
                                            "surface_Normal_force": U_surface_normal_force_chassis,
                                            "tire_motion": z_tire_motion,
                                            "time_evol": self.time_evol,
                                            "delta_time": delta_time, 
                                            "disp_indices": indices}
                
        return suspension_dynamics_dict
    