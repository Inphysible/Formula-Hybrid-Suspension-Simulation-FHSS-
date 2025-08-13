
import SuspensionBellcrankSimulatorCLASSES as C1
import SuspensionBellcrankSimulatorCLASSES3 as C3
import sys
import numpy as np 
import plotly.graph_objects as pgo
import plotly.io as pio
from plotly.subplots import make_subplots
import pandas as pd
import json
import copy
from colorama import init as initColors
from termcolor import colored
from Math_Functions import norm


#Loading the results of the MAIN portion of the simulator into a dictionary for graphing purposes.

initColors()

with open("../outputs/SuspensionDynamicGeometryFRONT.txt") as f:
    suspension_geometry_dict_front = json.load(f)
    suspension_geometry_dict_front_d = copy.deepcopy(suspension_geometry_dict_front)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_front[a] = np.array(suspension_geometry_dict_front[a])
    
with open("../outputs/SuspensionDynamicForcesFRONT.txt") as f:
    suspension_forces_dict_front = json.load(f)
for i, a in enumerate(suspension_forces_dict_front):
    suspension_forces_dict_front[a] = np.array(suspension_forces_dict_front[a])

    
with open("../outputs/SuspensionDynamicGeometryFRONTalt.txt") as f:
    suspension_geometry_dict_front_alt = json.load(f)
    suspension_geometry_dict_front_alt_d = copy.deepcopy(suspension_geometry_dict_front_alt)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_front_alt[a] = np.array(suspension_geometry_dict_front_alt[a])
    
with open("../outputs/SuspensionDynamicForcesFRONTalt.txt") as f:
    suspension_forces_dict_front_alt = json.load(f)
for i, a in enumerate(suspension_forces_dict_front_alt):
    suspension_forces_dict_front_alt[a] = np.array(suspension_forces_dict_front_alt[a])

    
with open("../outputs/SuspensionDynamicGeometryREAR.txt") as f:
    suspension_geometry_dict_rear = json.load(f)
    suspension_geometry_dict_rear_d = copy.deepcopy(suspension_geometry_dict_rear)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_rear[a] = np.array(suspension_geometry_dict_rear[a])
    
with open("../outputs/SuspensionDynamicForcesREAR.txt") as f:
    suspension_forces_dict_rear = json.load(f)
for i, a in enumerate(suspension_forces_dict_rear):
    suspension_forces_dict_rear[a] = np.array(suspension_forces_dict_rear[a])

    
with open("../outputs/SuspensionDynamicGeometryREARalt.txt") as f:
    suspension_geometry_dict_rear_alt = json.load(f)
    suspension_geometry_dict_rear_alt_d = copy.deepcopy(suspension_geometry_dict_rear_alt)
for i, a in enumerate(suspension_geometry_dict_front):
    suspension_geometry_dict_rear_alt[a] = np.array(suspension_geometry_dict_rear_alt[a])
    
with open("../outputs/SuspensionDynamicForcesREARalt.txt") as f:
    suspension_forces_dict_rear_alt = json.load(f)
for i, a in enumerate(suspension_forces_dict_rear_alt):
    suspension_forces_dict_rear_alt[a] = np.array(suspension_forces_dict_rear_alt[a])



suspension_parameters = np.loadtxt("../SuspensionGeometry.txt", delimiter=",")
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

suspension_geometry_dict_front_d['k'] = [suspension_geometry_dict_front_d['k']]
suspension_geometry_dict_front_alt_d['k'] = [suspension_geometry_dict_front_alt_d['k']]
suspension_geometry_dict_rear_d['k'] = [suspension_geometry_dict_rear_d['k']]
suspension_geometry_dict_rear_alt_d['k'] = [suspension_geometry_dict_rear_alt_d['k']]

while len(suspension_geometry_dict_front_d['k']) != len(suspension_geometry_dict_front_d['bot_aarm_arc']):
    suspension_geometry_dict_front_d['k'].append(suspension_geometry_dict_front_d['k'][0])
    suspension_geometry_dict_front_alt_d['k'].append(suspension_geometry_dict_front_alt_d['k'][0])
    suspension_geometry_dict_rear_d['k'].append(suspension_geometry_dict_rear_d['k'][0])
    suspension_geometry_dict_rear_alt_d['k'].append(suspension_geometry_dict_rear_alt_d['k'][0])
    

#This is to convert the geometric data from a dictionary to a Pandas DataFrame for use with Plotly. The dictionaries will still be used
#for dynamics (C3).
susp_data_front = pd.DataFrame.from_dict(suspension_geometry_dict_front_d, dtype=object)
susp_data_rear = pd.DataFrame.from_dict(suspension_geometry_dict_rear_d, dtype=object)
susp_data_front_alt = pd.DataFrame.from_dict(suspension_geometry_dict_front_alt_d, dtype=object)
susp_data_rear_alt = pd.DataFrame.from_dict(suspension_geometry_dict_rear_alt_d, dtype=object)

#Creates Wheel Geometry in each assembly's basis. This is here to avoid having extra code to turn the 3D-arrays into DataFrame elements, and to keep the Wheel class
    #in C1 to keep consistent with that class being designated for constrained components.
FLwheel, FLhub = C1.Wheel(suspension_geometry_dict_front['steering_rod_pivots'], suspension_geometry_dict_front['steering_rod_turning'], suspension_geometry_dict_front['top_aarm_arc'], \
         suspension_geometry_dict_front['bot_aarm_arc'], suspension_geometry_dict_front['upright_vector'], 0.02, 0.254).get_wheel_points(0)
    
FRwheel, FRhub = C1.Wheel(suspension_geometry_dict_front_alt['steering_rod_pivots'], suspension_geometry_dict_front_alt['steering_rod_turning'], suspension_geometry_dict_front_alt['top_aarm_arc'], \
         suspension_geometry_dict_front_alt['bot_aarm_arc'], suspension_geometry_dict_front_alt['upright_vector'], 0.02, 0.254).get_wheel_points(1)
    
RLwheel, RLhub = C1.Wheel(suspension_geometry_dict_rear['steering_rod_pivots'], suspension_geometry_dict_rear['steering_rod_turning'], suspension_geometry_dict_rear['top_aarm_arc'], \
         suspension_geometry_dict_rear['bot_aarm_arc'], suspension_geometry_dict_rear['upright_vector'], 0.02, 0.254).get_wheel_points(0)
    
RRwheel, RRhub = C1.Wheel(suspension_geometry_dict_rear_alt['steering_rod_pivots'], suspension_geometry_dict_rear_alt['steering_rod_turning'], suspension_geometry_dict_rear_alt['top_aarm_arc'], \
         suspension_geometry_dict_rear_alt['bot_aarm_arc'], suspension_geometry_dict_rear_alt['upright_vector'], 0.02, 0.254).get_wheel_points(1)
    
    
suspension_dynamics_dict_front = C3.Dynamics(suspension_geometry_dict_front, suspension_forces_dict_front).ExampleDynamicMotion1()

    ## suspension_dynamics_dict_front_alt = # suspension_dynamics_dict_front; the geometries of the suspension assemblies on either side of the centerline are identical.

suspension_dynamics_dict_rear = C3.Dynamics(suspension_geometry_dict_rear, suspension_forces_dict_rear).ExampleDynamicMotion1()

    ## suspension_dynamics_dict_rear_alt = # suspension_dynamics_dict_rear; the geometries of the suspension assemblies on either side of the centerline are identical.

    
# Create figure and figure objects
colors = ['black', "blue", "purple", "magenta", "orange", "cyan", "gray", "silver", 'black', ['dimgray']]
names = ['Front Left', 'Front Right', 'Rear Left', 'Rear Right']
hovertext = ['Front Left Wheel', 'Front Right Wheel', 'Rear Left Wheel', 'Rear Right Wheel']
data = []
traces = []
frames=[]
sliders = []
frames_exdy1 = []
frames_exdy2 = []
frames_exdy3 = []


fig = pgo.Figure()
for i in range(40):
    traces.append(i)
    if ((i+1)%10) == 0:
        data.append(pgo.Mesh3d(x=[], y=[], z=[],
                      legendgroup=i,
                      facecolor=colors[i],
                      showscale=False))
        colors.append(colors[i])
    else:
        data.append(pgo.Scatter3d(x=[], y=[], z=[],
                      legendgroup=i%10,
                      mode="lines+markers",
                      marker=dict(color=colors[i], size=10)))
        colors.append(colors[i])
fig.add_traces(data = data)

def FrameGen(alldata):
    # Frames
    for j in range(4):
        
        frames.append([pgo.Frame(data = [pgo.Scatter3d(x=np.array([alldata[j][0]['rod_points'][int(i)][0],alldata[j][0]['susp_rod_initial_points'][int(i)][0]]),
                                    y=np.array([alldata[j][0]['rod_points'][int(i)][1],alldata[j][0]['susp_rod_initial_points'][int(i)][1]]),
                                    z=np.array([alldata[j][0]['rod_points'][int(i)][2],alldata[j][0]['susp_rod_initial_points'][int(i)][2]]),\
                                    name=('Pushrod'),
                                    mode = ("lines+markers+text"),
                                    hovertext=('Magnitude: '+str(np.linalg.norm(np.array([np.array(alldata[j][0]['rod_points'][int(i)]) - np.array(alldata[j][0]['susp_rod_initial_points'][int(i)])]))))), 
                                         pgo.Scatter3d(x=np.array([alldata[j][0]['top_aarm_arc'][int(i)][0],alldata[j][0]['bot_aarm_arc'][int(i)][0]]),
                                    y=np.array([alldata[j][0]['top_aarm_arc'][int(i)][1],alldata[j][0]['bot_aarm_arc'][int(i)][1]]),
                                    z=np.array([alldata[j][0]['top_aarm_arc'][int(i)][2],alldata[j][0]['bot_aarm_arc'][int(i)][2]]),\
                                    name=('Kingpin Axis'),
                                    mode = ("lines+markers+text"),
                                    hovertext=('Magnitude: '+str(np.linalg.norm(np.array([np.array(alldata[j][0]['top_aarm_arc'][int(i)]) - np.array(alldata[j][0]['bot_aarm_arc'][int(i)])]))))), 
                                         pgo.Scatter3d(x=np.array([alldata[j][0]['top_aarm_arc'][int(i)][0],alldata[j][3][7][0],alldata[j][0]['top_aarm_arc'][int(i)][0],alldata[j][3][8][0],alldata[j][0]['center_point_top'][int(i)][0],alldata[j][3][7][0]]),
                                    y=np.array([alldata[j][0]['top_aarm_arc'][int(i)][1],alldata[j][3][7][1],alldata[j][0]['top_aarm_arc'][int(i)][1],alldata[j][3][8][1],alldata[j][0]['center_point_top'][int(i)][1],alldata[j][3][7][1]]),
                                    z=np.array([alldata[j][0]['top_aarm_arc'][int(i)][2],alldata[j][3][7][2],alldata[j][0]['top_aarm_arc'][int(i)][2],alldata[j][3][8][2],alldata[j][0]['center_point_top'][int(i)][2],alldata[j][3][7][2]]),\
                                    name=('Top AArm')), 
                                         pgo.Scatter3d(x=np.array([alldata[j][0]['bot_aarm_arc'][int(i)][0],alldata[j][3][10][0],alldata[j][0]['bot_aarm_arc'][int(i)][0],alldata[j][3][11][0],alldata[j][0]['center_point_bot'][int(i)][0],alldata[j][3][10][0]]),
                                    y=np.array([alldata[j][0]['bot_aarm_arc'][int(i)][1],alldata[j][3][10][1],alldata[j][0]['bot_aarm_arc'][int(i)][1],alldata[j][3][11][1],alldata[j][0]['center_point_bot'][int(i)][1],alldata[j][3][10][1]]),
                                    z=np.array([alldata[j][0]['bot_aarm_arc'][int(i)][2],alldata[j][3][10][2],alldata[j][0]['bot_aarm_arc'][int(i)][2],alldata[j][3][11][2],alldata[j][0]['center_point_bot'][int(i)][2],alldata[j][3][10][2]]),\
                                    name=('Bottom AArm'),
                                    surfacecolor=('blue')), 
                                         pgo.Scatter3d(x=np.array([alldata[j][0]['shock_points'][int(i)][0],alldata[j][3][4][0]]),
                                    y=np.array([alldata[j][0]['shock_points'][int(i)][1],alldata[j][3][4][1]]),
                                    z=np.array([alldata[j][0]['shock_points'][int(i)][2],alldata[j][3][4][2]]),\
                                    name=('Shock Absorber'),
                                    hovertext=('Displacement: '+str(alldata[j][0]['total_displacements'][int(i)]))), 
                                         pgo.Scatter3d(x=np.array([alldata[j][0]['steering_rack_turning'][int(i)][0],alldata[j][0]['steering_rack_turning'][0][0],alldata[j][0]['steering_rod_turning'][0][0][0]]),
                                    y=np.array([alldata[j][0]['steering_rack_turning'][0][1],alldata[j][0]['steering_rack_turning'][0][1],alldata[j][0]['steering_rod_turning'][int(i)][0][1]]),
                                    z=np.array([alldata[j][0]['steering_rack_turning'][0][2],alldata[j][0]['steering_rack_turning'][0][2],alldata[j][0]['steering_rod_turning'][int(i)][0][2]]),\
                                    name=('Tie Rod/Steering Rod'),
                                    mode = ("lines+markers+text"),
                                    hovertext=('Magnitude: '+str(np.linalg.norm(np.array([np.array(alldata[j][0]['steering_rack_turning'][0]) - np.array(alldata[j][0]['steering_rod_turning'][int(i)][0])])))+
                                               '\n Steering Rack Extension: '+str(np.linalg.norm(np.array([np.array(alldata[j][0]['steering_rack_turning'][0]) - np.array(alldata[j][0]['steering_rack_turning'][0])]))))),
                                         pgo.Scatter3d(x=np.array([alldata[j][3][1][0],alldata[j][0]['shock_points'][int(i)][0],alldata[j][3][1][0],alldata[j][0]['rod_points'][int(i)][0]]),
                                    y=np.array([alldata[j][3][1][1],alldata[j][0]['shock_points'][int(i)][1],alldata[j][3][1][1],alldata[j][0]['rod_points'][int(i)][1]]),
                                    z=np.array([alldata[j][3][1][2],alldata[j][0]['shock_points'][int(i)][2],alldata[j][3][1][2],alldata[j][0]['rod_points'][int(i)][2]]),\
                                    name=('Bellcrank'),
                                    surfacecolor=('blue')),
                                         pgo.Scatter3d(x=np.array([alldata[j][0]['steering_rod_pivots'][int(i)][0],alldata[j][0]['steering_rod_turning'][int(i)][0][0]]),
                                    y=np.array([alldata[j][0]['steering_rod_pivots'][int(i)][1],alldata[j][0]['steering_rod_turning'][int(i)][0][1]]),
                                    z=np.array([alldata[j][0]['steering_rod_pivots'][int(i)][2],alldata[j][0]['steering_rod_turning'][int(i)][0][2]]),\
                                    name=('Upright'),
                                    hovertext=('Magnitude: '+str(np.linalg.norm(np.array([np.array(alldata[j][0]['steering_rod_turning'][int(i)][0]) - np.array(alldata[j][0]['steering_rod_pivots'][int(i)])]))))),
                                         pgo.Scatter3d(x=[alldata[j][2][i][0]], y=[alldata[j][2][i][1]], z=[alldata[j][2][i][2]],
                                    name=('Wheel Hub'),
                                    hovertext=('Point on Upright Wheel rotates about')),
                                         pgo.Mesh3d(x=alldata[j][1][i][0], y=alldata[j][1][i][1], z=alldata[j][1][i][2],
                                    name=('Wheel'),
                                    hovertext=hovertext[j])], 
                           traces = traces[(j*10):((j+1)*10)],
                           name = f'{names[j]} Suspension Compression Increment {i}',
                           group = f'{names[j]} Suspension Compression Increments'
                          ) for i in range(int(alldata[j][3][-1][0]))])
    return frames

alldata = [[susp_data_front, FLwheel, FLhub, suspension_parameters_front], [susp_data_front_alt, FRwheel, FRhub, suspension_parameters_front_alt], \
                  [susp_data_rear, RLwheel, RLhub, suspension_parameters_rear], [susp_data_rear_alt, RRwheel, RRhub, suspension_parameters_rear_alt]]
frames = FrameGen(alldata)
fig_exdy1 = pgo.Figure(data = fig.data, frames=frames[0]+frames[1]+frames[2]+frames[3])
fig_exdy2 = pgo.Figure(data = fig.data, frames=frames[0]+frames[1]+frames[2]+frames[3])
fig_exdy3 = pgo.Figure(data = fig.data, frames=frames[0]+frames[1]+frames[2]+frames[3])
fig_manual = pgo.Figure(data = fig.data, frames=frames[0]+frames[1]+frames[2]+frames[3])

def frame_args(duration):
    return {
            "frame": {"duration": duration},
            "mode": "immediate",
            "fromcurrent": True, 
            "transition": {"duration": duration, "easing": "quadratic-in-out"},
            }
for i in range(4):
    sliders.append([
        {"transition": {"duration": 11, "easing": "cubic-in-out"},
         "pad": {"b": 10, "t": 10},
         "len": 0.25,
         "x": 0.04+(i*0.25),
         "y": 0,
         "currentvalue": {
            "font": {"size": 14},
            "prefix": str(names[i])+" Suspension Compression Increment ",
            "visible": True,
            "xanchor": "right"
        },
         "steps": [
                     {"args": [[f.name], frame_args(11)],
                      "label": str(j),
                      "method": "animate",
                      } for j, f in enumerate(frames[i])
                  ]
         },
            ])
    
for i in range(len(suspension_dynamics_dict_front['disp_indices'])):
    frames_exdy1.append(frames[0][suspension_dynamics_dict_front['disp_indices'][i]]) 
    frames_exdy1.append(frames[1][suspension_dynamics_dict_front['disp_indices'][i]]) 
    frames_exdy1.append(frames[2][suspension_dynamics_dict_rear['disp_indices'][i]]) 
    frames_exdy1.append(frames[3][suspension_dynamics_dict_rear['disp_indices'][i]])
    
    frames_exdy2.append(frames[0][suspension_dynamics_dict_front['disp_indices'][i]])
    
    frames_exdy3.append(frames[2][suspension_dynamics_dict_front['disp_indices'][i]])

updatemenus_exdy1 = {"buttons": [
    {
        "args": [frames_exdy1, {"frame": {"duration": 1},
                        "fromcurrent": True, "transition": {"duration": 0}}],
        "label": "Play",
        "method": "animate"
    },
    {
        "args": [[None], {"frame": {"duration": 0},
                          "mode": "immediate",
                          "transition": {"duration": 0}}],
        "label": "Pause",
        "method": "animate"
    }
                            ],
"direction": "left",
"pad": {"r": 10, "t": 47},
"showactive": False,
"type": "buttons",
"x": 0.04,
"xanchor": "right",
"y": 0,
"yanchor": "top"
        }
             
updatemenus_exdy2 = {"buttons": [
    {
        "args": [frames_exdy2, {"frame": {"duration": 1},
                        "fromcurrent": True, "transition": {"duration": 0}}],
        "label": "Play",
        "method": "animate"
    },
    {
        "args": [[None], {"frame": {"duration": 0},
                          "mode": "immediate",
                          "transition": {"duration": 0}}],
        "label": "Pause",
        "method": "animate"
    }
                            ],
"direction": "left",
"pad": {"r": 10, "t": 47},
"showactive": False,
"type": "buttons",
"x": 0.04,
"xanchor": "right",
"y": 0,
"yanchor": "top"
        }

updatemenus_exdy3 = {"buttons": [
    {
        "args": [frames_exdy3, {"frame": {"duration": 1},
                        "fromcurrent": True, "transition": {"duration": 0}}],
        "label": "Play",
        "method": "animate"
    },
    {
        "args": [[None], {"frame": {"duration": 0},
                          "mode": "immediate",
                          "transition": {"duration": 0}}],
        "label": "Pause",
        "method": "animate"
    }
                            ],
"direction": "left",
"pad": {"r": 10, "t": 47},
"showactive": False,
"type": "buttons",
"x": 0.04,
"xanchor": "right",
"y": 0,
"yanchor": "top"
        }

scene = {
    'xaxis': { "title": "x"},
    'yaxis': { "title": "y"},
    'zaxis': { "title": "z"},

}

fig_exdy1.update_layout(updatemenus=[updatemenus_exdy1],
         hovermode='closest',
         scene=scene,
         scene_aspectmode='data'
         )
fig_exdy2.update_layout(updatemenus=[updatemenus_exdy2],
         hovermode='closest',
         scene=scene,
         scene_aspectmode='data'
         )
fig_exdy3.update_layout(updatemenus=[updatemenus_exdy3],
         hovermode='closest',
         scene=scene,
         scene_aspectmode='data'
         )

fig_manual.update_layout(sliders=sliders[0]+sliders[1]+sliders[2]+sliders[3],
                         hovermode='closest',
                         scene=scene, 
                         scene_aspectmode='data')


print('----FRONT ASSEMBLY INFORMATION----')  
counter = 0

pushrod = np.round(np.array(list(map(norm, (suspension_geometry_dict_front["susp_rod_initial_points"] - suspension_geometry_dict_front["rod_points"])))),4)
steeringrod_steer = np.round(np.array(list(map(norm, (suspension_geometry_dict_front["steering_rod_turning"][0] - suspension_geometry_dict_front["steering_rack_turning"])))),4)
steeringrod_compress = np.round(np.array(list(map(norm, (suspension_geometry_dict_front["steering_rod_turning"][:,0] - suspension_geometry_dict_front["steering_rack_turning"][0])))),4)
steeringarm_steer = np.round(np.array(list(map(norm, (suspension_geometry_dict_front["steering_rod_turning"][0] - suspension_geometry_dict_front["steering_rod_pivots"][0])))),4)
steeringarm_compress = np.round(np.array(list(map(norm, (suspension_geometry_dict_front["steering_rod_turning"][:,0] - suspension_geometry_dict_front["steering_rod_pivots"])))),4)  

upright_check = np.all(np.round(suspension_geometry_dict_front["upright_vector_mag"],4) == \
                                       np.round(suspension_geometry_dict_front["upright_vector_mag"][0],4))
pushrod_check = np.all(pushrod == pushrod[0])
steering_rod_check = np.all(steeringrod_steer == steeringrod_compress)
steering_arm_check = np.all(steeringarm_steer == steeringarm_compress)
shock_absorber_forces_check = np.all(np.round((\
             (suspension_forces_dict_front["s_t_f_m"])**2 + \
                 (suspension_forces_dict_front["s_n_f_m"])**2)**(1/2),4) == suspension_forces_dict_front["f_m_s"])
pushrod_forces_check = np.all(np.round(((suspension_forces_dict_front["r_i_p_f_m"])**2 + \
             (suspension_forces_dict_front["r_t_f_m"])**2 + \
                 (suspension_forces_dict_front["r_n_f_m"])**2)**(1/2),4) == np.round(np.abs(suspension_forces_dict_front["u_d_f_m"]),4))

print((colored(" Upright Check:", on_color=('on_green' if upright_check else 'on_red'))) + " The magnitude of all the elements in upright vectors should be constant. We can "+ \
          "verify this by checking that the values in upright_vector_mag are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(upright_check)+'\033[0m')

print((colored(" Pushrod Check:", on_color=('on_green' if pushrod_check else 'on_red'))) + " The magnitude of all the elements in pushrod vectors should be constant. We can "+ \
          "verify this by checking that the norms of the vectors defining the pushrod are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(pushrod_check)+'\033[0m')
    
print((colored(" Steering Rod Check:", on_color=('on_green' if steering_rod_check else 'on_red'))) + " The magnitude of all the elements in steering rod vectors should be constant. We can "+ \
          "verify this by checking that the norms of the vectors defining the steering rod while the assembly is compressing and steering "+ \
              "are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(steering_rod_check)+'\033[0m')

print((colored(" Steering Arm Check:", on_color=('on_green' if steering_arm_check else 'on_red'))) + " The magnitude of all the elements in steering arm (moment arm applied to the upright by the steering rod) "+ \
      "vectors should be constant. We can "+ \
          "verify this by checking that the norms of the vectors defining the steering arm while the assembly is compressing and steering "+ \
              "are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(steering_arm_check)+'\033[0m')

print((colored(" Shock Absorber Forces Check:", on_color=('on_green' if shock_absorber_forces_check else 'on_red'))) + ' The Tangent and Normal components of the forces exerted on the bellcrank from normal force '+ \
      'by the shock absorber should sum to form the force vectors defining the shock absorber, the magnitudes of which should be equal (Note that the vectors '+ \
          'will not be tested due to symmetry issues with the components, and that the Binormal components of the force should be ZERO) - '+ \
          'Magnitudes: \033[4m'+str(shock_absorber_forces_check)+'\033[0m')
    
print((colored(" Pushrod Forces Check:", on_color=('on_green' if pushrod_forces_check else 'on_red'))) + ' The Tangent, Normal, and Binormal components of the forces exerted on the bellcrank from normal force '+ \
      'by the suspension rod should sum to form the force vectors defining the suspension rod, the magnitudes of which should be equal (Note that the vectors '+ \
          'will not be tested due to symmetry issues with the components) - '+ \
          'Magnitudes: \033[4m'+str(pushrod_forces_check)+'\033[0m')
    
if not upright_check \
    or not pushrod_check \
    or not steering_rod_check \
    or not steering_arm_check \
    or not shock_absorber_forces_check \
    or not pushrod_forces_check:
        counter += 1
        print("\n\n\033[2;30;41mFailed: One or more constrained components in the front have length conservation issues\033[0;0m")
        
print("\n\nThis is the max relative torsion force the bellcrank will experience during "+ \
                                                  "suspension compression:",np.round(np.max(np.maximum(suspension_forces_dict_front["r_i_p_f_m"], \
                                                                                            suspension_forces_dict_front["s_i_p_f_m"])),2),"N.")
print("This is the max relative torquing force the bellcrank experiences:",
                  np.round(np.max(np.maximum(suspension_forces_dict_front["r_t_f_m"],suspension_forces_dict_front["s_t_f_m"])),2),"N.")
print("This is the max relative translational force the bellcrank experiences:",
                  np.round(np.max(np.maximum(suspension_forces_dict_front["r_n_f_m"],suspension_forces_dict_front["s_n_f_m"])),2),"N.")
print("This is the max relative vertical load on the wheel when fully displaced:",
                  np.round(np.max(np.abs(suspension_forces_dict_front["v_f_m"])),2),"N.")
print("This is the max stroke distance of the shock:",np.max(np.abs(suspension_geometry_dict_front["total_displacements"])),"m.")
print("This is the initial toe angle:",np.round(suspension_geometry_dict_front["toe_angles"][0,0]*57.2957795,2),"degrees.")
print("This is the initial camber angle:",np.round(suspension_geometry_dict_front["camber_angles"][0][0]*57.2957795,2),"degrees.")
print("This is the initial castor angle:",np.round(suspension_geometry_dict_front["castor_angles"][0]*57.2957795,2),"degrees.")

print('\n\n\n----REAR ASSEMBLY INFORMATION----')

pushrod = np.round(np.array(list(map(norm, (suspension_geometry_dict_rear["susp_rod_initial_points"] - suspension_geometry_dict_rear["rod_points"])))),4)
steeringrod_steer = np.round(np.array(list(map(norm, (suspension_geometry_dict_rear["steering_rod_turning"][0] - suspension_geometry_dict_rear["steering_rack_turning"])))),4)
steeringrod_compress = np.round(np.array(list(map(norm, (suspension_geometry_dict_rear["steering_rod_turning"][:,0] - suspension_geometry_dict_rear["steering_rack_turning"][0])))),4)
steeringarm_steer = np.round(np.array(list(map(norm, (suspension_geometry_dict_rear["steering_rod_turning"][0] - suspension_geometry_dict_rear["steering_rod_pivots"][0])))),4)
steeringarm_compress = np.round(np.array(list(map(norm, (suspension_geometry_dict_rear["steering_rod_turning"][:,0] - suspension_geometry_dict_rear["steering_rod_pivots"])))),4)  

upright_check = np.all(np.round(suspension_geometry_dict_rear["upright_vector_mag"],4) == \
                                       np.round(suspension_geometry_dict_rear["upright_vector_mag"][0],4))
pushrod_check = np.all(pushrod == pushrod[0])
steering_rod_check = np.all(steeringrod_steer == steeringrod_compress)
steering_arm_check = np.all(steeringarm_steer == steeringarm_compress)
shock_absorber_forces_check = np.all(np.round(( \
             (suspension_forces_dict_rear["s_t_f_m"])**2 + \
                 (suspension_forces_dict_rear["s_n_f_m"])**2)**(1/2),4) == suspension_forces_dict_rear["f_m_s"])
pushrod_forces_check = np.all(np.round(((suspension_forces_dict_rear["r_i_p_f_m"])**2 + \
             (suspension_forces_dict_rear["r_t_f_m"])**2 + \
                 (suspension_forces_dict_rear["r_n_f_m"])**2)**(1/2),4) == np.round(np.abs(suspension_forces_dict_rear["u_d_f_m"]),4))

print((colored(" Upright Check:", on_color=('on_green' if upright_check else 'on_red'))) + " The magnitude of all the elements in upright vectors should be constant. We can "+ \
          "verify this by checking that the values in upright_vector_mag are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(upright_check)+'\033[0m')

print((colored(" Pushrod Check:", on_color=('on_green' if pushrod_check else 'on_red'))) + " The magnitude of all the elements in pushrod vectors should be constant. We can "+ \
          "verify this by checking that the norms of the vectors defining the pushrod are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(pushrod_check)+'\033[0m')
    
print((colored(" Steering Rod Check:", on_color=('on_green' if steering_rod_check else 'on_red'))) + " The magnitude of all the elements in steering rod vectors should be constant. We can "+ \
          "verify this by checking that the norms of the vectors defining the steering rod while the assembly is compressing and steering "+ \
              "are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(steering_rod_check)+'\033[0m')

print((colored(" Steering Arm Check:", on_color=('on_green' if steering_arm_check else 'on_red'))) + " The magnitude of all the elements in steering arm (moment arm applied to the upright by the steering rod) "+ \
      "vectors should be constant. We can "+ \
          "verify this by checking that the norms of the vectors defining the steering arm while the assembly is compressing and steering "+ \
              "are all equal at a high level, so this should be true: "+ \
                  '\033[4m'+str(steering_arm_check)+'\033[0m')

print((colored(" Shock Absorber Forces Check:", on_color=('on_green' if shock_absorber_forces_check else 'on_red'))) + ' The Tangent and Normal components of the forces exerted on the bellcrank from normal force '+ \
      'by the shock absorber should sum to form the force vectors defining the shock absorber, the magnitudes of which should be equal (Note that the vectors '+ \
          'will not be tested due to symmetry issues with the components, and that the Binormal components of the force should be ZERO) - '+ \
          'Magnitudes: \033[4m'+str(shock_absorber_forces_check)+'\033[0m')
    
print((colored(" Pushrod Forces Check:", on_color=('on_green' if pushrod_forces_check else 'on_red'))) + ' The Tangent, Normal, and Binormal components of the forces exerted on the bellcrank from normal force '+ \
      'by the suspension rod should sum to form the force vectors defining the suspension rod, the magnitudes of which should be equal (Note that the vectors '+ \
          'will not be tested due to symmetry issues with the components) - '+ \
          'Magnitudes: \033[4m'+str(pushrod_forces_check)+'\033[0m')
    
if not upright_check \
    or not pushrod_check \
    or not steering_rod_check \
    or not steering_arm_check \
    or not shock_absorber_forces_check \
    or not pushrod_forces_check:
        counter += 1
        print("\n\n\033[2;30;41mFailed: One or more constrained components in the rear have length conservation issues\033[0;0m")

print("\n\nThis is the max relative torsion force the bellcrank will experience during "+ \
                                                  "suspension compression:",np.round(np.max(np.maximum(suspension_forces_dict_rear["r_i_p_f_m"], \
                                                                                            suspension_forces_dict_rear["s_i_p_f_m"])),2),"N.")
print("This is the max relative torquing force the bellcrank experiences:",
                  np.round(np.max(np.maximum(suspension_forces_dict_rear["r_t_f_m"],suspension_forces_dict_rear["s_t_f_m"])),2),"N.")
print("This is the max relative translational force the bellcrank experiences:",
                  np.round(np.max(np.maximum(suspension_forces_dict_rear["r_n_f_m"],suspension_forces_dict_rear["s_n_f_m"])),2),"N.")
print("This is the max relative vertical load on the wheel when fully displaced:",
                  np.round(np.max(np.abs(suspension_forces_dict_rear["v_f_m"])),2),"N.")
print("This is the max stroke distance of the shock:",np.max(np.abs(suspension_geometry_dict_rear["total_displacements"])),"m.")
print("This is the initial toe angle:",np.round(suspension_geometry_dict_rear["toe_angles"][0,0]*57.2957795,2),"degrees.")
print("This is the initial camber angle:",np.round(suspension_geometry_dict_rear["camber_angles"][0][0]*57.2957795,2),"degrees.")
print("This is the initial castor angle:",np.round(suspension_geometry_dict_rear["castor_angles"][0]*57.2957795,2),"degrees.")

if counter > 0:
    sys.exit()

# fig_exdy1.write_html("Ex. of Animation of Dynamics (All Wheels).html", auto_open=True)
# fig_exdy2.write_html("Ex. of Animation of Dynamics (Front Wheel).html", auto_open=True)
# fig_exdy3.write_html("Ex. of Animation of Dynamics (Rear Wheel).html", auto_open=True)
fig_manual.write_html("../outputs/Animation of Dyanamics Manually (via sliders).html", auto_open=True)


fig = make_subplots(specs=[[{"secondary_y": True}, {"secondary_y": False}]], rows=1, cols=2, subplot_titles=("Shock Absorber Displacements & Chassis Force", "Chassis Position"))
fig.add_trace(pgo.Scatter(x=suspension_dynamics_dict_front['time_evol'], y=suspension_dynamics_dict_front['shock_motion']['z'], name="Shock Displacement", line=dict(color='blue')), secondary_y=False, row=1, col=1)
fig.add_trace(pgo.Scatter(x=suspension_dynamics_dict_front['time_evol'], y=suspension_dynamics_dict_front['chassis_force']['N'], name="Chassis Force", line=dict(color='red')), secondary_y=True, row=1, col=1)
fig.add_trace(pgo.Scatter(x=suspension_dynamics_dict_front['time_evol'], y=suspension_dynamics_dict_front['chassis_motion']['z']+suspension_parameters_front[10][2], name="Chassis Position", \
                          line=dict(color='grey')), row=1, col=2)
fig.update_xaxes(title_text = "Time, t (s)", row=1, col=1)
fig.update_xaxes(title_text = "Time, t (s)", row=1, col=2)
fig.update_yaxes(title_text = "Displacement (m)", secondary_y=False, row=1, col=1)
fig.update_yaxes(title_text = "Force (N)", secondary_y=True, row=1, col=1)
fig.update_yaxes(title_text = "Position, z (m)", row=1, col=2)
fig.update_layout(title_text="Front Suspension Response")
pio.renderers.default='browser'
fig.show()

fig2 = make_subplots(specs=[[{"secondary_y": True}, {"secondary_y": False}]], rows=1, cols=2, subplot_titles=("Shock Absorber Displacements & Chassis Force", "Chassis Position"))
fig2.add_trace(pgo.Scatter(x=suspension_dynamics_dict_rear['time_evol'], y=suspension_dynamics_dict_rear['shock_motion']['z'], name="Shock Displacement", line=dict(color='blue')), secondary_y=False, row=1, col=1)
fig2.add_trace(pgo.Scatter(x=suspension_dynamics_dict_rear['time_evol'], y=suspension_dynamics_dict_rear['chassis_force']['N'], name="Chassis Force", line=dict(color='red')), secondary_y=True, row=1, col=1)
fig2.add_trace(pgo.Scatter(x=suspension_dynamics_dict_rear['time_evol'], y=suspension_dynamics_dict_rear['chassis_motion']['z']+suspension_parameters_rear[10][2], name="Chassis Position", \
                          line=dict(color='grey')), row=1, col=2)
fig2.update_xaxes(title_text = "Time, t (s)", row=1, col=1)
fig2.update_xaxes(title_text = "Time, t (s)", row=1, col=2)
fig2.update_yaxes(title_text = "Displacement (m)", secondary_y=False, row=1, col=1)
fig2.update_yaxes(title_text = "Force (N)", secondary_y=True, row=1, col=1)
fig2.update_yaxes(title_text = "Position, z (m)", row=1, col=2)
fig2.update_layout(title_text="Rear Suspension Response")
pio.renderers.default='browser'
fig2.show()
