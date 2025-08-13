
import SuspensionBellcrankSimulatorCLASSES as C1
import SuspensionBellcrankSimulatorCLASSES2 as C2
import numpy as np 
import json
        
def MAIN(suspension_parameters_front, suspension_parameters_rear, alt=0):
    vec_1_arm = (suspension_parameters_front[3]-suspension_parameters_front[1], suspension_parameters_rear[3]-suspension_parameters_rear[1])  
    
    vec_1_shock = (suspension_parameters_front[5]-suspension_parameters_front[1], suspension_parameters_rear[5]-suspension_parameters_rear[1])           
    
    n_b = ((C1.cross(vec_1_arm[0],vec_1_shock[0])), (C1.cross(vec_1_arm[1],vec_1_shock[1])))
    
    center = (suspension_parameters_front[1], suspension_parameters_rear[1])
    
    vec_ref_len = (suspension_parameters_front[3]-suspension_parameters_front[2], suspension_parameters_rear[3]-suspension_parameters_rear[2])
    pushrod_ref_len = (np.sqrt(np.dot(vec_ref_len[0],vec_ref_len[0])), np.sqrt(np.dot(vec_ref_len[1],vec_ref_len[1])))
    
    shock_ref_point = (suspension_parameters_front[4], suspension_parameters_rear[4])
    
    shock_point1 = (suspension_parameters_front[5], suspension_parameters_rear[5])
    shock_ref_len = (np.sqrt(np.dot(shock_point1[0]-shock_ref_point[0],shock_point1[0]-shock_ref_point[0])), \
        np.sqrt(np.dot(shock_point1[1]-shock_ref_point[1],shock_point1[1]-shock_ref_point[1])))
    
    k = (suspension_parameters_front[17][0], suspension_parameters_rear[17][0])
    
    upright_point_top = (suspension_parameters_front[6], suspension_parameters_rear[6])
    
    upright_point_bot = (suspension_parameters_front[9], suspension_parameters_rear[9])
    if alt == 1:
        n_aarm_top = (suspension_parameters_front[8]-suspension_parameters_front[7], suspension_parameters_rear[8]-suspension_parameters_rear[7])
        
        n_aarm_bot = (suspension_parameters_front[11]-suspension_parameters_front[10], suspension_parameters_rear[11]-suspension_parameters_rear[10])
    else:
        n_aarm_top = (suspension_parameters_front[7]-suspension_parameters_front[8], suspension_parameters_rear[7]-suspension_parameters_rear[8])
        
        n_aarm_bot = (suspension_parameters_front[10]-suspension_parameters_front[11], suspension_parameters_rear[10]-suspension_parameters_rear[11])
    
    Vec_1_top_alt1 = (suspension_parameters_front[7]-suspension_parameters_front[6], suspension_parameters_rear[7]-suspension_parameters_rear[6])
    unit_n_aarm_top = (n_aarm_top[0]/(np.sqrt(np.dot(n_aarm_top[0],n_aarm_top[0]))), n_aarm_top[1]/(np.sqrt(np.dot(n_aarm_top[1],n_aarm_top[1]))))
    Vec_1_top_alt2 = (unit_n_aarm_top[0]*(np.dot(unit_n_aarm_top[0],Vec_1_top_alt1[0])), unit_n_aarm_top[1]*(np.dot(unit_n_aarm_top[1],Vec_1_top_alt1[1])))
    center_point_top = (suspension_parameters_front[7] - Vec_1_top_alt2[0], suspension_parameters_rear[7] - Vec_1_top_alt2[1])
    
    Vec_1_bot_alt1 = (suspension_parameters_front[10]-suspension_parameters_front[9], suspension_parameters_rear[10]-suspension_parameters_rear[9])
    unit_n_aarm_bot = (n_aarm_bot[0]/(np.sqrt(np.dot(n_aarm_bot[0],n_aarm_bot[0]))), n_aarm_bot[1]/(np.sqrt(np.dot(n_aarm_bot[1],n_aarm_bot[1]))))
    Vec_1_bot_alt2 = (unit_n_aarm_bot[0]*(np.dot(unit_n_aarm_bot[0],Vec_1_bot_alt1[0])), unit_n_aarm_bot[1]*(np.dot(unit_n_aarm_bot[1],Vec_1_bot_alt1[1]))) 
    center_point_bot = (suspension_parameters_front[10] - Vec_1_bot_alt2[0], suspension_parameters_rear[10] - Vec_1_bot_alt2[1])
    
    which_arm = (suspension_parameters_front[12][0], suspension_parameters_rear[12][0])
    
    steering_rod_point = (suspension_parameters_front[13], suspension_parameters_rear[13])
    
    steering_rack_point1 = (suspension_parameters_front[14], suspension_parameters_rear[14])
    
    vec = (steering_rod_point[0]-steering_rack_point1[0], steering_rod_point[1]-steering_rack_point1[1])
    steering_rod_length = (np.sqrt(np.dot(vec[0],vec[0])), np.sqrt(np.dot(vec[1],vec[1])))
    
    steering_rack_point_final = (suspension_parameters_front[15], suspension_parameters_rear[15])
    
    centerline_point = (suspension_parameters_front[16], suspension_parameters_rear[16])
    
    shock_forces_front = np.ones(int(suspension_parameters_rear[18][0]))
    shock_forces_rear = np.ones(int(suspension_parameters_rear[18][0]))
    
    
    suspension_geometry_dict1 = C2.Vars(upright_point_top[0], center_point_top[0], upright_point_bot[0], center_point_bot[0], \
                                                        n_aarm_top[0], n_aarm_bot[0], suspension_parameters_front[2], int(suspension_parameters_front[18][0]), which_arm[0], \
                                                            steering_rod_length[0], steering_rod_point[0], steering_rack_point1[0], \
                                                                steering_rack_point_final[0], centerline_point[0]).Run()                            
    suspension_geometry_dict2, suspension_forces_dict2 = C1.Vars(vec_1_shock[0], n_b[0], center[0], shock_ref_len[0], suspension_geometry_dict1["susp_rod_initial_points"], shock_ref_point[0], k[0], \
                    vec_1_arm[0], pushrod_ref_len[0], shock_forces_front).Run()
        
    suspension_geometry_dict3 = C2.Vars(upright_point_top[1], center_point_top[1], upright_point_bot[1], center_point_bot[1], \
                                                        n_aarm_top[1], n_aarm_bot[1], suspension_parameters_rear[2], int(suspension_parameters_rear[18][0]), which_arm[1], \
                                                            steering_rod_length[1], steering_rod_point[1], steering_rack_point1[1], \
                                                                steering_rack_point_final[1], centerline_point[1]).Run()                            
    suspension_geometry_dict4, suspension_forces_dict4 = C1.Vars(vec_1_shock[1], n_b[1], center[1], shock_ref_len[1], suspension_geometry_dict3["susp_rod_initial_points"], shock_ref_point[1], k[1], \
                    vec_1_arm[1], pushrod_ref_len[1], shock_forces_rear).Run()
    suspension_geometry_dict_front = {**suspension_geometry_dict1, **suspension_geometry_dict2}
    suspension_geometry_dict_rear = {**suspension_geometry_dict3, **suspension_geometry_dict4}

    return suspension_geometry_dict_front, suspension_geometry_dict_rear, suspension_forces_dict2, suspension_forces_dict4

if __name__ == '__main__':
    
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
            
    suspension_parameters_front_alt[15][1] = suspension_parameters_front_alt[15][1]+(2*np.abs(suspension_parameters_front_alt[15][1] - suspension_parameters_front_alt[14][1]))
    
    suspension_geometry_dict_front, suspension_geometry_dict_rear, suspension_forces_dict2, suspension_forces_dict4 \
        = MAIN(suspension_parameters_front, suspension_parameters_rear)
    suspension_geometry_dict_front_alt, suspension_geometry_dict_rear_alt, suspension_forces_dict2_alt, suspension_forces_dict4_alt \
        = MAIN(suspension_parameters_front_alt, suspension_parameters_rear_alt, 1)
        
    
    for i, a in enumerate(suspension_geometry_dict_front):
        suspension_geometry_dict_front[str(a)] = suspension_geometry_dict_front[str(a)].tolist()
        suspension_geometry_dict_rear[str(a)] = suspension_geometry_dict_rear[str(a)].tolist()
        suspension_geometry_dict_front_alt[str(a)] = suspension_geometry_dict_front_alt[str(a)].tolist()
        suspension_geometry_dict_rear_alt[str(a)] = suspension_geometry_dict_rear_alt[str(a)].tolist()
    for i, a in enumerate(suspension_forces_dict2):
        suspension_forces_dict2[str(a)] = suspension_forces_dict2[str(a)].tolist()
        suspension_forces_dict4[str(a)] = suspension_forces_dict4[str(a)].tolist()
        suspension_forces_dict2_alt[str(a)] = suspension_forces_dict2_alt[str(a)].tolist()
        suspension_forces_dict4_alt[str(a)] = suspension_forces_dict4_alt[str(a)].tolist()
    with open("../outputs/SuspensionDynamicGeometryFRONT.txt", 'w') as g1, open("../outputs/SuspensionDynamicGeometryREAR.txt", 'w') as g2, \
        open("../outputs/SuspensionDynamicGeometryFRONTalt.txt", 'w') as g3, open("../outputs/SuspensionDynamicGeometryREARalt.txt", 'w') as g4:
        g1.truncate(0)
        g1.write(json.dumps(suspension_geometry_dict_front))
        g2.truncate(0)
        g2.write(json.dumps(suspension_geometry_dict_rear))
        g3.truncate(0)
        g3.write(json.dumps(suspension_geometry_dict_front_alt))
        g4.truncate(0)
        g4.write(json.dumps(suspension_geometry_dict_rear_alt))
    g1.close()
    g2.close()
    g3.close()
    g4.close()
    with open("../outputs/SuspensionDynamicForcesFRONT.txt", 'w') as f1, open("../outputs/SuspensionDynamicForcesREAR.txt", 'w') as f2, \
        open("../outputs/SuspensionDynamicForcesFRONTalt.txt", 'w') as f3, open("../outputs/SuspensionDynamicForcesREARalt.txt", 'w') as f4:
        f1.truncate(0)
        f1.write(json.dumps(suspension_forces_dict2))
        f2.truncate(0)
        f2.write(json.dumps(suspension_forces_dict4))
        f3.truncate(0)
        f3.write(json.dumps(suspension_forces_dict2_alt))
        f4.truncate(0)
        f4.write(json.dumps(suspension_forces_dict4_alt))
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    for i, a in enumerate(suspension_geometry_dict_front):
        suspension_geometry_dict_front[str(a)] = np.array(suspension_geometry_dict_front[str(a)])
        suspension_geometry_dict_rear[str(a)] = np.array(suspension_geometry_dict_rear[str(a)])
        suspension_geometry_dict_front_alt[str(a)] = np.array(suspension_geometry_dict_front_alt[str(a)])
        suspension_geometry_dict_rear_alt[str(a)] = np.array(suspension_geometry_dict_rear_alt[str(a)])
    
    for i, a in enumerate(suspension_forces_dict2):
        suspension_forces_dict2[str(a)] = np.array(suspension_forces_dict2[str(a)])
        suspension_forces_dict4[str(a)] = np.array(suspension_forces_dict4[str(a)])
        suspension_forces_dict2_alt[str(a)] = np.array(suspension_forces_dict2_alt[str(a)])
        suspension_forces_dict4_alt[str(a)] = np.array(suspension_forces_dict4_alt[str(a)])
        
    import SuspensionBellcrankSimulatorPLOTS
        