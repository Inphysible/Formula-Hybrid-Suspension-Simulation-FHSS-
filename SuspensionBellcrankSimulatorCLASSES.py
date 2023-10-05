
import numpy as np
from scipy import interpolate
from scipy.optimize import root
from Math_Functions import cross, norm, div, multi, Disk#, Cylinder

class circle():

    def __init__(self, Vec_1, n, Center_vec, reference_length, initial_points):

        self._Vec_1 = Vec_1
        self._n = n
        self._Center_vec = Center_vec
        self._reference_length = reference_length
        self._initial_points = initial_points
        self._unit_n = self._n/(np.sqrt(np.dot(self._n,self._n)))
        self._Vec_2 = np.cross(self._unit_n, self._Vec_1)
        return

    def function_circle(self, theta):

        f = self._Vec_1*np.cos(theta) + self._Vec_2*np.sin(theta) + self._Center_vec
        return f

    def residual(self, theta, initial_points):
        b = self.function_circle(theta) - initial_points
        return self._reference_length - norm(b)

    def solve_theta(self, prev_theta):

        theta_vector = np.zeros(len(self._initial_points))
        for i, a in enumerate(self._initial_points):
            f = lambda theta : self.residual(theta, a)
            theta_iterated = root(f, prev_theta[i], method='lm')
            if not theta_iterated.success:
                print("failed")

            theta_vector[i] = theta_iterated.x
        #NORMALIZING THETA_VECTOR VALUES TO BE BETWEEN 0 AND 2pi RADIANS
        for i in range(len(theta_vector)):
            if (2*np.pi)<np.absolute(theta_vector[i]):
                while (2*np.pi)<np.absolute(theta_vector[i]):
                    if theta_vector[i]<1:
                        theta_vector[i] += (2*np.pi)
                    if 1<theta_vector[i]:
                        theta_vector[i] -= (2*np.pi)
                else:
                    break
        return theta_vector
    
class shock(circle):
    
    def __init__(self, Vec_1, n ,Center_vec, reference_length, initial_points, shock_reference_point, shock_reference_length, k):
        self._shock_reference_point = shock_reference_point
        self._shock_reference_length = shock_reference_length
        self._k = k
        super().__init__(Vec_1, n ,Center_vec, reference_length, initial_points)
        
    def get_shock_points(self, theta_vector):
        return np.array(list(map(self.function_circle, theta_vector)))
    
    def get_total_displacements_shock(self, shock_points):
        
        point_to_circle = shock_points - self._shock_reference_point
        lengths = np.asarray(list(map(norm, point_to_circle)))
        
        return np.asarray(lengths - self._shock_reference_length)
    
    def get_relative_force_magnitudes_shock(self, shock_forces, adjusted_shock_points):
         
        point_to_circle = adjusted_shock_points - self._shock_reference_point
        lengths = np.asarray(list(map(norm, point_to_circle)))
        unit_point_to_circle = np.asarray(list(map(div, point_to_circle, lengths)))
    
        return np.asarray(list(map(multi, shock_forces, unit_point_to_circle)))
    
    def shock_tangent_forces(self, f_m_s_v, shock_points):
        
        data = shock_points-self._Center_vec
        distance = np.zeros([len(data)])
        s_t_f = np.zeros([len(data), 3])
        
        for i, a in enumerate(data):
            chord = data[0] - a
            distance[i] = norm(chord)

        data = data.T
        tck,u = interpolate.splprep(data, u=distance, s=0)
        yderv = interpolate.splev(u,tck,der=1)

        tangents = np.zeros([len(yderv[0]), 3])
            
        for i in range(len(tangents)):
            tangents[i] = np.array([yderv[0][i],yderv[1][i],yderv[2][i]])
            s_t_f[i] = tangents[i]*((np.dot(f_m_s_v[i],tangents[i]))/((norm(tangents[i]))**2))
            s_t_f[i][np.isnan(s_t_f[i])] = 0
            
        return s_t_f
    
    def shock_tangent_force_magnitudes(self, s_t_f):
        return np.array(list(map(norm, s_t_f)))
            
    def shock_normal_forces(self, f_m_s_v, shock_points):
        
        s_p = shock_points-self._Center_vec
        s_n_f = np.zeros([len(s_p), 3])
        
        for i, a in enumerate(s_p):
            s_n_f[i] = a*(np.dot(f_m_s_v[i],a)/(np.dot(a,a)))
        
        return s_n_f
    
    def shock_normal_force_magnitudes(self, s_n_f):
        return np.array(list(map(norm, s_n_f)))
    
    def shock_in_plane_forces(self, f_m_s_v):
        
        s_i_p_f = np.zeros_like(f_m_s_v)
        
        for i, a in enumerate(f_m_s_v):
            s_i_p_f[i] = self._unit_n*((np.dot(a,self._unit_n))/((norm(self._unit_n))**2))
        
        return s_i_p_f
    
    def shock_in_plane_force_magnitudes(self, s_i_p_f):
        return np.array(list(map(norm, s_i_p_f)))

class suspension_rod(circle):

    def __init__(self, Vec_1, n, Center_vec, reference_length, initial_points, prev_theta=0.1):
        super().__init__(Vec_1, n ,Center_vec, reference_length, initial_points)
        
        if type(prev_theta) == float:
            if prev_theta == 0.1:
                prev_theta = np.zeros([len(self._initial_points)])
                prev_theta[:] = 0.1
                
        self._theta_vector = self.solve_theta(prev_theta)

    def get_rod_points(self):
        return np.array(list(map(self.function_circle, self._theta_vector))), self._theta_vector
        
    def get_force_magnitudes_arm(self, f_m_s_v, rod_points, shock_points):
        
        rod_radial_vec = rod_points - self._Center_vec
        shock_radial_vec = shock_points - self._Center_vec
        tau_bell_vector = np.asarray(list(map(cross, shock_radial_vec, f_m_s_v)))
        f_m_a_v = -np.asarray(list(map(div, (np.asarray(list(map(cross, rod_radial_vec, tau_bell_vector)))), \
            ((np.asarray(list(map(norm, rod_radial_vec))))**2))))

        
        return np.array(list(map(norm, f_m_a_v)))

    def unit_direction(self, rod_points):
        
        directional_vec = rod_points - self._initial_points
        
        return np.array(list(map(div, directional_vec, (np.array(list(map(norm, directional_vec)))))))

    def out_of_plane_vectors(self, rod_points, u_d_v):

        out_of_plane_vec = np.zeros([len(self._theta_vector), 3])
        rod_radial_vec = rod_points - self._Center_vec
        unit_rod_radial_vec = np.array(list(map(div, rod_radial_vec, np.array(list(map(norm, rod_radial_vec))))))
        
        for i, a in enumerate(u_d_v):
            t = self._Vec_1*(np.dot(a,self._Vec_1)/(np.dot(self._Vec_1,self._Vec_1)))
            n = self._Vec_2*(np.dot(a,self._Vec_2)/(np.dot(self._Vec_2,self._Vec_2)))
            comps = n+t
            out_of_plane_vec[i] = comps/norm(comps)

        return out_of_plane_vec, unit_rod_radial_vec   

    def in_plane_vectors(self, u_d_v):
        
        in_plane_vec = np.zeros([len(self._theta_vector), 3])
        t = np.zeros([3])
        
        for i, a in enumerate(u_d_v):
            t = self._unit_n*((np.dot(a,self._unit_n))/(np.dot(self._unit_n,self._unit_n)))
            in_plane_vec[i] = t/norm(t)
        in_plane_vec[np.isnan(in_plane_vec)] = 0
                
        return in_plane_vec

    def rod_tangent_forces(self, f_m_a, rod_points):
        
        data = rod_points-self._Center_vec
        distance = np.zeros([len(data)])
        r_t_f = np.zeros([len(data), 3])
        r_t_v = np.zeros([len(data), 3])
        
        for i, a in enumerate(data):
            chord = data[0] - a
            distance[i] = norm(chord)

        data = data.T
        tck,u = interpolate.splprep(data, u=distance, s=0)
        yderv = interpolate.splev(u,tck,der=1)
    
        tangents = np.zeros([len(yderv[0]), 3])

        for i in range(len(tangents)):
            tangents[i] = np.array([yderv[0][i],yderv[1][i],yderv[2][i]])
            r_t_f[i] = f_m_a[i]*tangents[i] 
            r_t_v[i] = tangents[i]/norm(tangents[i])
            
        return r_t_f, r_t_v
            
    def rod_tangent_force_magnitudes(self, r_t_f):
        return np.array(list(map(norm, r_t_f)))

    
    def out_of_plane_force_magnitudes(self, r_t_f, o_p_v):
        
        v_mag = np.zeros([len(self._theta_vector)])
        
        for i, a in enumerate(o_p_v):
            v_mag[i] = (np.dot(r_t_f[i],r_t_f[i]))/np.dot(r_t_f[i],a)
        v_mag[np.isnan(v_mag)] = 0
        return v_mag
    
    def out_of_plane_forces(self, o_p_v, o_p_f_m):
        return np.array(list(map(multi, o_p_f_m, o_p_v)))

        
    def rod_normal_forces(self, o_p_f, rod_points):
        
        r_p = rod_points-self._Center_vec
        r_n_f = np.zeros([len(rod_points), 3])
        
        for i, a in enumerate(r_p):
            r_n_f[i] = a*(np.dot(o_p_f[i],a)/(np.dot(a,a)))
       
        return r_n_f
    
    def rod_normal_force_magnitudes(self, r_n_f):
        return np.array(list(map(norm, r_n_f)))
    
    def unit_direction_force_magnitudes(self, u_d_v, o_p_f, o_p_f_m):
        
        w_mag = np.zeros([len(self._theta_vector)])
        
        for i, a in enumerate(u_d_v):
            w_mag[i] = (o_p_f_m[i]**2)/(np.dot(a,o_p_f[i]))
        w_mag[np.isnan(w_mag)] = 0
        return w_mag
    
    def unit_direction_forces(self, u_d_v, u_d_f_m):
        return np.array(list(map(multi, u_d_f_m, u_d_v)))
    
    def in_plane_force_magnitudes(self, u_d_f_m, o_p_f_m):
       
       u_mag = np.zeros([len(self._theta_vector)])
       
       u_mag = ((u_d_f_m**2)-(o_p_f_m**2))**(1/2)
       u_mag[np.isnan(u_mag)] = 0
       return u_mag
    
    def in_plane_forces(self, r_i_p_v, r_i_p_f_m):
        #Note the "in_plane" functions correspond to the binormal vectors
        #in the 2d plane of the bellcrank, which causes torsion.
        return np.array(list(map(multi, r_i_p_f_m, r_i_p_v)))
    
    def get_vertical_displacement_force_magnitudes(self, u_d_f):
        return u_d_f[:,2]
    
class Wheel():
    
    def __init__(self, steering_rod_pivots, steering_rod_turning, top_aarm_arc, bot_aarm_arc, upright_vector, wheel_width, wheel_radius):
        
        self._steering_rod_pivots = np.array(steering_rod_pivots)
        self._steering_rod_turning = np.array(steering_rod_turning)
        self._top_aarm_arc = np.array(top_aarm_arc)
        self._bot_aarm_arc = np.array(bot_aarm_arc)
        self._upright_vector = np.array(upright_vector)
        self._wheel_width = np.array(wheel_width)
        self._wheel_radius = np.array(wheel_radius)
        return 
    
    def get_wheel_points(self, alt):
        
        steering = np.zeros([len(self._top_aarm_arc), 3])
        c = np.linspace(len(self._top_aarm_arc), len(self._top_aarm_arc), len(self._top_aarm_arc))
        h = np.linspace(self._wheel_width, self._wheel_width, len(self._top_aarm_arc))
        
        for i, a in enumerate(self._steering_rod_turning):
            steering[i] = a[0]
        V1 = steering - self._steering_rod_pivots
        
        V2 = self._top_aarm_arc - self._steering_rod_pivots
        A = ((self._upright_vector/2) + self._bot_aarm_arc)
        if alt == 0:
            V3 = np.array(list(map(cross, V2, V1)))
        elif alt == 1:
            V3 = np.array(list(map(cross, V1, V2)))
        V3mag = np.array(list(map(norm, V3)))
        A = ((np.array(list(map(div, V3, V3mag))))*h[0]) + A
        R = np.linspace(self._wheel_radius, self._wheel_radius, len(self._top_aarm_arc))
                
        arr = np.array(list(map(Disk,V1,V2,A,R,c)))
        
        return arr, A
    
class Vars():
    
    def __init__(self, vec_1_shock, n, center, shock_ref_len, ip, shock_ref_point, k,
                                                             vec_1_arm, ref_len, shock_forces):
        self._vec_1_shock = vec_1_shock
        self._vec_1_arm = vec_1_arm
        self._n = n
        self._center = center
        self._shock_ref_len = shock_ref_len
        self._ip = ip
        self._shock_ref_point = shock_ref_point
        self._k = k
        self._shock_ref_point = shock_ref_point
        self._ref_len = ref_len
        self._shock_forces = shock_forces
        return
    
    def Run(self):
        geometry_vars_dict = {"k":self._k,
                     "c_arm":[],
                     "rod_points":[],
                     "shock_points":[],
                     "u_d_v":[],
                     "o_p_v":[],
                     "r_t_v":[],
                     "r_n_v":[],
                     "r_i_p_v":[],
                     "total_displacements":[],
                     }
        
        forces_dict = {"f_m_s_v":[],
                       "f_m_s":self._shock_forces,
                       "f_m_a":[],
                       "o_p_f":[],
                       "u_d_f":[],
                       "r_i_p_f":[],
                       "o_p_f_m":[],
                       "u_d_f_m":[],
                       "r_i_p_f_m":[], 
                       "v_f_m":[],  
                       "r_t_f":[],
                       "r_t_f_m":[], 
                       "r_n_f":[],
                       "r_n_f_m":[],
                       "s_t_f":[],
                       "s_t_f_m":[], 
                       "s_n_f":[],
                       "s_n_f_m":[], 
                       "s_i_p_f":[],
                       "s_i_p_f_m":[]
                       }
    
            
        shock_obj = shock(self._vec_1_shock, self._n, self._center, self._shock_ref_len, self._ip, 
                                                  self._shock_ref_point, self._shock_ref_len, self._k)
        arm_obj = suspension_rod(self._vec_1_arm, self._n, self._center, self._ref_len, self._ip)
    
    
        geometry_vars_dict["c_arm"] = (arm_obj.solve_theta(np.zeros(len(self._ip)))) 
        geometry_vars_dict["rod_points"], theta = (arm_obj.get_rod_points()) #ALWAYS KEEP ENABLED FOR GRAPHING PURPOSES
        geometry_vars_dict["shock_points"] = (shock_obj.get_shock_points(geometry_vars_dict["c_arm"])) #ALWAYS KEEP ENABLED FOR GRAPHING PURPOSES
        geometry_vars_dict["total_displacements"] = (shock_obj.get_total_displacements_shock(geometry_vars_dict["shock_points"]))
        
        geometry_vars_dict["u_d_v"] = (arm_obj.unit_direction(geometry_vars_dict["rod_points"]))
        geometry_vars_dict["o_p_v"], geometry_vars_dict["r_n_v"] = (arm_obj.out_of_plane_vectors(geometry_vars_dict["rod_points"], geometry_vars_dict["u_d_v"]))
        geometry_vars_dict["r_i_p_v"] = (arm_obj.in_plane_vectors(geometry_vars_dict["u_d_v"]))

        forces_dict["f_m_s_v"] = (shock_obj.get_relative_force_magnitudes_shock(forces_dict["f_m_s"], geometry_vars_dict["shock_points"]))
        forces_dict["f_m_a"] = (arm_obj.get_force_magnitudes_arm(forces_dict["f_m_s_v"], geometry_vars_dict["rod_points"], geometry_vars_dict["shock_points"]))
        forces_dict["r_t_f"], geometry_vars_dict["r_t_v"] = (arm_obj.rod_tangent_forces(forces_dict["f_m_a"], geometry_vars_dict["rod_points"]))
        forces_dict["o_p_f_m"] = (arm_obj.out_of_plane_force_magnitudes(forces_dict["r_t_f"], geometry_vars_dict["o_p_v"]))
        forces_dict["o_p_f"] = (arm_obj.out_of_plane_forces(geometry_vars_dict["o_p_v"], forces_dict["o_p_f_m"]))
        forces_dict["u_d_f_m"] = (arm_obj.unit_direction_force_magnitudes(geometry_vars_dict["u_d_v"], forces_dict["o_p_f"], forces_dict["o_p_f_m"]))
        forces_dict["u_d_f"] = (arm_obj.unit_direction_forces(geometry_vars_dict["u_d_v"], forces_dict["u_d_f_m"]))
        forces_dict["r_i_p_f_m"] = (arm_obj.in_plane_force_magnitudes(forces_dict["u_d_f_m"], forces_dict["o_p_f_m"]))
        forces_dict["r_i_p_f"] = (arm_obj.in_plane_forces(geometry_vars_dict["r_i_p_v"], forces_dict["r_i_p_f_m"]))
        forces_dict["v_f_m"] = (arm_obj.get_vertical_displacement_force_magnitudes(forces_dict["u_d_f"]))
        forces_dict["r_t_f_m"] = (arm_obj.rod_tangent_force_magnitudes(forces_dict["r_t_f"])) 
        forces_dict["r_n_f"] = (arm_obj.rod_normal_forces(forces_dict["o_p_f"], geometry_vars_dict["rod_points"]))
        forces_dict["r_n_f_m"] = (arm_obj.rod_normal_force_magnitudes(forces_dict["r_n_f"])) 
        forces_dict["s_t_f"] = (shock_obj.shock_tangent_forces(forces_dict["f_m_s_v"], geometry_vars_dict["shock_points"]))
        forces_dict["s_t_f_m"] = (shock_obj.shock_tangent_force_magnitudes(forces_dict["s_t_f"])) 
        forces_dict["s_n_f"] = (shock_obj.shock_normal_forces(forces_dict["f_m_s_v"], geometry_vars_dict["shock_points"]))
        forces_dict["s_n_f_m"] = (shock_obj.shock_normal_force_magnitudes(forces_dict["s_n_f"]))
        forces_dict["s_i_p_f"] = (shock_obj.shock_in_plane_forces(forces_dict["f_m_s_v"]))
        forces_dict["s_i_p_f_m"] = (shock_obj.shock_in_plane_force_magnitudes(forces_dict["s_i_p_f"]))
        
        return geometry_vars_dict, forces_dict