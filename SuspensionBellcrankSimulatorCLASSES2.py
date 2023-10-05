
import SuspensionBellcrankSimulatorCLASSES as C1
import numpy as np
from scipy.optimize import root
from Math_Functions import cross, norm, div

class A_arms():
    
    def __init__(self, upright_point_top, upright_point_bot, n_top, n_bot, center_point_top, center_point_bot, susp_rod_ip1, points_count, which_arm):
        
        self._Vec_top = upright_point_top - center_point_top
        self._Vec_bot = upright_point_bot - center_point_bot
        self._upright_point_top = upright_point_top 
        self._upright_point_bot = upright_point_bot
        self._n_top = n_top
        self._n_bot = n_bot
        self._center_point_top = center_point_top
        self._center_point_bot = center_point_bot
        self._points_count = points_count
        if int(which_arm) == 0:
            self._Vec_susp_rod_alt_top = susp_rod_ip1 - self._center_point_top
            unit_n_top = self._n_top/(np.sqrt(np.dot(self._n_top,self._n_top)))
            self._Vec_susp_rod_alt2_top = unit_n_top*(np.dot(unit_n_top,self._Vec_susp_rod_alt_top))
            self._center_point_susp_rod = self._center_point_top + self._Vec_susp_rod_alt2_top
            self._Vec_susp_rod = susp_rod_ip1 - self._center_point_susp_rod
            self._n_arm = self._n_top
        elif int(which_arm) == 1:
            self._Vec_susp_rod_alt_bot = susp_rod_ip1 - self._center_point_bot
            unit_n_bot = self._n_bot/(np.sqrt(np.dot(self._n_bot,self._n_bot)))
            self._Vec_susp_rod_alt2_bot = unit_n_bot*(np.dot(unit_n_bot,self._Vec_susp_rod_alt_bot))
            self._center_point_susp_rod = self._center_point_bot + self._Vec_susp_rod_alt2_bot
            self._Vec_susp_rod = susp_rod_ip1 - self._center_point_susp_rod
            self._n_arm = self._n_bot
        return

    def top_A_arm_arc(self):
        
        theta = np.linspace(0,2*np.pi,10000)
        vecs = []
        vecs_solved = np.zeros([self._points_count,3])
        aarm_circle = C1.circle(self._Vec_top, self._n_top, self._center_point_top, 0, np.zeros([1,3]))
        
        for i, a in enumerate(theta):
            vecs.append(aarm_circle.function_circle(a))
            theta_sweep = a
            if (vecs[i][2] - vecs[0][2]) >= 0.0520:
                break
        theta_vector = np.linspace(0,theta_sweep,self._points_count)
        
        vecs_solved = np.array(list(map(aarm_circle.function_circle, theta_vector)))  
        
        return vecs_solved, theta_vector
           
    def bot_A_arm_arc(self, top_A_arm_arc):
        
        l_vec = self._upright_point_top - self._upright_point_bot
        upright_length = norm(l_vec)
        
        return C1.suspension_rod(self._Vec_bot, self._n_bot, self._center_point_bot, upright_length, top_A_arm_arc).get_rod_points()
    
    def upright_vectors(self, top_arc, bot_arc):
        
        upright_vectors = top_arc-bot_arc
        upright_lengths = np.array(list(map(norm, upright_vectors)))
        unit_upright_vectors = np.array(list(map(div, upright_vectors, upright_lengths)))
            
        return upright_vectors, unit_upright_vectors
    #Note: the function below describes the arc the suspension rod takes during compression, which serves as the initial
    #points (ip) for that class.
    def suspension_rod_initial_points(self, A_arm_arc_theta):
        
        vecs = np.zeros([self._points_count,3])
        susprod_circle = C1.circle(self._Vec_susp_rod, self._n_arm, self._center_point_susp_rod, 0, np.zeros([1,3]))
        vecs = np.array(list(map(susprod_circle.function_circle, A_arm_arc_theta)))
        
        return vecs
        
class Steering_rod():

    def __init__(self, upright_vectors, unit_upright_vectors, top_A_arm_arc, steering_rod_length, steering_rod_point, steering_rack_point1, steering_rack_point_final):
        self._upright_vector = upright_vectors
        self._unit_upright_vector = unit_upright_vectors
        self._top_A_arm_arc = top_A_arm_arc
        self._steering_rod_length = steering_rod_length
        self._steering_rod_point = steering_rod_point
        self._steering_rack_point1 = steering_rack_point1
        self._steering_rack_point_final = steering_rack_point_final
        self._steering_rod_ip = np.zeros([len(self._unit_upright_vector),3])
        self._steering_rod_ip[:,0] = np.linspace(self._steering_rack_point1[0],self._steering_rack_point_final[0],len(self._unit_upright_vector))
        self._steering_rod_ip[:,1] = np.linspace(self._steering_rack_point1[1],self._steering_rack_point_final[1],len(self._unit_upright_vector))
        self._steering_rod_ip[:,2] = np.linspace(self._steering_rack_point1[2],self._steering_rack_point_final[2],len(self._unit_upright_vector))
        return
    
    def steering_rod_turning(self):
        
        dummy = 0
        theta_vectors = []
        s_r_t = []
        theta_vector = np.zeros(len(self._unit_upright_vector))
        vecs = np.zeros([len(self._unit_upright_vector),3])
        datum_basis_horiz = np.zeros([len(self._unit_upright_vector),3])
        n_basis = np.zeros([len(self._unit_upright_vector),3,3])
        n = np.zeros([len(self._unit_upright_vector),3])
        Vec_1 = np.zeros([len(self._unit_upright_vector),3])
        
        vec = self._top_A_arm_arc[0] - self._steering_rod_point
        n[0] = self._unit_upright_vector[0]*((np.dot(vec,self._unit_upright_vector[0])))
        n[0:] = norm(n[0])*self._unit_upright_vector
        s_pivot_points = self._top_A_arm_arc - n
        
        Vec_1[0] = self._steering_rod_point - s_pivot_points[0]
        
        datum_basis_horiz[:] = np.array([1,0,0])
        horiz = np.array(list(map(cross, n, datum_basis_horiz)))
        horiz_alt = np.array(list(map(cross, horiz, n)))
        n_basis = np.array([np.array(list(map(div, horiz_alt, (np.array(list(map(norm, horiz_alt))))))), \
                            np.array(list(map(div, horiz, (np.array(list(map(norm, horiz))))))), \
                            np.array(list(map(div, n, (np.array(list(map(norm, n)))))))])
            
        nToVec_1_basis_find = lambda t: norm(((n_basis[0,0] * np.cos(t)) + (n_basis[0,1] * np.sin(t))) - (Vec_1[0]/(norm(Vec_1[0]))))
        nToVec_1_basis = root(nToVec_1_basis_find, 0)
        unit_Vec_1 = ((n_basis[0,:] * np.cos(nToVec_1_basis.x)) + (n_basis[1,:] * np.sin(nToVec_1_basis.x)))
        Vec_1 = norm(Vec_1[0])*unit_Vec_1

        for i, a in enumerate(s_pivot_points):
            vecs[i], dummy = C1.suspension_rod(Vec_1[i], n[i], a, self._steering_rod_length, np.array([self._steering_rack_point1]), np.array([dummy])).get_rod_points()
            
            S_R_T, theta_vector = C1.suspension_rod((vecs[i] - a), n[i], a, self._steering_rod_length, self._steering_rod_ip, theta_vector).get_rod_points()
            s_r_t.append(S_R_T)
            theta_vectors.append(theta_vector)
            
        return np.asarray(s_r_t), self._steering_rod_ip, s_pivot_points, np.asarray(theta_vectors)
    
#Now with the previous class' outputs, we can determine the castor, toe, camber, and steering angles all at once.
class Angles():
    
    def __init__(self, upright_vectors, unit_upright_vectors, steering_rod_turning, steering_rod_pivots, centerline_point):
        self._upright_vector = upright_vectors
        self._unit_upright_vector = unit_upright_vectors
        self._steering_rod_turning = steering_rod_turning
        self._s_pivot_points = steering_rod_pivots
        self._centerline_vec = centerline_point/(float(norm(centerline_point)))
        self._alt_centerline_vec = (self._centerline_vec+np.array([(4/5),(3/5),0]))/(np.sqrt(np.dot(self._centerline_vec+np.array([(4/5),(3/5),0]),self._centerline_vec+np.array([(4/5),(3/5),0]))))
        self._centerline_vert_vec = np.abs(cross(self._centerline_vec,self._alt_centerline_vec))
        self._centerline_n_vec = np.abs(cross(self._centerline_vert_vec,self._centerline_vec))
        return
        
    def toe_angles(self):
        
        toe_angles = []
        
        for i in enumerate(self._steering_rod_turning):
            arr = np.zeros([len(i[1])])
            for j, a in enumerate(i[1]):
                b_0 = a - self._s_pivot_points[i[0]]
                b = b_0 - np.array([0,0,b_0[2]])
                if b[1] > 0:
                    arr[j] = -(np.arccos((np.dot((b/(norm(b))),self._centerline_vec))))
                    
                elif b[1] < 0:
                    arr[j] = (np.arccos((np.dot((b/(norm(b))),self._centerline_vec))))
            toe_angles.append(arr)
        
        return np.asarray(toe_angles)
        
    def camber_angles(self):
        
        camber_angles = []
        centerline_n_vec_mag = np.sqrt(np.dot(self._centerline_n_vec,self._centerline_n_vec))
        
        for i in enumerate(self._steering_rod_turning):
            arr = np.zeros([len(i[1])])
            for j, a in enumerate(i[1]):
                b_0 = a - self._s_pivot_points[i[0]]
                b_1 = cross(self._unit_upright_vector[i[0]], b_0)
                b = b_1 - np.array([b_1[0],0,0])
                arr[j] = (-(np.arccos((np.dot(b,self._centerline_n_vec))/ \
                                ((norm(b)) * (centerline_n_vec_mag)))))
            camber_angles.append(arr)
            
        return np.asarray(camber_angles)
        
    def castor_angles(self):
        
        castor_angles = np.zeros([len(self._upright_vector)])
        centerline_vert_vec_mag = norm(self._centerline_vert_vec)
        
        for j, a in enumerate(self._upright_vector):
            b = a - np.array([0,a[1],0])
            if b[0] > 0:
                castor_angles[j] = -(np.arccos((np.dot(b,self._centerline_vert_vec))/((norm(b)) * \
                                                                (centerline_vert_vec_mag))))
            elif b[0] < 0:     
                castor_angles[j] = np.arccos((np.dot(b,self._centerline_vert_vec))/((norm(b)) * \
                                                                (centerline_vert_vec_mag)))
            
        return castor_angles
        
    def steering_angles(self, toe_angles):
        
        steering_angles = np.zeros([len(toe_angles)])
        
        for j, a in enumerate(toe_angles):
            if a[0] > a[1]:
                steering_angles[j] = np.abs(a[-1] - a[0])
            elif a[0] < a[1]:
                steering_angles[j] = np.abs(a[-1] + a[0])
                
        return steering_angles
                    
class Vars():
    
   def __init__(self, upright_point_top, center_point_top, upright_point_bot, center_point_bot, n_top, n_bot, susp_rod_ip1, \
                    points_count, which_arm, steering_rod_length, steering_rod_point, steering_rack_point1, \
                        steering_rack_point_final, centerline_point):
       
        self._upright_point_top = upright_point_top
        self._upright_point_bot = upright_point_bot
        self._n_top = n_top
        self._n_bot = n_bot
        self._center_point_top = center_point_top
        self._center_point_bot = center_point_bot
        self._susp_rod_ip1 = susp_rod_ip1
        self._points_count = points_count
        self._which_arm = which_arm
        self._steering_rod_length = steering_rod_length
        self._steering_rod_point = steering_rod_point
        self._steering_rack_point1 = steering_rack_point1
        self._steering_rack_point_final = steering_rack_point_final
        self._steering_rod = self._steering_rod_point-self._steering_rack_point1
        self._centerline_point = centerline_point
        return
    
   def Run(self):        
       vars_dict = {"susp_rod_initial_points":[],
                    "top_aarm_arc":[],
                    "center_point_top":np.linspace(self._center_point_top,self._center_point_top,self._points_count),
                    "bot_aarm_arc":[],
                    "center_point_bot":np.linspace(self._center_point_bot,self._center_point_bot,self._points_count),
                    "bot_aarm_arc_angles":[],
                    "upright_vector":[],
                    "upright_vector_mag":[],
                    "steering_rod_turning":[],
                    "steering_rod_pivots":[],
                    "steering_rack_turning":[],
                    "toe_angles":[],
                    "camber_angles":[],
                    "castor_angles":[],
                    "steering_angles":[]
                    }


       A_arm_obj = A_arms(self._upright_point_top, self._upright_point_bot, self._n_top, self._n_bot, self._center_point_top, self._center_point_bot, self._susp_rod_ip1, self._points_count, self._which_arm)
       
       
       top_arc, top_arc_theta_vector = A_arm_obj.top_A_arm_arc()
       bot_arc, bot_arc_theta_vector = A_arm_obj.bot_A_arm_arc(top_arc)
       vars_dict["upright_vector"], unit_upright_vector = A_arm_obj.upright_vectors(top_arc, bot_arc)
       if self._which_arm == 0:
           vars_dict["susp_rod_initial_points"] = (A_arm_obj.suspension_rod_initial_points(top_arc_theta_vector)) 
       elif self._which_arm == 1:
           vars_dict["susp_rod_initial_points"] = (A_arm_obj.suspension_rod_initial_points(bot_arc_theta_vector)) 
       
       
       Steering_rod_obj = Steering_rod(vars_dict["upright_vector"], unit_upright_vector, top_arc, self._steering_rod_length, self._steering_rod_point, self._steering_rack_point1, self._steering_rack_point_final)
       
       
       vars_dict["steering_rod_turning"], vars_dict["steering_rack_turning"], vars_dict["steering_rod_pivots"], s_theta_vectors = \
           Steering_rod_obj.steering_rod_turning()
       
       
       Angles_obj = Angles(vars_dict["upright_vector"], unit_upright_vector, vars_dict["steering_rod_turning"], vars_dict["steering_rod_pivots"], self._centerline_point)
       
       
       vars_dict["toe_angles"] = Angles_obj.toe_angles()
       vars_dict["camber_angles"] = Angles_obj.camber_angles()
       vars_dict["castor_angles"] = Angles_obj.castor_angles()
       vars_dict["steering_angles"] = Angles_obj.steering_angles(vars_dict["toe_angles"])
       vars_dict["upright_vector_mag"] = np.array(list(map(norm, vars_dict["upright_vector"])))
       vars_dict["top_aarm_arc"] = top_arc
       vars_dict["bot_aarm_arc"] = bot_arc
       vars_dict["bot_aarm_arc_angles"] = bot_arc_theta_vector
       
       if True:
           print('Success')
    
       return vars_dict
   