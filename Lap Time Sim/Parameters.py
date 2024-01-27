import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd #Read motor parameters from Excel file

class FS_car:
    def __init__(self, gear_ratio, car_mass, grip_limit, wheel_rad, drag_coeff, frontal_area):
        self.g = gear_ratio
        self.rw = wheel_rad
        self.F_limit = grip_limit
        self.mass = car_mass
        
        #Defining drag characteristis
        self.drag = drag_coeff
        self.A = frontal_area


class motor:

    def __init__(self, speed_max, trans_freq, normal_torque, normal_power, peak_torque, peak_power, obj_name):
        self.w3 = 2*np.pi*speed_max/60  #Frequencies in revolutions per second instead of minute
        self.w1 = 2*np.pi*trans_freq/60
        #Key points on frequency graph
        self.T = normal_torque
        self.T_peak = peak_torque
        self.P = normal_power
        self.P_peak = peak_power

        #Motor Identifier for labelling the correct graph
        self.name = obj_name
    
    @property
    def speeds(self):
        return np.linspace(0, self.w3, resolution)
    
    @property
    def torques(self):
        #Constant Torque up to speed for torque decay region
        const_region = np.empty(math.floor((self.w1/self.w3)*resolution))
        const_region.fill(self.T)

        #Torque decay as 1/w in constant power region
        decay_region = [(self.T*self.w1)/w for w in self.speeds[(math.ceil((self.w1/self.w3)*resolution)):]]
        decay_reg = np.array(decay_region)

        return np.concatenate((const_region, decay_reg)) 
         
    @property
    def powers(self):
        rise_region = np.array([w*self.T for w in self.speeds[0:(math.floor((self.w1/self.w3)*resolution))]])

        const_region = np.zeros(resolution - math.ceil((self.w1/self.w3)*resolution))
        const_region.fill(self.T*self.w1) #Power defined by current rotation speed - rather than electromagnetic true power output

        return np.concatenate((rise_region,const_region)) 
    
    def torque(self, w):
        return self.torques[math.floor((w/self.w3)*resolution)]
    

#_____Car properties______    

rc = 23 #Diameter of engine-side sprocket in mm
rs = 100 #Diameter of Axle sprocket in mm
gear_ratio = rs/rc
car_mass = 240 #Dry mass in Kg
grip_mu = 4 #Maximum coefficinet of friction for car grip limit
wheel_rad = 0.225 #Wheel total diameter in m - assuming 100 contact patch @ 0 deg. camber

#Ask aero to roughly model these values
drag_coeff = 0.728
frontal_area = 2  #Approcimate total frontal area in m^2

FBR23 = FS_car(gear_ratio, car_mass, grip_mu*(car_mass*9.81*0.5),wheel_rad, drag_coeff, frontal_area)
#Note - the grip limit assumes the car is on the limit of wheely-ing - front wheels slam back onto ground after movement 
   
    
#______Motor properties______
resolution = 100 #Affects number of increments the speed-range of the motor is split into    

"""Requirement: generate a number of motor torque-speed characteristics based on varying peak torques and transition frequencies"""


test_torques = np.linspace(60, 80, 4)
#test_powers = np.array([70000]) #Will not require parameter anyways as T*w = Power
w1_variable = np.linspace(1000,5000,5)

motor_list = {}
#motors_list = [motor1, motor2, motor3, motor4, motor5, motor6, motor7, motor8, motor9, motor10, motor11, motor12, motor13, motor14, motor15, motor16, motor17, motor18, motor19, motor20]

n = 1 #Count up each motor instance
for T_max in test_torques:
    for W1 in w1_variable: #The peak_power and torque figures are arbitrarily set to not need to modify __init__ function of the class!
        motor_list[f"motor{n}"] = motor(18000, W1, T_max, T_max*W1*np.pi*2/60, 1.2*T_max, 1.2*(T_max*W1*np.pi*2/60), "Motor{}".format(n))         
      #  motors_list[n-1] = motor(18000, W1, T_max, T_max*W1*np.pi*2/60, 1.2*T_max, 1.2*(T_max*W1*np.pi*2/60), "Motor{}".format(n))
        n +=1

print(motor_list['motor3'].P_peak)





