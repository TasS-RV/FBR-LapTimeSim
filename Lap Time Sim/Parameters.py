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
resolution = 1000 #Affects number of increments the speed-range of the motor is split into    

motor_file = "Motors_def.xlsx"
df = pd.ExcelFile(motor_file)

# Assuming you have a sheet named "Sheet1" in your Excel file
# Access a specific table (e.g., Table1) within that sheet
table_name = 'Sheet1'
table_df = df.parse(table_name)

# Defining up to 3 motors based on contents of Excel file 
motor1 = motor(table_df['w3'][0], table_df['w1'][0], table_df['T_nominal'][0], table_df['P_nominal'][0], table_df['T_peak'][0], table_df['P_peak'][0], table_df['Motor'][0])
motor2 = motor(table_df['w3'][1], table_df['w1'][1], table_df['T_nominal'][1], table_df['P_nominal'][1], table_df['T_peak'][1], table_df['P_peak'][1], table_df['Motor'][1])
motor3 = motor(table_df['w3'][2], table_df['w1'][2], table_df['T_nominal'][2], table_df['P_nominal'][2], table_df['T_peak'][2], table_df['P_peak'][2], table_df['Motor'][2])

motors_list = [motor1, motor2, motor3]






