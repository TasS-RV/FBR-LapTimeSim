import numpy as np
from scipy.interpolate import interp1d
import scipy.optimize as opt
import matplotlib.pyplot as plt
import csv

from scipy.interpolate import splrep, BSpline

class powertrain():
    def __init__(self):
        self.engine_data = 'R6_power.csv' # data file with engine rpm, torque, power in columns
        self.final_drive = 1 # final drive ratio
        self.ratios_0 = [1] #Allows custom gear ratio to be applied - mainly needed for EV and convenience (TS) - "tidy up later"
        
        self.ratios = self.final_drive/(1*np.array(self.ratios_0)) #gear ratios
        self.shift_time = 0.5 #shift time in seconds

        #setting up parameters to be set by engine() method
        self.f_power = None # interp1d function for engine power = func(rpm)
        self.f_torque = None # interp1d function for engine torque = func(rpm)
        self.max_rpm = 0
        self.min_rpm = 0
        
        self.ic = True  #True if car is combustion, false if EV (combustion by default to avoid errors)

        self.engine()
        #Inclusion of parameter - launch RPM, non zero for CV, 0 for 2-step condition of EV 
        
    def update(self):
        '''updates engine interp1d functions and gear ratios'''
        self.engine()
        self.ratios = self.final_drive/1*np.array(self.ratios_0) #(TS) gear ratios - updated by input of custom Transmission ratios

    def engine(self):
        '''reads in engine data from a csv'''
        path = self.engine_data
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter=',')
            headers = next(reader)
            data = np.array(list(reader)).astype(float)
        
        #Number of divisions to split x_values into for spline smoothing
        n = 10
        #interpolator makes data into linear piecewise functions for power and torque
        
        # if self.ic == True:
        #     self.f_power = interp1d(data[:,0],data[:,2])
        #     self.f_torque = interp1d(data[:,0],data[:,1])
        
        # elif self.ic == False: #For Motors - sharp drop-off can result in ValueError - smoothing curve using spline tool
        #     power_curve = BSpline(*splrep(data[:,0], data[:,2], s=round(len(data[:,0])/n)))(data[:,0])  #data[:,0] stores all the x-axis values of RPM
        #     torque_curve = BSpline(*splrep(data[:,0], data[:,1], s=round(len(data[:,0])/n)))(data[:,0])  
        #     #Reconverts back to smoother interpolator function

        #     self.f_power = interp1d(data[:,0], power_curve)
        #     self.f_torque = interp1d(data[:,0], torque_curve)
        
        #Comment the above if-elif block to revert to the previous default interpolator functions
        self.f_power = interp1d(data[:,0],data[:,2])
        self.f_torque = interp1d(data[:,0],data[:,1])
        


        self.max_rpm = data[-1,0]
        self.min_rpm = data[0,0]
'''(TS)Must modify the interpolator function to be able to accomodate Power and Torque curves - currently this may be getting generated incorrectly'''


class car():
    def __init__(self):
        #mass parameters
        self.x_wb = 1.55 # wheelbase length in m
        self.m = 310 # mass in kg
        self.cg = (1.55*0.525,0.264) #coordinates of center of gravity (x,z) in m

        #aero parameters
        self.k_df = 100*9.81/(0.5*1.225*25**2) #downforce coefficient, with reference area = 1m^2
        self.k_drag = 0.872 #drag coefficient, with reference area = 1m^2
        self.cp = (0.5*1.55,0.3) #center of pressure in m

        #power, braking
        self.mu = 1.6 # coefficient of friction on tyres
        self.powertrain = powertrain()

        #accel, braking velocity distributions
        self.brake_v = None #interp1d function of velocity against distance while braking v = fun(x)
        self.brake_v_inv = None #interp1d function of distance while braking against velocity x = fun(v)
        self.accel_v = None #interp1d function of velocity against distance while accelerating v = fun(x)
        self.accel_v_inv = None #interp1d function of distance while accelerating against velocity x = fun(v)
        self.update()

    def update(self):
        '''recalculates the interp1d functions'''
        self.braking(plot=False)
        self.acceleration(verbose=False)

    def a_brake(self, v):
        '''returns the braking acceleration for a given velocity'''
        return (0 - (self.m*9.81 +0.5*1.225*v**2*self.k_df)*self.mu - 0.5*1.225*(v*abs(v))*self.k_drag)/self.m

    def braking(self, plot=False):
        '''calculates braking from 100ms^-1 down to 0, sets the interpolator for braking velocity'''
        v = np.linspace(100,0,500)
        x = np.zeros(len(v))
        dv = v[1]-v[0]
        #This is where the integration of the velocity curve happens
        for i in range(1,len(v)):
            x[i] = x[i-1] + 0.5*(v[i]+v[i-1])/self.a_brake(0.5*(v[i]+v[i-1])) * dv
        #now we'll add on some linear sections on either end to expand the input range of the interpolators
        v_ext = np.append(np.array([100]),v)
        v_ext = np.append(v_ext,np.array(-100))
        x_ext = np.append(np.array([-500]),x)
        x_ext = np.append(x_ext, np.array([x[-1]+300]))
        #this then just plots the curve if you want to see it
        if plot:
            print(x[-1])
            plt.plot(x_ext,v_ext)
            #plt.plot(x_ext,np.sqrt(100**2-2*9.81*self.mu*x_ext))
            plt.show()
        #setting the interpolators
        f_v_brake = interp1d(x_ext,v_ext)
        self.brake_v = f_v_brake
        f_v_brake_inv = interp1d(v_ext,x_ext)
        self.brake_v_inv = f_v_brake_inv

    def engine_aero_accel(self, v, info=False):
        '''This function calculates the aceleration for a given speed. info returns data about what's the limiting acceleration.'''
        drag_accn = self.k_drag*0.5*1.225*v**2/self.m

        gear = 0    
        rpm = v*self.powertrain.ratios*60/(2*np.pi*0.2445) #conversion between car speed and engine rpm
        #assume a launch with 2-step and thus engine rpm of 5000 rpm
        if np.min(rpm) <= 5000: 
            engine_accn = self.powertrain.f_torque(5000)*self.powertrain.ratios[0]/(0.2445*self.m) #data[0,1]
            if info:
                gear = -1   #Simply for convenience - shows we are grip limited
        # else return highest accelleration achievable in any of the gears 

#_______________________________________________ Launch clause - based on slow speed, expectation on clutch

        else:
            gear_accns=np.zeros(len(rpm))
            for i in range(len(rpm)):
                if rpm[i]<=self.powertrain.max_rpm:  #rpm in gear i - acceleration in ith gear (finds gear, for whichever gives max acc.)
                    gear_accns[i]=self.powertrain.f_power(rpm[i])/(v*self.m)
            gear = np.argmax(gear_accns)+1  #+1 due to indexing
            engine_accn = np.max(gear_accns)
            #if we're on the rev limiter then constant speed
            if rpm[-1]>=self.powertrain.max_rpm:
                gear = 5
                engine_accn = drag_accn  #Forces to 0 (as drag_acc. is subtracted off)

        #this is the function derived in nuclino documentation
        a_grip = self.mu/self.m*(self.m*9.81*self.cg[0] + 0.5*1.225*v**2*(self.k_df*self.cp[0]+self.k_drag*(self.cp[1]-self.cg[1])))/(self.x_wb-self.mu*self.cg[1])-0.5*1.225*v**2*self.k_drag/self.m #mu*R/m, R=(downforce+rear tyre load due to inertia)
        if info:
            if engine_accn < a_grip:  
                return engine_accn - drag_accn, gear
            else:
                return a_grip, 0    #Nothing else - 0 simply just to not indicate a gear, and shwo you are grip limited
        else:
            return min(engine_accn- drag_accn, a_grip)  #Worst case scenario (taking into account weight transfer to rear tyres etc...)

    def acceleration(self, verbose=False):
        '''Calculates 500m of acceleration from a standing start and sets the interpolator for acceleration'''
        x = np.linspace(0,500,2000)
        dx = x[1]-x[0]
        v=np.zeros(len(x))
        a = v.copy()  
        v[0]=0
        limit=np.zeros(len(x)) # this will keep track of what's the limiting acceleration - grip=0, 1-5=1st-5th gear
        shift_x = 0
        #integrating over x
        for i in range(1,len(x)):
                accel, limit[i] = self.engine_aero_accel(v[i-1], info=True)
                if limit[i] == limit[i-1] + 1 and limit[i-1] > 0:
                    #a gear shift has taken place
                    shift_x = x[i] + v[i-1]*self.powertrain.shift_time #assuming no significant speed change during shift
                if x[i] < shift_x:
                    #we are still during the shift
                    drag_accn = self.k_drag*0.5*1.225*v[i-1]**2/self.m
                    accel = 0
                v[i] = np.sqrt(v[i-1]**2 + 2*dx*accel)
                a[i] = accel
        
        #plotting function
        if verbose:
            t=np.zeros(len(x))
            for i in range(1,len(x)):
                t[i] = t[i-1] + 2/(v[i-1]+v[i])*dx
            fig,ax = plt.subplots(2,2)
            fig.set_tight_layout(True)
            ax[0,0].plot(x,v)
            ax[0,0].plot(x,10*limit)
            ax[0,0].set(xlabel='x/m', ylabel='v/ms^-1', title='Time = {}'.format(round(t[-1],3)))
            ax[0,1].plot(t,x)
            ax[0,1].plot(t,10*limit)
            ax[0,1].set(xlabel='t/s', ylabel='x/m')
            rpm = []
            for i in range(len(x)):
                if limit[i] -1 <=0:
                    gear = 0
                else:
                    gear = limit[i]-1
                rpm.append(v[i]*self.powertrain.ratios[int(gear)]*60/(2*np.pi*0.2445)) 
            ax[1,0].plot(x,rpm)
            ax[1,0].plot(x,limit)
            ax[1,0].set(xlabel='x/m', ylabel='engine rpm')
            ax[1,1].plot(t,a)
            ax[1,1].plot(t,limit)
            ax[1,1].set(xlabel='t/s', ylabel='a/ms^-2')
            plt.show()

        #extend the range of the functions by adding on some linear sections
        #this section is just to allow optimisers like a newton-raphson solver to converge on v=0, x=0
        x_ext = np.append(np.array([-500]),x)
        v_ext = np.append(np.array([-50]),v)
        # assume after 500m, the velocity is constant
        x_ext = np.append(x_ext,np.array([1000]))
        v_ext = np.append(v_ext,np.array(v[-1]))

        #set the interpolators
        f_v_accel = interp1d(x_ext,v_ext)
        self.accel_v = f_v_accel
        f_v_accel_inv = interp1d(v_ext,x_ext)
        self.accel_v_inv = f_v_accel_inv

    def cornering_speed(self, r_turn):
        '''returns the grip limited speed of the car around a corner of radius r_turn'''
        if r_turn < 2*self.m/(self.mu*1.225*self.k_df+0.0000001): #0.0000001 here just to avoid divide by zero error when df=0
            vmax = np.sqrt(self.mu*self.m*9.81/(self.m/r_turn-0.5*self.mu*1.225*self.k_df)) #as defined in nuclino docs
        else:
            vmax = 90 #if we're above critical radius just set it to be our "no upper limit" speed
        return np.min([90,vmax]) # cap it to 90m/s
    