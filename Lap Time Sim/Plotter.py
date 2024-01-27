from Parameters import FS_car, motor, motors_list, FBR23
import math
import matplotlib.pyplot as plt
import numpy as np

#The following is to bypass maximum recursion depth error in the s_net() computation
#Allows up to 10s acceleration time - with dt (timestep) of 0.005s
import os, sys
sys.setrecursionlimit(10000) 


#_____Plotting and Save functions_____

graphs_folder = "Graphs"
if not os.path.exists(graphs_folder):
    os.makedirs(graphs_folder)


def plot_save(motor_name, file_loc):
    '''
    Saves each plot - unique identifier is motor name, with peak power/ torque characteristics to distinguish model from datasheet.
    '''
    plot_filename = os.path.join(graphs_folder, "{}.png".format(motor_name.name))
    plt.savefig(plot_filename)


def graph_plot(motor_, displacements, acc_times):
    '''
    General plotting function following iterative calculation of displacement:

    -Plots the following graphs: Power, Torque (both over motor speed range) and Displacement-time graph
    -Contains a minor, janky fix for mismatched dimensions in the Power-Torque curves, but should not affect the dynamics calculations.
    -States the time to completion for the Acceleration event.
    '''
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 6)) #Refresh figure plane to 

    #Bit of a janky fix for mixmatched dimensions, but due to high resolution, speed increments are small enough to neglect removal of highest speeds
    if (len(motor_.powers) or len(motor_.torques)) != len(motor_.speeds):
        array_lowest = min(len(motor_.torques), len(motor_.speeds), len(motor_.powers))
        torques = motor_.torques[:array_lowest-1] 
        speeds = motor_.speeds[:array_lowest-1]
        powers = motor_.powers[:array_lowest-1]
    
    else:
        torques = motor_.torques
        speeds = motor_.speeds
        powers = motor_.powers

    # Plot torque curve with 1/w decay in const. power region
    axes[0].plot(speeds, torques, label='Peak torque: {:.2f} Nm'.format(motor_.T_peak), color='blue')
    axes[0].set_title('Torque curve --- Max freq: {:.2f}'.format(motor_.w3))
    axes[0].legend()
    

    #Plot power curve with linear increase up to end of const. torque region
    axes[1].plot(speeds, powers, label='Peak power: {:.2f} W'.format(motor_.P_peak), color='red')
    axes[1].set_title('Power curve --- Peak torque/power transition frequency = {:.2f}'.format(motor_.w1))
    axes[1].legend()
    
    #Plot displacement-time curve for acceleration event
    axes[2].plot(acc_times, displacements, label='Time elapsed: {:.3f}s'.format(acc_times[-1]), color='green')
    axes[2].set_title('Displacement-Time curve (Acc. event). Motor: {}'.format(motor_.name))
    axes[2].legend()
    
    #Ensuring graphs fit with no axes/ text overlaps
    plt.tight_layout()
    plot_save(motor_, graphs_folder)
   
    #plt.show() #Can comment this out and only check folder





#____Selection of vehicle and (motor) powertrain____
car = FBR23
loss_factor = 0.88 #Using an automotive standard of 12% drivetrain losses - usually correction factor applied in dyno runs [confirm source]
dt = 0.01 #Increments in seconds


def acc(w_rpm, motor_, car):
    '''
    -Function takes in definition of motor and car parameters
    -Takes in motor spin frequency in RPM
    -Calculates net acceleration - based on motor drive torque (varies with spin frequency) and subtracts
    aerodynamic drag (scales with v^2)
    -Acceleration output given in m/s^2
    '''
    w = w_rpm*2*np.pi/60 #Motor speed given in RPM - conversion to rad/s
    vel = w*car.rw/car.g  #Car speed in m/s

    drag_force = (car.drag*car.A*1.2*(vel)**2)/2
    motor_force = motor_.torque(w)*car.g/car.rw    
    
    #Sanity check at appropiate speeds
    #print(f"{vel*3600/1000} km/hr")
     
    #Work done against aerodynamic drag
    drag_work_done = (loss_factor*motor_.P_peak - drag_force*vel) #Resistive powerloss = F_drag*velocity

    if (motor_force-drag_force)/car.mass > 0 and drag_work_done > 0:
        return (motor_force-drag_force)/car.mass #F=ma applied 
    else: 
        #Drag force matches either: 1) Maximum driving torque from motor, or 2) Work done by motor - no further acceleration
        return 0 



def s_net(vel, s_list, t_list, motor_, car):
    '''
    -Function takes in definition motor and car parameters, and initial velocity and time elapsed
    -Takes in vehicle velocity in m/s
    -Increments with s = ut + 0.5at^2 formula
    -Updates the instantaneous velocity with either ds/dt or u_prev + current_acceleration*dt
    -Increments the displacement (so far) with ds steps, and increments the total time elapsed using dt
    '''
    w_rpm = vel*car.g*60/(2*np.pi*(car.rw)) #Converting from /s speed back to RPM shaft rotation speed
    acc_current = acc(w_rpm, motor_, car) #Acceleration based on current speed
    ds = vel*dt + (1/2)*acc_current*dt**2 #s_n - current displacement increment
    #print("\n{}".format(vel))
     

    vel = ds/dt  
    vel = vel + acc_current*dt 
    #time_tot += dt #Incrementing with set timestep
    s_list.append(s_list[-1]+ds)
    t_list.append(t_list[-1]+dt) 
    
    if  s_list[-1] > 75: #Acceleration track finished - final displacement and completion time
        return s_list, t_list 
    else:
        return s_net(vel, s_list, t_list, motor_, car) 
   

'''
-Iterating through each instance of motor - with varying motor properties.
-Iteratively computing the time taken to complete 75m of track - acc() and s_net() functions used
-Also prints Maximum speed (either drag limited, or limited by rotational speed of the motor itself)
'''


for motor_ in motors_list:

    displacements_list = [0]
    times_list  = [0]

    w_list = motor_.speeds
    displacements, acc_times = s_net(0, displacements_list, times_list, motor_, car)
    
    #Pre-check in text window/ terminal: if sensible times for covering 75m track
    print("Displacement: {:.2f} m with Time taken: {:.2f} s".format(displacements[-1], acc_times[-1]))
    
    #Rotation frequency limited top speed: ceiling of 85% of max motor speed for safety
    rot_max = 0.85*(motor_.w3)*car.rw/car.g  
    #Equilibruim between motor power and resistive losses (aerodynamic only) OR limited by max. motor rotation speed
    
    top_speed = min(math.pow((2*motor_.P_peak*loss_factor)/(car.drag*car.A*1.2), (1/3)), rot_max)
    top_speed = top_speed*(3600/1000)
    
    print("Maximum speed: {:.1f} km/hr \n".format(top_speed))
    
    graph_plot(motor_, displacements, acc_times)









