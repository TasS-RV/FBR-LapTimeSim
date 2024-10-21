from Tracks import *
from Parameters import motors_list, I_rms, V_max


def tp_curve_check(car):
    """
    Plots the power and torque curves of the car's powertrain.

    Parameters:
    -----------
    car : object
        The car object containing the powertrain information.
    """
    rpm_values = np.linspace(car.powertrain.min_rpm, car.powertrain.max_rpm, 1000)
    power_y = car.powertrain.f_power(rpm_values)
    torque_y = car.powertrain.f_torque(rpm_values)
    
    fig, ax1 = plt.subplots()

    # Plot torque on the primary y-axis
    ax1.plot(rpm_values, torque_y, 'b-', label='Torque (Nm)')
    ax1.set_xlabel('Motor Shaft Speed (RPM)')
    ax1.set_ylabel('Torque (Nm)', color='b')
    ax1.tick_params(axis='y', labelcolor='b')

    # Create a secondary y-axis for power
    ax2 = ax1.twinx()
    ax2.plot(rpm_values, power_y / 1000, 'r-', label='Power (kW)')
    ax2.set_ylabel('Power (kW)', color='r')
    ax2.tick_params(axis='y', labelcolor='r')

    # Add title and grid
    plt.title('Power and Torque Curves')
    fig.tight_layout()  # Adjust layout to prevent overlap
    plt.grid(True)
    plt.show()


def motor_optimise(trackfile, CarInstance, motors_dict, motors_sub = []):
    """
    Optimizes the motor selection for the given car instance on a specified track.

    This function iterates over a subset of motors, simulates the car's performance on the track,
    and identifies the motor that results in the fastest lap time.

    Parameters:
    -----------
    trackfile : str
        Name of the .dxf file containing the track path.
    CarInstance : object
        An instance of the car class containing the car's parameters and powertrain. This will be FBRev in this case.
    motors_dict : dict
        A dictionary containing - 'n'th motors to extract from the motors_list dictionary. If this is not passed in, by default it will run all motors.
    motors_sub : list, optional
        A list of motor identifiers to be considered for optimization. If empty, all motors in motors_dict are considered.

    Returns:
    --------
    - Prints the best motor and its corresponding lap time and energy consumption.

    Notes:
    ------
    - The function uses a global variable 'laps' which multiplies the Power consumption/ lap.
    - The function reads motor data from CSV files named "Motor{n}.csv", where n is the motor number.
    - The function assumes the existence of 'process_track' and 'read_track' functions for track processing.
    """
    global laps

    if len(motors_sub) == 0:
        motors_subsample = motors_list #If no motors are specified, run all of them - from the Parameters.py file
    else:
        motors_subsample = {key: motors_dict[key] for key in motors_sub if key in motors_dict}
    
    fastest_time = 1000 #Arbitrarily high lap-time: impossible

    # Comment this out to do the iterative run
    for n, motor in enumerate(motors_subsample, 1):
        CarInstance.powertrain.engine_data = f"Motor{n}.csv"  
        CarInstance.powertrain.update() #Updates file reading
        CarInstance.update()
        time = process_track(read_track(trackfile),CarInstance, verbose=False)[0] 
        
        if time < fastest_time:
            fastest_time = time
            best_motor = n #nth motor - required for auto-generating the track performance for the fastest motor

    CarInstance.powertrain.engine_data = f"Motor{best_motor}.csv"  
    CarInstance.powertrain.update() 
    CarInstance.update()

    print("Best motor is: {} with a lap-time of: {:.2f}s.\n Total energy consumption: {:.2f} kWh".format(best_motor, fastest_time, process_track(read_track(trackfile),CarInstance, verbose=False)[1]/(1e3)*laps/3600))



#______ This block generates the FBRev - with identical characteristics to the FBR23 drive train except:

FBRev = CarData.car()
FBRev.m=320 #25 kg Aero kit + 15 kg surplus electronics  
"""
- Rear gear teeth/ front gear teeth: currently approx. 31/ 11, wants to be higher for greater torque.
- The gear_ratios must be entered as an array - even if this is an EV, it must = [1]. For the R6 engine, it will be the 5 speed
gear ratios from the CRANK to the ENGINE-SIDE SPROCKET.

- The engine_data must be set to the Motor{n}.txt data file. This is the file in the following format:
rpm (w) | torque (Nm) | power (W)
0       | 0           |  0
Set this by changing motor_num to the desired motor number.

Data is obtained from Inetic IEV180 motor datasheet. Found in: https://app.nuclino.com/FBR/EV/EV5---Motor-and-Motor-Test-Bench-0b117bf6-17bf-40d6-9a3a-6842e60d8159

- FBRev.powertrain.final_drive = g --> g is the gear ratio between the motor sprocket and the wheel sprocket (on the differential).
- FBRev.powertrain.ic = False --> This is an electric vehicle, so the internal combustion engine is not present.
"""
FBRev.powertrain.ic = False
FBRev.powertrain.final_drive = 12
gear_ratios = [1] 
motor_num = 9
FBRev.powertrain.engine_data = f"Motor{motor_num}.csv"


FBRev.powertrain.update() 
FBRev.update()
'''
1. To check the Motor torque-speed curve, run the tp_curve_check function with the FBRev object as the argument.
2. Set the track file to the desired track file - options include: FSA_track.dxf, boring_track.dxf, chicane.dxf.
3. Set the number of laps - this will vary depending on track length, but total distance should be 22km.
'''
trackfile = "FSA_track.dxf"
show_plot = True
laps = 30

#tp_curve_check(FBRev) #Uncomment to check the torque-speed curve of the motor

#Will be an array of motor names, which are esentially the keys in the motors_list dictionary.
#["motor2"]


if __name__ == "__main__":
    # Code below runs the motor optimisation function, which iterates over each Motor Torque-Power-speed curve and finds the fastest lap time.
    motor_optimise(trackfile, FBRev, motors_list, motors_sub = ["motor1","motor2","motor3","motor4"])

    # Code below is to run a single instance
    results = process_track(read_track(trackfile),FBRev, verbose=False)
    energy_consumed = results[1]/(1e6)
    fastest_time = results[0]

    process_track(read_track(trackfile),FBRev, verbose= show_plot)
    print(f"____Configuration information and results:_____\n\n"
        f"Motor file: {motor_num}\n"
        f"Track file: {trackfile}\n"
        f"Fastest Lap time: {fastest_time}\n"
        f"Energy consumption/ lap: {energy_consumed} MJ\n")
    


'''
Need to also implement a way to obtain the drive cycle, as both a plot and a .csv file.
'''