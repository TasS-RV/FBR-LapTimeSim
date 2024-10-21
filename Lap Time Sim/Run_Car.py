from Tracks import *
from Parameters import motors_list, I_rms, V_max


def tp_curve_check(car):
    rpm_values = np.linspace(car.powertrain.min_rpm,car.powertrain.max_rpm, 1000)
    power_y = FBR27.powertrain.f_power(rpm_values)
    torque_y = FBR27.powertrain.f_torque(rpm_values)
    plt.plot(rpm_values, power_y/1000, '-', rpm_values, torque_y, '-')
    plt.show()


def motor_optimise(trackfile, CarInstance, motors_dict, motors_sub = []):
    # Purposefully truncates the number of motors, to minimise the time required for the runs to setup - this is if we don't necessarily want to run the optimisation code.
    
    if len(motors_sub) == 0:
        motors_subsample = motors_list #If no motors are specified, run all of them - from the Parameters.py file
    else:
        motors_subsample = {key: motors_dict[key] for key in motors_sub if key in motors_dict}
    
    fastest_time = 1000 #Arbitrarily high lap-time: impossible

    # Comment this out to do the iterative run
    for n, motor in enumerate(motors_subsample, 1):
        FBRev.powertrain.engine_data = f"Motor{n}.csv"  
        FBRev.powertrain.update() #Updates file reading
        FBRev.update()
        time = process_track(read_track(trackfile),FBRev, verbose=False)[0] 
        
        if time < fastest_time:
            fastest_time = time
            best_motor = n #nth motor - required for auto-generating the track performance for the fastest motor

    FBRev.powertrain.engine_data = f"Motor{best_motor}.csv"  
    FBRev.powertrain.update() 
    FBRev.update()

    print("Best motor is: {} with a lap-time of: {:.2f}s.\n Total energy consumption: {:.2f} MJ".format(best_motor, fastest_time, process_track(read_track(trackfile),FBRev, verbose=False)[1]/(1e6)))



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
motor_num = 4
FBRev.powertrain.engine_data = f"Motor{motor_num}.csv"


'''
1. To check the Motor torque-speed curve, run the tp_curve_check function with the FBRev object as the argument.
2. Set the track file to the desired track file - options include: FSA_track.dxf, boring_track.dxf, chicane.dxf.
'''
trackfile = "FSA_track.dxf"
#tp_curve_check(FBRev) #uncomment

FBRev.powertrain.update() 
FBRev.update()

#Will be an array of motor names, which are esentially the keys in the motors_list dictionary.
#["motor2"]


if __name__ == "__main__":
    results = process_track(read_track(trackfile),FBRev, verbose=False)
    energy_consumed = results[1]/(1e6)
    fastest_time = results[0]

    print(f"____Configuration information and results:_____\n\n"
        f"Motor file: {motor_num}\n"
        f"Track file: {trackfile}\n"
        f"Fastest Lap time: {fastest_time}\n"
        f"Energy consumption: {energy_consumed} MJ\n")
    


'''
Need to also implement a way to obtain the drive cycle, as both a plot and a .csv file.
'''