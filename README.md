# FBR-LapTimeSim
LapTime Simulator for Full Blue Racing. Simplified Physics model whcih also includes aerodynamic drag, alongside grip limited cornering. 

## Getting Started

This branch is mainly oriented towards the EV motor performance modelling:

```bash
# Clone the repository
git clone https://github.com/TasS-RV/FBR-LapTimeSim

# Navigate into the project directory
cd ./"Lap Time Sim"

# Install dependencies
pip install -r requirements.txt

# Run the main script
python Run_Car.py
```
This will produce the following results:

1. Coloured map with Velocity profile of the vehicle over time on the track, at each Sector of the track.
2. Duty cycle - this is again the velocity profile, alongside a Power vs time graph. Note, for knee-points of the Torque graph, it is NOT uncommon to see the Power graph looking like a square wave pulse, where the car is almost constantly power limited rather than torque limited.
3. Duty cycle - this is the 'Current' profile. The current going in can also be negative during braking, which corresponds to regenerative braking. Note, for now the I_rms and V_max are set to 200 A and 400V respectively as some rough parameters of the Battery Pack.

The duty cycles can be found in: the datafile and the graphs
DutyCycle_Motor9_FSA_track.csv
DutyCycle_Motor9_FSA_track.png

As an example, the above graphs and datafiles can be found, as the run was done with the Motor9.csv.  You can change which motor to use, by simply changing the motor_num in the Run_Car.py file.
```python
FBRev.powertrain.final_drive = 12
gear_ratios = [1] 
motor_num = 9
FBRev.powertrain.engine_data = f"Motor{motor_num}.csv"
```
