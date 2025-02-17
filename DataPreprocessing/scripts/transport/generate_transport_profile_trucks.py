import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# ToDos:
# Cr, c1, c2, P_aux ... : Do they have paper as source or figure out. Find P_aux from comparison with known power data.
# Check if plausible that electric consumption so much lower than ICEV engine power.
# How to allow for timer? Relative to arrival time?
# weekend_work_l2_14


# ADD MISSING ROWS/TIMESTEPS IN DRAYAGE FILE:
"""
# Data usually missing breaks! Use average of timesteps before and after missing data
input_file_path = 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/Fleet DNA Drayage Representative_original.csv'
output_file_path = 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/Fleet DNA Drayage Representative_interpolated.csv'
df = pd.read_csv(input_file_path, index_col=0)

# Find missing indices
expected_indices = list(range(30957))
missing_indices = list(set(expected_indices) - set(df.index))
print(missing_indices[-1])

# Create list of lists where consecutive missing indices are grouped.
# Example: Missing indices: 3,4,7,12,13,14 => [[3,4],[7],[12,13,14]]
input_list = missing_indices

result_list_of_lists = []
current_group = []

for i, num in enumerate(input_list):
    if i == 0:
        # For the first element, start a new sublist
        current_group.append(num)
    elif num - input_list[i - 1] == 1:
        # If the current number is consecutive to the previous one, add to the current sublist
        current_group.append(num)
    else:
        # If not consecutive, start a new sublist
        result_list_of_lists.append(current_group.copy())
        current_group = [num]
result_list_of_lists.append(current_group)
print(result_list_of_lists)

for missing_block in result_list_of_lists:
    for time in missing_block:
        df.loc[time] = (df.loc[missing_block[0]-1] + df.loc[missing_block[-1]+1]) / 2
        df.loc[time] = df.loc[time].round(2)

df_sorted = df.sort_index()
df_sorted.to_csv(output_file_path, sep=',', decimal='.')
"""


# PROCEDURE TO OBTAIN START TIME DISTRIBUTION
"""
# Read trip start times from "Collection of activity data from on-road heavy-duty ..." (https://ww2.arb.ca.gov/sites/default/files/classic/research/apr/past/13-301.pdf), Appendix E, for line haul out-of-state, line haul in-state, drayage north, drayage South, food distribution, beverage distribution, Local moving, urban buses
# This gives "start_time_distribution_raw.csv"
# Isolate start-up peaks manually for vehicle classes that have many, short trips per day (assuming each vehicle stops same number of times within one hour: means that step increase means additional shift starters): Drayage North (7am peak), Drayage South (6 and 17 peak; BESTES BEISPIEL), Beverage distribution (3-9 peak), Local moving (7 peak), Buses (4am peak)
# Keep profiles for vehicle classes that have fewer, longer trips or operate around the clock: Line haul out-of-state, Line haul in-state, Food distribution
# This gives "start_time_distribution_adjusted.csv"

# Now: Normalize each column (transport mpde) to 1. Sum up line-haul in-state and local moving (will get regional haul driving cycle). Sum up Drayage north and south (Drayage cycle). Sum up Food and Beverage distribution (Delivery cycle). Normalize again and round. Change last row so that sum is 1.
# Store this as "start_time_distribution.csv"
path = 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/start_time_distribution_adjusted.csv'
start_times = pd.read_csv(path, index_col=0)
start_times = start_times.div(start_times.sum(axis=0), axis=1)
start_times["Veh_Road_LongHaul_Elec"] = start_times["Line Haul - Out of State"]
start_times["Veh_Road_TransitBus_Elec"] = start_times["Urban Buses"]
start_times["Veh_Road_RegHaul_Elec"] = start_times["Local Moving"] + start_times["Line Haul - In State"]
start_times["Veh_Road_Drayage_Elec"] = start_times["Drayage - Northern CA"] + start_times["Drayage - Southern CA"]
start_times["Veh_Road_LocDelTruck_Elec"] = start_times["Food Distribution"] + start_times["Beverage Distribution"]
#start_times = start_times.div(start_times.sum(axis=0), axis=1).round(3)
start_times = start_times.div(start_times.sum(axis=0), axis=1)
#start_times.iloc[-1] = 1 - start_times.cumsum().iloc[-1] # Change last row, so that sum is 1.

# Discard unused columns for clarity
start_times = start_times[["Veh_Road_LongHaul_Elec","Veh_Road_RegHaul_Elec","Veh_Road_Drayage_Elec","Veh_Road_LocDelTruck_Elec","Veh_Road_TransitBus_Elec"]]
print(start_times)

# Plot
ax = start_times.plot(kind='bar', rot=0)
#for col in start_times.columns:
#    ax.plot(start_times.index, start_times[col], marker='o', linestyle='-', linewidth=1, markersize=1)
plt.xlabel('Time [hr]')
plt.ylabel('Fraction of trip starts [-]')
plt.title('Start time distributions')
#plt.legend(title='Legend Title')
plt.show()

# Save
start_times.to_csv('C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/start_time_distribution.csv')
"""


# FROM DRIVE CYCLE, GIVEN VEHICLE SPECS, COMPUTE CONSUMPTION CURVE AND SOC CURVE. APPLY CHARGING REGIME. GET CHARGING LOAD CURVE. SAVE drive_cycle AS CSV.

# Constants, vehicle data for vehicle model, paths to driving cycles
g = 9.8066 # [m/s2]
roh_air = 1.2256 # [kg/m3]

vehicle_data = {
        
        'Veh_Road_RegHaul_Elec': {
            # Based on BYD 6F truck (class 6)
            # Source 1: https://en.byd.com/truck/class-6-truck/
            # Source 2: Bayindirli, Cihan; 2016; The determination of aerodynamic drag coefficient offf truck and trailer model
            # Source 3: https://californiahvip.org/wp-content/uploads/2021/02/MY20-6F-ZE-200115.pdf
            'm': 11793.4, # Total mass of vehicle [kg]; Assumption: not a function of time; Source 1: 26k lbs
            'C_r': 1.75, # Rolling resistance param
            'c_1': 0.0328, # Rolling resistance param
            'c_2': 4.575, # Rolling resistance param
            'A_f': 5.72, # Front area of vehicle [m2]; Source 3: 94.5 in * 93.9 in = 5.72 m2
            'C_D': 0.6, # Aerodynamic drag coefficient: Source 2
            'P_aux': 0, # Power of auxiliary systems [kW]
            'eta_DL': 0.94, # Driveline efficiency
            'eta_EM': 0.95, # Electric motor efficiency
            'eta_BAT': 0.95, # Battery system efficiency
            'alpha': 0.22, # Fitting parameter for regenerative breaking efficiency, in [0,1]
            'SOC_t0': 0.8, # Assume fully charged vehicle at begin of drive cycle [-]
            'Cap_B': 281, # Battery capacity [kWh]; Source: 1

            'path_drive_cycle': 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/Fleet DNA Regional-Haul Representative_.csv',
            'min_break_for_charging': 900, # [s], 15 min break
            'specified_charge_duration': 1800, # [s], 30 min
            'P_charge': 120, # [kW]
            },
        'Veh_Road_TransitBus_Elec': {
            # Based on BYD K8M (not longest, not shortest)
            # Source 1: https://en.byd.com/bus/k8m/ 
            # Source 2: Bayindirli, Cihan; 2016; The determination of aerodynamic drag coefficient offf truck and trailer model
            'm': 17134.68, # Total mass of vehicle [kg]; Assumption: not a function of time; Source: 1, assumes half-occupied
            'C_r': 1.75, # Rolling resistance param
            'c_1': 0.0328, # Rolling resistance param
            'c_2': 4.575, # Rolling resistance param
            'A_f': 8.82, # Front area of vehicle [m2]; Source 1: 102 in x 134 in = 8.82 m2
            'C_D': 0.6, # Aerodynamic drag coefficient; Source 2
            'P_aux': 0, # Power of auxiliary systems [kW]
            'eta_DL': 0.94, # Driveline efficiency
            'eta_EM': 0.95, # Electric motor efficiency
            'eta_BAT': 0.95, # Battery system efficiency
            'alpha': 0.22, # Fitting parameter for regenerative breaking efficiency, in [0,1]
            'SOC_t0': 0.8, # Assume fully charged vehicle at begin of drive cycle [-]
            'Cap_B': 391, # Battery capacity [kWh]

            'path_drive_cycle': 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/Fleet DNA Transit Bus Representative_.csv',
            'min_break_for_charging': 300, # [s], 5 min break
            'specified_charge_duration': 1800, # [s], 30 min
            'P_charge': 150, # [kW] Source 1
            },

        'Veh_Road_LongHaul_Elec': {
            # Based on Tesla Semi
            # Source 1: https://www.tesla.com/semi
            # Source 2: Tesla Semi unveiling event
            'm': 37194.57, # Total mass of vehicle [kg]; Assumption: not a function of time, Source 1: GCVW: 82k lbs
            'C_r': 1.75, # Rolling resistance param
            'c_1': 0.0328, # Rolling resistance param
            'c_2': 4.575, # Rolling resistance param
            'A_f': 2.81, # Front area of vehicle [m2]
            'C_D': 0.36, # Aerodynamic drag coefficient; Source 2: Tesla Semi unveiling event
            'P_aux': 0, # Power of auxiliary systems [kW]
            'eta_DL': 0.94, # Driveline efficiency
            'eta_EM': 0.95, # Electric motor efficiency
            'eta_BAT': 0.95, # Battery system efficiency
            'alpha': 0.22, # Fitting parameter for regenerative breaking efficiency, in [0,1]
            'SOC_t0': 0.8, # Assume fully charged vehicle at begin of drive cycle [-]
            'Cap_B': 900, # Battery capacity [kWh] # Npt disclosed, estimates based on CEO announcements

            'path_drive_cycle': 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/Fleet DNA Long-Haul Representative_.csv',
            'min_break_for_charging': 900, # [s] 900 15 min break
            'specified_charge_duration': 1800, # [s], 30 min
            'P_charge': 750, # [kW]; Tesla Megacharger, not disclosed yet
            },
        
        'Veh_Road_LocDelTruck_Elec': {
            # Based on Rivian "Delivery 500" and original paper.
            #'m1': 3082.61, # Tare mass of vehicle [kg]; Paper: 5500; Source: Rivian "Delivery 500": GVWR = 9350 lbs - Payload = 2734 lbs = 3082.61 kg
            #'m2': 620.06, # Payload of vehicle [kg]; Assumption: not a function of time and at half load; Paper: 2500; Source: Rivian "Delivery 500": 2734 lbs / 2 = 620.06 kg
            'm': 3702.67, # Total mass of vehicle [kg], see above for more
            'C_r': 1.75, # Rolling resistance param; Paper: 1.75
            'c_1': 0.0328, # Rolling resistance param; Paper: 0.0328
            'c_2': 4.575, # Rolling resistance param; Paper: 4.575
            'A_f': 7.13, # Front area of vehicle [m2]; Paper: 2.81; Source: Rivian "Delivery 500": 114.7 in * 96.4 in = 7.13 m2
            'C_D': 0.316, # Aerodynamic drag coefficient; Paper: 0.316
            'P_aux': 0, # Power of auxiliary systems [kW]; Paper: Not specified.
            'eta_DL': 0.94, # Driveline efficiency; Paper: 0.94
            'eta_EM': 0.95, # Electric motor efficiency; Paper: 0.95
            'eta_BAT': 0.95, # Battery system efficiency; Paper: 0.95
            'alpha': 0.22, # Fitting parameter for regenerative breaking efficiency, in [0,1]; Paper: 0.22
            'SOC_t0': 0.8, # Assume fully charged vehicle at begin of drive cycle [-]
            'Cap_B': 135, # Battery capacity [kWh]; Paper: ca. 50 kWh; Source: Rivian EDV: https://www.motortrend.com/news/checking-in-on-the-rivian-amazon-edv-electric-van/

            'path_drive_cycle': 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/Fleet DNA Local Delivery Representative_.csv',
            'min_break_for_charging': 900, # [s] ,900, 1e7 longer than a day, dont want LocDelTruck to charge
            'specified_charge_duration': 1800, # [s], 30 min
            'P_charge': 100, # [kW] Source: Rivian "Delivery 500", Wikipedia says 150kW
            # CAPEX: 83000 USD minus 7500 USD tax credit
            },
        'Veh_Road_Drayage_Elec': {
            # Based on BYD 8TT (3rd generation)
            # Source 1: https://en.byd.com/news/byd-and-einride-sign-largest-ever-order-for-heavy-duty-battery-electric-trucks-outside-of-asia/
            # Source 2: https://californiahvip.org/wp-content/uploads/2020/09/2019-BYD-8TT-Cut-Sheet-190306.pdf 
            # Source 3: Bayindirli, Cihan; 2016; The determination of aerodynamic drag coefficient offf truck and trailer model
            'm': 47627, # Total mass of vehicle [kg]; Assumption: not a function of time; Source 1: BYD: GCWR: 47627kg
            'C_r': 1.75, # Rolling resistance param
            'c_1': 0.0328, # Rolling resistance param
            'c_2': 4.575, # Rolling resistance param
            'A_f': 7.86, # Front area of vehicle [m2]; Source 2: BYD: 100.4 in * 121.3 in = 7.86 m2
            'C_D': 0.6, # Aerodynamic drag coefficient; Source 3
            'P_aux': 0, # Power of auxiliary systems [kW]
            'eta_DL': 0.94, # Driveline efficiency
            'eta_EM': 0.95, # Electric motor efficiency
            'eta_BAT': 0.95, # Battery system efficiency
            'alpha': 0.22, # Fitting parameter for regenerative breaking efficiency, in [0,1]
            'SOC_t0': 0.8, # Assume fully charged vehicle at begin of drive cycle [-]
            'Cap_B': 563, # Battery capacity [kWh]; Source 1

            'path_drive_cycle': 'C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/Fleet DNA Drayage Representative_interpolated.csv',
            'min_break_for_charging': 900, # [s] 15 min
            'specified_charge_duration': 1800, # [s], 30 min
            'P_charge': 185, # [kW]
            },
        }

"""
def identify_breaks(input_vector):

    # Returns vector that is 0 for all time steps with no break and the length of the break for all time steps that are part of a break.

    # Create list of tuples, with: (integer, how often that integer occurs in a row)
    # [1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 0, 0, 0, 5, 5]
    # => [(1, 2), (2, 3), (3, 1), (4, 4), (0,3), (5, 2)]
    series_lengths = []
    current_integer = input_vector[0]
    current_length = 1

    for i in range(1, len(input_vector)):
        if input_vector[i] == input_vector[i - 1]:
            current_length += 1
        else:
            series_lengths.append((current_integer, current_length))
            current_integer = input_vector[i]
            current_length = 1

    series_lengths.append((current_integer, current_length))  # Handle the last series

    # Regenerate a vector that counts how long the current series of zeroes is
    # [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 0] 
    result_vector = []

    for tuple_item in series_lengths:
        value, count = tuple_item

        if value != 0:
            result_vector.extend([0] * count)
        else:
            result_vector.extend([count] * count)

    return result_vector

# From driving cycle compute load curve
for v in vehicle_data:

    # FIRST: FROM DRIVING CYCLE COMPUTE POWER DEMAND AND SOC
    drive_cycle = pd.read_csv(vehicle_data[v]['path_drive_cycle'], index_col=0)

    # Convert vehicle speed from mph to m/s
    drive_cycle['Speed (m/s)'] = drive_cycle['Speed (mph)'] * 0.45

    # Compute acceleration for each time step and add to dataframe
    # Take difference to previous element and divide by time interval 
    time_interval = 1 # Data is provided in 1 second intervals
    drive_cycle['Acceleration (m/s2)'] = drive_cycle['Speed (m/s)'].diff() / time_interval # Returns NaN for first time step since no "previous" speed
    drive_cycle['Acceleration (m/s2)'].fillna(0, inplace=True) # Replace NaN in first line by 0

    # Compute total vehicle mass
    #m = vehicle_data[v]['m1'] + vehicle_data[v]['m2'] # Total vehicle mass [kg]
    m = vehicle_data[v]['m']

    # Compute instantaneous power at wheels [kW] for every time step and add to dataframe (see Fiori et al., 2018)
    drive_cycle['P_W (kW)'] = (
        m * drive_cycle['Acceleration (m/s2)']
        + m * g * np.cos(drive_cycle['Grade (rise/run)']) * vehicle_data[v]['C_r'] / 1000 * (vehicle_data[v]['c_1'] * drive_cycle['Speed (m/s)'] + vehicle_data[v]['c_2'])
        + 0.5 * roh_air * vehicle_data[v]['A_f'] * vehicle_data[v]['C_D'] * drive_cycle['Speed (m/s)'] ** 2
        + m * g * np.sin(drive_cycle['Grade (rise/run)'])
    ) * drive_cycle['Speed (m/s)'] / 1000

    # Compute Regenerative breaking efficiency [-] for every time step and add to data frame
    drive_cycle['eta_RB'] = np.where(drive_cycle['Acceleration (m/s2)'] < 0, np.exp(-1 * vehicle_data[v]['alpha'] / drive_cycle['Acceleration (m/s2)'].abs()), 0)

    # Compute Instantaneous traction power [kW] for every time step and add to data frame
    drive_cycle['P_T (kW)'] = np.where(drive_cycle['P_W (kW)'] < 0, (drive_cycle['P_W (kW)'] + vehicle_data[v]['P_aux'])*vehicle_data[v]['eta_DL']*vehicle_data[v]['eta_EM']*vehicle_data[v]['eta_BAT']*drive_cycle['eta_RB'], (drive_cycle['P_W (kW)'] + vehicle_data[v]['P_aux'])/(vehicle_data[v]['eta_DL']*vehicle_data[v]['eta_EM']*vehicle_data[v]['eta_BAT'])) 
    
    # Compute (cummulative) Energy consumption [kWh] up to timestep t (in the original model the energies per kilometer are computed)
    drive_cycle['EC_cum (kWh)'] = drive_cycle['P_T (kW)'].where(drive_cycle['P_T (kW)'] > 0).fillna(0).cumsum() / 3600

    # Compute (cummulative) Energy recovery [kWh] up to timestep t (in the original model the energies per kilometer are computed)
    drive_cycle['ER_cum (kWh)'] = drive_cycle['P_T (kW)'].where(drive_cycle['P_T (kW)'] < 0).fillna(0).cumsum() / 3600
    
    # Compute state of charge [-] (in the original model the formula looks a little different since EC and ER are computed per kilometer)
    drive_cycle['SOC'] = vehicle_data[v]['SOC_t0'] - (drive_cycle['EC_cum (kWh)'] - drive_cycle['ER_cum (kWh)']) / vehicle_data[v]['Cap_B']
    SOC_copy = drive_cycle['SOC'].copy(deep=True)

    
    # SECOND: APPLY A CHARGING REGIME AND GET CONSUMPTION CURVE
    upper_SOC_limit = 0.8 # 0.8 # do not charge beyond this SOC level
    lower_SOC_limit = 0.2 # need to charge if dropping below this SOC level

    # Charge whenever resting
    # Find breaks longer than e.g. 5 min where vehicle does not drive: Find consecutive series of at least vehicle_data[v]['min_break_for_charging'] rows, i.e., seconds, with "speed" == 0. Mark those with True, all other with False.
    drive_cycle['break_length'] = identify_breaks(drive_cycle['Speed (m/s)'].to_list())
    
    # Find breaks long enough for charging
    # Replace all integers that are below "break_length" (e.g., 4s) with False, all other with True.
    # => Returns vector that is False for all time steps with no or too-short break and True for break timesteps.
    min_break_length = vehicle_data[v]['min_break_for_charging']
    drive_cycle['charging_break (bool)'] = [False if isinstance(item, int) and item < min_break_length else True for item in drive_cycle['break_length']]
    print(drive_cycle['charging_break (bool)'].sum())

    drive_cycle['SOC_wBreakCharging'] = drive_cycle['SOC']
    drive_cycle['charging_power_existing_breaks (kW)'] = 0
    drive_cycle['charging_power_depletion_breaks (kW)'] = 0
    drive_cycle['charging_power_after_cycle (kW)'] = 0
    charge_rate = vehicle_data[v]['P_charge'] # for now assuming max rate as charging rate, could use modulation etc here. For that include into iteration.
    modus = "specified_time" # once battery empty, charge for "specified_time", "until_full", or "as_necessary"
    #modus = "until_full"
    
    #for index, row in drive_cycle.iterrows():
    for index in range(172800): # two days to be long enough, iterrows() does not work since adding rows in loop
        print(index)
        alist=drive_cycle.index.tolist()
        if index not in drive_cycle.index.tolist(): # if sequence over, break
                break   

        # If this timestep qualifies for break charging:
        if drive_cycle.loc[index, 'charging_break (bool)']:
        #if row['charging_break (bool)']:

            # Compute new SOC if charging would occur:
            if index == 0:
                SOC_wBreakCharging = (vehicle_data[v]['SOC_t0'] * vehicle_data[v]['Cap_B'] + charge_rate * 1 / 3600) / vehicle_data[v]['Cap_B']
            else:
                SOC_wBreakCharging = (drive_cycle.at[index-1, 'SOC_wBreakCharging'] * vehicle_data[v]['Cap_B'] + charge_rate * 1 / 3600) / vehicle_data[v]['Cap_B']
            delta_SOC = SOC_wBreakCharging - drive_cycle.at[index, 'SOC_wBreakCharging']

            # If battery would not be overfull after receiving this charge, then charge:
            if SOC_wBreakCharging < upper_SOC_limit:
                # Update the SOC after charging for one second:
                drive_cycle.at[index, 'SOC_wBreakCharging'] = SOC_wBreakCharging
                # Update all following SOCs by this charge:
                drive_cycle.loc[(index + 1):, 'SOC_wBreakCharging'] += delta_SOC
                # Update the charging power during that time step:
                drive_cycle.at[index, 'charging_power_existing_breaks (kW)'] = charge_rate
            continue
        
        if drive_cycle.at[index, 'SOC_wBreakCharging'] <= lower_SOC_limit:
            if modus == "specified_time":
                # For every second in imposed charging break, insert/shift speed, add charge amount, ...
                for break_second in range(vehicle_data[v]['specified_charge_duration']):
                    
                    # Compute new SOC if charging would occur:
                    if index == 0:
                        SOC_wBreakCharging = (vehicle_data[v]['SOC_t0'] * vehicle_data[v]['Cap_B'] + charge_rate * 1 / 3600) / vehicle_data[v]['Cap_B']
                    else:
                        SOC_wBreakCharging = (drive_cycle.at[(index+break_second)-1, 'SOC_wBreakCharging'] * vehicle_data[v]['Cap_B'] + charge_rate * 1 / 3600) / vehicle_data[v]['Cap_B']
                    delta_SOC = SOC_wBreakCharging - drive_cycle.at[(index+break_second)-1, 'SOC_wBreakCharging']
                    
                    if SOC_wBreakCharging < upper_SOC_limit: 
                        rows_to_copy = drive_cycle.loc[(index+break_second):].copy(deep=True)
                        drive_cycle.loc[max(drive_cycle.index) + 1] = pd.Series() # Add new row

                        # Take care of current timestep: Set all kinetic/params to zero, add SOC, add charge rate
                        drive_cycle.loc[(index+break_second), :] = 0
                        drive_cycle.at[(index+break_second), 'SOC_wBreakCharging'] = SOC_wBreakCharging
                        drive_cycle.at[(index+break_second), 'charging_power_depletion_breaks (kW)'] = charge_rate
                        #drive_cycle.at[(index+break_second), 'SOC'] = drive_cycle.loc[index, 'SOC'] # Carry along the SOC a
                        # Take care of delayed timesteps: Increment SOC for all subsequent timesteps by charge amount
                        rows_to_copy['SOC_wBreakCharging'] += delta_SOC
                        drive_cycle.loc[((index+break_second)+1):] = rows_to_copy.values
                continue
            if modus == "until_full":
                # For every second in imposed charging break, insert/shift speed, add charge amount, ...
                battery_full = False
                second_counter = 0
                while battery_full == False:
                    
                    # Compute new SOC if charging would occur:
                    if index == 0:
                        SOC_wBreakCharging = (vehicle_data[v]['SOC_t0'] * vehicle_data[v]['Cap_B'] + charge_rate * 1 / 3600) / vehicle_data[v]['Cap_B']
                    else:
                        SOC_wBreakCharging = (drive_cycle.at[(index+second_counter)-1, 'SOC_wBreakCharging'] * vehicle_data[v]['Cap_B'] + charge_rate * 1 / 3600) / vehicle_data[v]['Cap_B']
                    delta_SOC = SOC_wBreakCharging - drive_cycle.at[(index+second_counter), 'SOC_wBreakCharging']
                    
                    if SOC_wBreakCharging < upper_SOC_limit: 
                        rows_to_copy = drive_cycle.loc[(index+second_counter):].copy(deep=True)
                        drive_cycle.loc[max(drive_cycle.index) + 1] = pd.Series() # Add new row
                        # Take care of current timestep: Set all kinetic/params to zero, add SOC, add charge rate
                        drive_cycle.loc[(index+second_counter), :] = 0
                        drive_cycle.at[(index+second_counter), 'SOC_wBreakCharging'] = SOC_wBreakCharging
                        drive_cycle.at[(index+second_counter), 'charging_power_depletion_breaks (kW)'] = charge_rate
                        #drive_cycle.at[(index+break_second), 'SOC'] = drive_cycle.loc[index, 'SOC'] # Carry along the SOC a
                        # Take care of delayed timesteps: Increment SOC for all subsequent timesteps by charge amount
                        rows_to_copy['SOC_wBreakCharging'] += delta_SOC
                        drive_cycle.loc[((index+second_counter)+1):] = rows_to_copy.values
                        second_counter += 1
                    else:
                        battery_full = True
                continue
        
        # After end of cycle, charge vehicle to full.
        if (index+1) not in drive_cycle.index.tolist(): # if last second, i.e., last element in dataframe
            print("Last index: ", index)

            battery_full = False
            while battery_full == False: 
                # Compute new SOC if charging would occur:
                SOC_wBreakCharging = (drive_cycle.at[index, 'SOC_wBreakCharging'] * vehicle_data[v]['Cap_B'] + charge_rate * 1 / 3600) / vehicle_data[v]['Cap_B']
                delta_SOC = SOC_wBreakCharging - drive_cycle.at[index, 'SOC_wBreakCharging']
                
                if SOC_wBreakCharging < upper_SOC_limit:
                    drive_cycle.loc[index + 1] = pd.Series() # Add new row
                    # Take care of new timestep: Set all kinetic/params to zero, add SOC, add charge rate
                    drive_cycle.loc[index + 1, :] = 0
                    drive_cycle.at[index + 1, 'SOC_wBreakCharging'] = SOC_wBreakCharging
                    drive_cycle.at[index + 1, 'charging_power_after_cycle (kW)'] = charge_rate
                    index += 1
                    print(SOC_wBreakCharging)
                else:
                    battery_full = True


    drive_cycle['charging_power_total (kW)'] = drive_cycle['charging_power_existing_breaks (kW)'] + drive_cycle['charging_power_depletion_breaks (kW)'] + drive_cycle['charging_power_after_cycle (kW)']
    

    # THIRD STEP: SAVE AND PLOT
    
    drive_cycle.to_csv('C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/' + v + '_Results.csv')

    #plt.subplot(2,1,1)
    #plt.plot(drive_cycle.index, drive_cycle["Speed (mph)"])
    #plt.xlabel('Time')
    #plt.ylabel('Speed')
    #plt.title(v + ': Speed over Time')
    ##plt.show()

    #plt.subplot(2,1,2)
    #plt.plot(drive_cycle.index, drive_cycle["SOC"], label='w/ charging',color='red', linestyle='--')
    #plt.plot(drive_cycle.index, drive_cycle["SOC_wBreakCharging"], label='w/o charging')
    #plt.xlabel('Time')
    #plt.ylabel('SOC')
    #plt.title(v + ': SOC over Time')
    #plt.tight_layout()
    #plt.show()

    # Creating subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    # Plotting in the first subfigure
    ax1.plot(drive_cycle.index, drive_cycle["Speed (mph)"], color='blue')
    ax1.set_ylabel("Speed (mph)")
    ax1.legend()

    # Plotting in the second subfigure with two y-axes
    #ax2.plot(drive_cycle.index, drive_cycle["SOC"], label='w/o charging',color='red', linestyle='--')
    ax2.plot(SOC_copy.index, SOC_copy, label='w/o charging',color='red', linestyle='--')
    #ax2.plot(drive_cycle.index, drive_cycle["SOC"])
    ax2.plot(drive_cycle.index, drive_cycle["SOC_wBreakCharging"], label='w/ charging')
    ax2.set_ylabel('SOC [-]')
    ax2.set_xlabel('Time [s]')
    ax2.legend(loc='lower left')

    ax2_right = ax2.twinx()
    ax2_right.plot(drive_cycle.index, drive_cycle['charging_power_existing_breaks (kW)'], label='Charging power (existing breaks)', color='green', linewidth=0.5)
    ax2_right.plot(drive_cycle.index, drive_cycle['charging_power_depletion_breaks (kW)'], label='Charging power (added breaks)', color='orange', linewidth=0.5)
    ax2_right.plot(drive_cycle.index, drive_cycle['charging_power_after_cycle (kW)'], label='Charging power (after cycle)', color='red', linewidth=0.5)
    ax2_right.set_ylabel('Charging power (kW)')
    ax2_right.legend(loc='upper right')

    # Set limits for right y-axis of lower plot
    ax2_right.set_ylim(-50, 1000)

    # Title for the entire plot
    plt.suptitle(v)

    # Adjust layout for better spacing
    plt.tight_layout()

    # Display the plot
    plt.show()

"""

# SAMPLE VEHICLES AND OBTAIN AVERAGE CHARGING PROFILE

# Set sample size (divide generated profile by this number later)
number_of_sampled_vehicles_per_class = 1000 # 1000

def sample_from_distribution(prob_distribution, size):
    return np.random.choice(len(prob_distribution), size=size, p=prob_distribution)

for v in vehicle_data:
    # Load from csv or continue using objects from above code.
    start_times = pd.read_csv('C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/start_time_distribution.csv', index_col=0)
    drive_cycle = pd.read_csv('C:\\Users/mheyer/Code_Stanford/BRIDGES_for_CA/Data/RawData/transport/' + v + '_Results.csv', index_col=0)

    # Sample from the distribution (start_hours lists start hours [0,23] for each sampled vehicle as list.)
    start_hours = sample_from_distribution(start_times[v], number_of_sampled_vehicles_per_class)
    print(start_hours)

    # Aggregate second-spaced-load profile into hourly load profile (DataFrame with hours as index and a column "charging_power_total (kW)")
    # Consider modeling at an secondly/minutely resolution, e.g., by using a uniform distribution to randomly sample start time within one hour
    drive_cycle["Time (hours)"] = drive_cycle.index
    drive_cycle["Time (hours)"] = drive_cycle["Time (hours)"] // 3600 # Floor division, divide and discard remainer.
    drive_cycle_hourly_load = drive_cycle.groupby("Time (hours)")['charging_power_total (kW)'].agg(lambda x: 1 * x.sum() / 3600).reset_index() # Secondly kW value * 1 and sum => kJ per hour => / 3600 => kW over hour
    drive_cycle_hourly_load.set_index("Time (hours)", inplace=True) 
    print(drive_cycle_hourly_load)

    # Position the representative hourly load profile (drive_cycle_hourly_load) over the day (according to "start_hours" sampled from the distribution) for a number of vehicles (1000).
    vehicle_profiles = pd.DataFrame({"start_hours": start_hours, "charging_power_profile_hourly_perCar (kW)": [drive_cycle_hourly_load['charging_power_total (kW)'].tolist()] * len(start_hours)}) # Has a row for each vehicle, first column is start_hour, second column is charging power hourly load (as list)
    vehicle_profiles = vehicle_profiles.rename_axis('vehicle_id')
    #.reset_index().set_index('vehicle_id') # Give the index axis a label
    def append_zeros(row): # Function that adds zeros to the beginning of the representative drive cycle, therefore shifting it back in the day
        return [0] * row['start_hours'] + row['charging_power_profile_hourly_perCar (kW)']
    vehicle_profiles["charging_power_profile_hourly_perCar (kW)"] = vehicle_profiles.apply(append_zeros, axis=1)
    print(vehicle_profiles)

    # Sum up the charging loads of the 1000 vehicles and divide by 1000 to get average charging load per added vehicle
    max_length = vehicle_profiles["charging_power_profile_hourly_perCar (kW)"].apply(len).max() # Find the maximum length of lists in the column
    padded_lists = vehicle_profiles["charging_power_profile_hourly_perCar (kW)"].apply(lambda x: x + [0] * (max_length - len(x))) # Pad the lists with zeros to make them of equal length
    array_of_lists = np.array(padded_lists.tolist()) # Convert the column of lists to a NumPy array
    charging_power_profile_hourly_total = np.sum(array_of_lists, axis=0) # Sum the arrays element-wise along the columns
    charging_power_profile_hourly_average = charging_power_profile_hourly_total / number_of_sampled_vehicles_per_class # Get average profile per car
    charging_power_profile_hourly_average = pd.DataFrame({'charging_power_profile_hourly_average (kW)': charging_power_profile_hourly_average})
    charging_power_profile_hourly_average = charging_power_profile_hourly_average.rename_axis('Time (hours)')
    #.reset_index().set_index('Time (hours)', inplace=True)
    print(charging_power_profile_hourly_average)

    # Flip hours beyond one day to the start of the day
    charging_power_profile_h_av_day1 = charging_power_profile_hourly_average.head(24)
    charging_power_profile_h_av_day2 = charging_power_profile_hourly_average.tail(len(charging_power_profile_hourly_average)-24).reset_index(drop=True)
    extended_column = charging_power_profile_h_av_day2['charging_power_profile_hourly_average (kW)'].tolist() + [0] * (24-len(charging_power_profile_h_av_day2)) # Pad carryover to 24h
    charging_power_profile_h_av_day2 = pd.DataFrame({"charging_power_profile_hourly_average (kW)": extended_column}) # Pad carryover to 24h
    charging_power_profile_h_av_sum = charging_power_profile_h_av_day1 + charging_power_profile_h_av_day2


    # Plot average charging power load profile
    # Creating subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))

    # Plotting in the first subfigure
    ax1.plot(charging_power_profile_hourly_average.index, charging_power_profile_hourly_average['charging_power_profile_hourly_average (kW)'], color='blue')
    ax1.set_ylabel("Charging load (kW)")
    ax1.legend()

    # Plotting in the second subfigure with two y-axes
    ax2.plot(charging_power_profile_h_av_sum.index, charging_power_profile_h_av_sum['charging_power_profile_hourly_average (kW)'], label='combined', marker='s', markersize=2, markeredgecolor='black', markerfacecolor='black', color='red')
    ax2.plot(charging_power_profile_h_av_day1.index, charging_power_profile_h_av_day1['charging_power_profile_hourly_average (kW)'], label='regular', marker='s', markersize=2, markeredgecolor='black', markerfacecolor='black', linewidth=0.75, color='green', linestyle='--')
    ax2.plot(charging_power_profile_h_av_day2.index, charging_power_profile_h_av_day2['charging_power_profile_hourly_average (kW)'], label='carryover', marker='s', markersize=2, markeredgecolor='black', markerfacecolor='black', linewidth=0.75, color='orange', linestyle='--')
    ax2.set_ylabel("Charging load (kW)")
    ax2.set_xlabel('Time [h]')
    ax2.legend(loc='upper right')

    # Title for the entire plot
    plt.suptitle(v)

    # Adjust layout for better spacing
    plt.tight_layout()

    # Display the plot
    plt.show()










