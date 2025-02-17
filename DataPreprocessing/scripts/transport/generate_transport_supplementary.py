import pandas as pd
import numpy as np



def generate_transport_supplementary(emfac_filepath, output_filepath):

    """ 
    Creates and populates the transport supplementary input file for BRIDGES.
    This file holds different scenarios for VMT projections, vehicle fuel 
    economy projections, and fuel carbon intensity projections.
    
    Parameters:
 
    """
    
    number_of_nodes = 16

    ##########################################
    ###### EMFAC Scenario ####################
    ##########################################

    # This scenario is treated as base case scenario for VMT, CI, and Fuel Economy.

    H2_vehicle_MJ_per_mile = 1.818 # MJ/mi ... 66 mpge (Source: Pathways report) = 66 mi/GGE / (1 kgH2/GGE * 120 MJ/kgH2) and find reciprocal (Kehrwert). Source for 1 kgH2/GGE: http://large.stanford.edu/courses/2016/ph240/kountz2/ 
    NG_vehicle_MJ_per_mile = 3.857 # MJ/mi ... 31 mpge (Source: Pathways report) for CNG: 31 mi/GGE / (2.57 kgCNG/GGE / 2.794e-3 kgCNG/gal * 0.13 MJ/gal) = 31 / 119.578. And find reciprocal (Kehrwert). Sources: https://afdc.energy.gov/fuels/equivalency-methodology (=> 5.66 lbCNG/GGE, 0.0458 lb/scf) and 0.13 MJ/galCNG (CARB, see below).
    H2_vehicle_kgCO2_p_MJ = 0.255      # kgCO2e/MJ For SMR lifecycle. Alternatively 0 tailpipe + electricity emissions (not linked yet). For H2 only tailpipe emissions (in line with other EMFAC values which seem to limit to tailpipe emissions). This is a ToDo.
    NG_vehicle_kgCO2_p_MJ = 0.122  # kgCO2e/MJ 129 kgCO2/mmBTU = 0.122 kg/MJ for NG (entire lifecycle, https://www.nrel.gov/docs/fy16osti/64267.pdf). 

    energy_densities_MJ_per_x = {
        "Gasoline": 115.83, # MJ/gal. Use value for CaRFG (=E10) not CARBOB since most is E10: https://www.eia.gov/todayinenergy/detail.php?id=26092
        "Diesel": 134.47, # MJ/gal
        "Hydrogen": 120, # MJ/kg
        "Natural Gas": 0.13, #  MJ/gal. Source: Table 4 in https://ww2.arb.ca.gov/sites/default/files/barcu/regact/2011/lcfs2011/frooalapp.pdf (2011 rulemaking, since EMFAC has volume for Fuel consumption (1000 gal/year)). 0.98 MJ/scf == :7.48 gal/scf ==> 0.31 MJ/gal. Using value for CNG, since more common for cars, consider switch to LNG for trucks.
    } 
    # Source: Table 4 in https://ww2.arb.ca.gov/our-work/programs/low-carbon-fuel-standard/lcfs-regulation > 2020 rulemaking.
    # On gasoline: CARBOB makes up the petroleum fraction of California reformulated gasoline (CaRFG) before any fuel oxygenate is added; CARFG is essentially 90 percent CARBOB blended with 10 percent ethanol by volume.


    df = pd.read_csv(emfac_filepath, skiprows=7, header=0)

    # Drop all vehicle classes besides LDA, LDT1, and LDT2
    df = df[df['Vehicle Category'].isin(['LDA', 'LDT1', 'LDT2'])]
    # For light-duty: EMFAC-Emissions: LDA, LDT1, LDT2  =>  Use for statewide Population, TotalVMT, CVMT/EVMT, EnergyConsumption for each FUEL, for each year 2019, and 2015-2050 (also available for counties if this resolution is desired)
    #                                  Units: Energy Consumption: kWh/year, Emissions: tons/year, Fuel consumption: 1000 gal/year
    #                                  Has trips per year! - good for calibration!
    #                 EMFAC-Fleet: "P, Passenger Cars","T1, Light-duty trucks (GVWR <6000 lbs, ETW ≤3750 lbs)","T2, Light-duty trucks (GVWR <6000 lbs, ETW 3751–5750 lbs)" => used for vehicle population on COUNTY level

    # Group by Fuel and Year (looses info on "Model Year") and sum other columns.
    df = df.groupby(['Fuel', 'Calendar Year']).sum().reset_index()

    # Compute VMT [mi/veh/year]:
    df["VMT_[mi/yr/vehicle]"] = df["Total VMT"] / df["Population"]

    # Compute vehicle economy, i.e., consumption [MJ/veh/mi]:
    def compute_energy_consumption(row):
        if row["Fuel"] == "Diesel":
            return (row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Diesel"]) / row["CVMT"] # *1000 because number is in "1000 gal/year"
        elif row["Fuel"] == "Gasoline":
            return (row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Gasoline"]) / row["CVMT"]
        elif row["Fuel"] == "Natural Gas":
            return (row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Natural Gas"]) / row["CVMT"]
        elif row["Fuel"] == "Plug-in Hybrid":
            return (row["Energy Consumption"] * 3.6 + row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Gasoline"]) / (row["CVMT"] + row["EVMT"])
        elif row["Fuel"] == "Electricity":
            return (row["Energy Consumption"] * 3.6) / row["EVMT"] # x3.6 converts kWh to MJ
    df['VehicleFuelEconomy_[MJ/vehicle/mile]'] = df.apply(compute_energy_consumption, axis=1)

    # Compute fuel carbon intensity
    def compute_carbon_intensity(row):
        if row["Fuel"] == "Diesel":
            return row["CO2_TOTEX"] * 1000 / (row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Diesel"]) # *1000 because emissions are in ton per year; *1000 because number is in "1000 gal/year"
        elif row["Fuel"] == "Gasoline":
            return row["CO2_TOTEX"] * 1000 / (row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Gasoline"])
        elif row["Fuel"] == "Natural Gas":
            return row["CO2_TOTEX"] * 1000 / (row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Natural Gas"])
        elif row["Fuel"] == "Plug-in Hybrid":
            return row["CO2_TOTEX"] * 1000 / (row["Energy Consumption"] * 3.6 + row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Gasoline"])
        elif row["Fuel"] == "Electricity":
            return row["CO2_TOTEX"] * 1000 / (row["Energy Consumption"] * 3.6) # x3.6 converts kWh to MJ
    df['CarbonIntensity_TransportFuel_[kgCO2e/MJ]'] = df.apply(compute_carbon_intensity, axis=1)

    """
    #########################################
    #### Compute transport fuel prices ######
    #########################################

    # Copy the generated values to FuelCostLookUp_wTransport.csv

    # DATA: 
    # Year:             2020    2021    2022    2023
    # Diesel:           3.377   4.164   6.028   5.357   [$/gal]     EIA. Weekly Retail Gasoline and Diesel Prices (Dollars per Gallon, Including Taxes). https://www.eia.gov/dnav/pet/pet_pri_gnd_dcus_sca_w.htm . No 2 Diesel. Average of weekly averages.
    # Gasoline:         3.050   4.013   5.311   4.773   [$/gal]     EIA. Weekly Retail Gasoline and Diesel Prices (Dollars per Gallon, Including Taxes). https://www.eia.gov/dnav/pet/pet_pri_gnd_dcus_sca_w.htm . Regular Reformulated Gasoline. Average of weekly averages.
    # Hydrogen:         16      16      19      27.5    [$/kg]      S&P Global Commodity Insights. Retail California hydrogen pump price. Average of begin and end of year price. https://www.spglobal.com/commodityinsights/en/market-insights/latest-news/energy-transition/012324-logistical-woes-and-high-pump-prices-stall-california-h2-market-development#:~:text=Platts%20last%20assessed%20California's%20retail,S%26P%20Global%20Commodity%20Insights%20data. 
    # Natural Gas:      2.175   2.233   2.680   2.988   [$/GasolineGallonEquivalent]    DoE. Average retail fuel prices in the US. https://afdc.energy.gov/fuels/prices.html . CNG Retail price (US wide). Average of first days of four quarters.
    # How to project, how is done for other fuels? NG is random values from a year.
    # For projection, check CARBs prediction in case LCFS is established: https://ww2.arb.ca.gov/sites/default/files/2023-09/lcfs_sria_2023_0.pdf (Standardized Regulatory Impact Assessment SRIA for 2023 amendments)

    # To convert $/GGE to $/MJ for CNG: 2.175 $/GGE / (2.57 kgCNG/GGE / 2.794e-3 kgCNG/gal * 0.13 MJ/gal) = 2.175 / 119.578. Sources: https://afdc.energy.gov/fuels/equivalency-methodology (=> 5.66 lbCNG/GGE, 0.0458 lb/scf) and 0.13 MJ/galCNG (CARB, see below).

    # Assuming constant as of 2024. $/MJ.
    Diesel_Transport_Cost = [3.377,4.164,6.028,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357,5.357]
    Diesel_Transport_Cost_USD_per_MJ = [round(x / 134.47, 4) for x in Diesel_Transport_Cost]
    print(Diesel_Transport_Cost_USD_per_MJ)
    Gasoline_Transport_Cost = [3.050,4.013,5.311,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773,4.773]
    Gasoline_Transport_Cost_USD_per_MJ = [round(x / 115.83, 4) for x in Gasoline_Transport_Cost]
    print(Gasoline_Transport_Cost_USD_per_MJ)
    Hydrogen_Transport_Cost = [16,16,19,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5,27.5]
    Hydrogen_Transport_Cost_USD_per_MJ = [round(x / 120, 4) for x in Hydrogen_Transport_Cost]
    print(Hydrogen_Transport_Cost_USD_per_MJ)
    CNG_Transport_Cost = [2.175,2.233,2.680,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988,2.988]
    CNG_Transport_Cost_USD_per_MJ = [round(x / 119.578, 4) for x in CNG_Transport_Cost]
    print(CNG_Transport_Cost_USD_per_MJ)
    def compute_hybrid_cost(row): # *1000 because emissions are in ton per year; *1000 because number is in "1000 gal/year"
        if row["Fuel"] == "Plug-in Hybrid":                                                                                     # In year 2020 pick 0st element
            return (row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Gasoline"] * Gasoline_Transport_Cost_USD_per_MJ[(row["Calendar Year"]-2020)]) / (row["Energy Consumption"] * 3.6 + row["Fuel Consumption"] * 1000 * energy_densities_MJ_per_x["Gasoline"]) # In the numerator there is no summand: "row["Energy Consumption"] * 3.6 * $/MJelec" since we only care about non-electric part. 
    Hybrid_Transport_Cost_USD_per_MJ = df[ (df["Fuel"]=="Plug-in Hybrid") & (df["Calendar Year"]>=2020) ].apply(compute_hybrid_cost, axis=1).round(4).tolist() # Non-electric part of cost per MJ(elec+conventional).
    print(Hybrid_Transport_Cost_USD_per_MJ)
    # Copy these to FuelCostLookUp_wTransport.csv
    """

    ##########################################
    ###### Extending the EMFAC Base Case #####
    ##########################################

    # The EMFAC dataset (as base case dataset) lacks the following data:
    #       (1) VMT, CI, Economy for Plug-in Hybrids in 2000-2009. 
    #               Solution: Take the 2010 value for previous years. However, this data is most likely not required for any analysis.
    #       (2) VMT, CI, Economy for Natural Gas and Hydrogen Vehicles in all years. 
    #               Solution for VMT: Assume same VMT as Gasoline Vehicle for NG and H2 vehicles.
    #               Solution for Fuel Economy: Assume same values as in pathways report for NG and H2 vehicles. Assume them constant over time.
    #               Solution for CI: 0 kg/MJ for H2 (only tailpipe) and 129 kgCO2/mmBTU = 0.122 kg/MJ for NG (entire lifecycle, https://www.nrel.gov/docs/fy16osti/64267.pdf). For H2 only tailpipe emissions (in line with other EMFAC values which seem to limit to tailpipe emissions). This is a ToDo.

    df = df[['Fuel', 'Calendar Year', 'VMT_[mi/yr/vehicle]', 'VehicleFuelEconomy_[MJ/vehicle/mile]', 'CarbonIntensity_TransportFuel_[kgCO2e/MJ]']]

    # See (1)
    for year in range(2000,2010): # i.e., 2000-2009
        added_rows = {  "Fuel": "Plug-in Hybrid",
                        "Calendar Year": year,
                        "VMT_[mi/yr/vehicle]":                       df.loc[((df['Fuel'] == 'Plug-in Hybrid') & (df['Calendar Year'] == 2010)), "VMT_[mi/yr/vehicle]"].values[0],
                        'VehicleFuelEconomy_[MJ/vehicle/mile]':      df.loc[((df['Fuel'] == 'Plug-in Hybrid') & (df['Calendar Year'] == 2010)), 'VehicleFuelEconomy_[MJ/vehicle/mile]'].values[0],
                        'CarbonIntensity_TransportFuel_[kgCO2e/MJ]': df.loc[((df['Fuel'] == 'Plug-in Hybrid') & (df['Calendar Year'] == 2010)), 'CarbonIntensity_TransportFuel_[kgCO2e/MJ]'].values[0]}
        df = df.append(added_rows, ignore_index=True)

    # See (2) (Hydrogen)
    for year in range(2000,2051): # i.e., 2000-2050
        added_rows = {  "Fuel": "Hydrogen",
                        "Calendar Year": year,
                        "VMT_[mi/yr/vehicle]":                       df.loc[((df['Fuel'] == 'Gasoline') & (df['Calendar Year'] == year)), "VMT_[mi/yr/vehicle]"].values[0],
                        'VehicleFuelEconomy_[MJ/vehicle/mile]':      H2_vehicle_MJ_per_mile,
                        'CarbonIntensity_TransportFuel_[kgCO2e/MJ]': H2_vehicle_kgCO2_p_MJ}
        df = df.append(added_rows, ignore_index=True)

    # See (2) (Natural Gas)
    for year in range(2000,2051): # i.e., 2000-2050
        added_rows = {  "Fuel": "Natural Gas",
                        "Calendar Year": year,
                        "VMT_[mi/yr/vehicle]":                       df.loc[((df['Fuel'] == 'Gasoline') & (df['Calendar Year'] == year)), "VMT_[mi/yr/vehicle]"].values[0],
                        'VehicleFuelEconomy_[MJ/vehicle/mile]':      NG_vehicle_MJ_per_mile,
                        'CarbonIntensity_TransportFuel_[kgCO2e/MJ]': NG_vehicle_kgCO2_p_MJ}
        df = df.append(added_rows, ignore_index=True)


    ##########################################
    ###### LCFS Scenario #####################
    ##########################################

    # Provide alternative CI values compared to EMFAC. I assume that the LCFS standards only effect gasoline and diesel. For the other fuels the EMFAC base case values are used for CI values.

    index = [yr for yr in range(2011,2051)]
    #print(index)
    LCFS_scenario_data = {
        #                                      2011                                2015                                         2020                                         2025                                         2030                                         2035                                         2040                                         2045                                         2050
        "LCFS_current_gasoline_[kgCO2e/MJ]":  [0.09561, 0.09537, 0.09796, 0.09796, 0.09796, 0.09650, 0.09502, 0.09355, 0.09323, 0.09198, 0.09074, 0.08950, 0.08825, 0.08701, 0.08577, 0.08452, 0.08328, 0.08204, 0.08080, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955, 0.07955],
        "LCFS_current_diesel_[kgCO2e/MJ]":    [0.09447, 0.09424, 0.09705, 0.09705, 0.09705, 0.09997, 0.09844, 0.09691, 0.09417, 0.09292, 0.09166, 0.09041, 0.08915, 0.08789, 0.08664, 0.08538, 0.08413, 0.08287, 0.08162, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036, 0.08036],
        # Proposed: -30% by 2030, -90% by 2045 (2010 baseline), one-time -5% step in 2025. (Available as absolut numbers in regulation text and percentage reduction in FSOR.)
        #           Percent reduction (2010 baseline):                                                                                                              12.5%    18.75%                                       30%                                          52.5%                                        75%                                          90%      90%
        "LCFS_proposed_gasoline_[kgCO2e/MJ]": [0.09561, 0.09537, 0.09796, 0.09796, 0.09796, 0.09650, 0.09502, 0.09355, 0.09323, 0.09198, 0.09074, 0.08950, 0.08825, 0.08701, 0.08055, 0.07832, 0.07609, 0.07386, 0.07163, 0.06940, 0.06494, 0.06048, 0.05602, 0.05155, 0.04709, 0.04263, 0.03817, 0.03371, 0.02924, 0.02478, 0.02181, 0.01883, 0.01586, 0.01288, 0.00991, 0.00991, 0.00991, 0.00991, 0.00991, 0.00991],
        "LCFS_proposed_diesel_[kgCO2e/MJ]":   [0.09447, 0.09424, 0.09705, 0.09705, 0.09705, 0.09997, 0.09844, 0.09691, 0.09417, 0.09292, 0.09166, 0.09041, 0.08915, 0.08789, 0.08593, 0.08355, 0.08117, 0.07879, 0.07641, 0.07403, 0.06927, 0.06451, 0.05975, 0.05499, 0.05023, 0.04547, 0.04071, 0.03595, 0.03119, 0.02644, 0.02326, 0.02009, 0.01692, 0.01374, 0.01057, 0.01057, 0.01057, 0.01057, 0.01057, 0.01057],
        # Alternative 1: -28% by 2030, -90% by 2045 (2010 baseline), one-time -3% step in 2025. (Available as percentage reduction in FSOR.)                        12.4%
        "LCFS_alternative1_[percent_reduct]": [                                                                                                                     0.1240 , 0.1680 , 0.1900 , 0.2130 , 0.2350 , 0.2580 , 0.2800 , 0.3270 , 0.3740 , 0.4210 , 0.4680 , 0.5150 , 0.5620 , 0.6090 , 0.6560 , 0.7030 , 0.7500 , 0.7800 , 0.8100 , 0.8400 , 0.8700 , 0.9000 , 0.9000 , 0.9000 , 0.9000 , 0.9000 , 0.9000],
        # Alternative 2: -35% by 2030, -90% by 2045 (2010 baseline), one-time -5% step in 2025. (Available as percentage reduction in FSOR.)                        12.4%
        "LCFS_alternative2_[percent_reduct]": [                                                                                                                     0.1240 , 0.1860 , 0.2190 , 0.2520 , 0.2850 , 0.3170 , 0.3500 , 0.3900 , 0.4300 , 0.4700 , 0.5100 , 0.5500 , 0.5900 , 0.6300 , 0.6700 , 0.7100 , 0.7500 , 0.7800 , 0.8100 , 0.8400 , 0.8700 , 0.9000 , 0.9000 , 0.9000 , 0.9000 , 0.9000 , 0.9000],
        }

    # From the absolut CI numbers in the regulation text and percentage reduction in the FSOR we can infer that the 2010 baseline is 0.09944 kgCO2e/MJ (gasoline) and 0.10044 kgCO2e/MJ diesel.
    LCFS_scenario_data["LCFS_alternative1_gasoline_[kgCO2e/MJ]"] = LCFS_scenario_data["LCFS_proposed_gasoline_[kgCO2e/MJ]"][:13] + [0.09944 * (1-a) for a in LCFS_scenario_data["LCFS_alternative1_[percent_reduct]"]]
    LCFS_scenario_data["LCFS_alternative1_diesel_[kgCO2e/MJ]"]   = LCFS_scenario_data["LCFS_proposed_diesel_[kgCO2e/MJ]"][:13] + [0.10044 * (1-a) for a in LCFS_scenario_data["LCFS_alternative1_[percent_reduct]"]]
    LCFS_scenario_data["LCFS_alternative2_gasoline_[kgCO2e/MJ]"] = LCFS_scenario_data["LCFS_proposed_gasoline_[kgCO2e/MJ]"][:13] + [0.09944 * (1-a) for a in LCFS_scenario_data["LCFS_alternative2_[percent_reduct]"]]
    LCFS_scenario_data["LCFS_alternative2_diesel_[kgCO2e/MJ]"]   = LCFS_scenario_data["LCFS_proposed_diesel_[kgCO2e/MJ]"][:13] + [0.10044 * (1-a) for a in LCFS_scenario_data["LCFS_alternative2_[percent_reduct]"]]

    # Delete the percentage rows
    LCFS_scenario_data.pop("LCFS_alternative1_[percent_reduct]")
    LCFS_scenario_data.pop("LCFS_alternative2_[percent_reduct]")

    # Add CI values for Plug-in Hybrid, Electricity, Natural Gas and Hydrogen.
    # (For Elec, NG and H2 use the EMFAC values).
    df_aux_gasoline = df[df["Fuel"] == "Gasoline"]
    df_aux_diesel = df[df["Fuel"] == "Diesel"]
    df_aux_hybrid = df[df["Fuel"] == "Plug-in Hybrid"]
    df_aux_natgas = df[df["Fuel"] == "Natural Gas"]
    df_aux_hydrogen = df[df["Fuel"] == "Hydrogen"]
    df_aux_elec = df[df["Fuel"] == "Electricity"]
    
    LCFS_scenario_data["LCFS_current_elec_[kgCO2e/MJ]"] =           [df_aux_elec.loc[(df_aux_elec["Calendar Year"]==year), "CarbonIntensity_TransportFuel_[kgCO2e/MJ]"].values[0] for year in index]
    LCFS_scenario_data["LCFS_current_natgas_[kgCO2e/MJ]"] =         [df_aux_natgas.loc[(df_aux_natgas["Calendar Year"]==year), "CarbonIntensity_TransportFuel_[kgCO2e/MJ]"].values[0] for year in index]
    LCFS_scenario_data["LCFS_current_hydrogen_[kgCO2e/MJ]"] =       [df_aux_hydrogen.loc[(df_aux_hydrogen["Calendar Year"]==year), "CarbonIntensity_TransportFuel_[kgCO2e/MJ]"].values[0] for year in index]
    LCFS_scenario_data["LCFS_proposed_elec_[kgCO2e/MJ]"] =          LCFS_scenario_data["LCFS_current_elec_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_proposed_natgas_[kgCO2e/MJ]"] =        LCFS_scenario_data["LCFS_current_natgas_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_proposed_hydrogen_[kgCO2e/MJ]"] =      LCFS_scenario_data["LCFS_current_hydrogen_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_alternative1_elec_[kgCO2e/MJ]"] =      LCFS_scenario_data["LCFS_current_elec_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_alternative1_natgas_[kgCO2e/MJ]"] =    LCFS_scenario_data["LCFS_current_natgas_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_alternative1_hydrogen_[kgCO2e/MJ]"] =  LCFS_scenario_data["LCFS_current_hydrogen_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_alternative2_elec_[kgCO2e/MJ]"] =      LCFS_scenario_data["LCFS_current_elec_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_alternative2_natgas_[kgCO2e/MJ]"] =    LCFS_scenario_data["LCFS_current_natgas_[kgCO2e/MJ]"]
    LCFS_scenario_data["LCFS_alternative2_hydrogen_[kgCO2e/MJ]"] =  LCFS_scenario_data["LCFS_current_hydrogen_[kgCO2e/MJ]"]

    # (For Hybrid: Assumption: CI_EMFAC_hybrid / CI_EMFAC_gasoline = const, i.e., ration of electric and gas miles unchanged. => CI_LCFS_hybrid = const * CI_LCFS_gasoline)
    ratiolist = [x / y for x,y in [(df_aux_hybrid.loc[(df_aux_hybrid["Calendar Year"]==year), "CarbonIntensity_TransportFuel_[kgCO2e/MJ]"].values[0],
                                    df_aux_gasoline.loc[(df_aux_gasoline["Calendar Year"]==year), "CarbonIntensity_TransportFuel_[kgCO2e/MJ]"].values[0]) for year in range(2011,2051)]] # i.e., 2011-2050
    LCFS_scenario_data["LCFS_current_hybrid_[kgCO2e/MJ]"] =         [x * y for x, y in zip(ratiolist, LCFS_scenario_data["LCFS_current_gasoline_[kgCO2e/MJ]"])] # elementwise multiplication of lists
    LCFS_scenario_data["LCFS_proposed_hybrid_[kgCO2e/MJ]"] =        [x * y for x, y in zip(ratiolist, LCFS_scenario_data["LCFS_proposed_gasoline_[kgCO2e/MJ]"])]
    LCFS_scenario_data["LCFS_alternative1_hybrid_[kgCO2e/MJ]"] =    [x * y for x, y in zip(ratiolist, LCFS_scenario_data["LCFS_alternative1_gasoline_[kgCO2e/MJ]"])]
    LCFS_scenario_data["LCFS_alternative2_hybrid_[kgCO2e/MJ]"] =   [x * y for x, y in zip(ratiolist, LCFS_scenario_data["LCFS_alternative2_gasoline_[kgCO2e/MJ]"])]

    df_LCFS_scenario_data = pd.DataFrame(LCFS_scenario_data, index)



    ##########################################
    ###### VMT Scenario ######################
    ##########################################

    # Scoping Plan - BAU scenario: -4% (2019 baseline) to 2045
    # Scoping Plan - Proposed scenario: -25% (2019) by 2030, -30% (2019) by 2045
    # Assumption: to be observed over all vehicle classes.

    # CECs California human population prediction as of Jan 2024. 2022 is last historic year. Obtained xlsx from Dimitri: TN254253_20240131T141810_CED 2023 Baseline Forecast - Total State.xlsx
    CA_population = {
        2019: 39529566,
        2030: 39430871,
        2045: 40259016, # data only goes to 2040. So continued with last growth rate of 2039-2040: 0.000759657 [-]
    }

    #df_VMT_scenarios = df[(df["Calendar Year"] == 2019) | (df["Calendar Year"] == 2030) | (df["Calendar Year"] == 2045)]
    df_VMT_scenarios = df[df["Calendar Year"] == 2019]
    def multiply_population(row):
        return row["VMT_[mi/yr/vehicle]"] / CA_population.get(row["Calendar Year"], 0)
    df_VMT_scenarios["VMT_[mi/yr/vehicle/person]"] = df_VMT_scenarios.apply(multiply_population, axis=1)
    df_VMT_scenarios["VMT_BAU_2045_[mi/yr/vehicle]"] = df_VMT_scenarios["VMT_[mi/yr/vehicle/person]"]*0.96*CA_population[2045] # VMT value for each fuel in 2045: 4% reduction compared to 2019 base line
    df_VMT_scenarios["VMT_ScopingPlan_2030_[mi/yr/vehicle]"] = df_VMT_scenarios["VMT_[mi/yr/vehicle/person]"]*0.75*CA_population[2030] # VMT value for each fuel in 2030: 25% reduction compared to 2019 base line
    df_VMT_scenarios["VMT_ScopingPlan_2045_[mi/yr/vehicle]"] = df_VMT_scenarios["VMT_[mi/yr/vehicle/person]"]*0.70*CA_population[2045] # VMT value for each fuel in 2045: 30% reduction compared to 2019 base line
    df_VMT_scenarios.set_index('Fuel', inplace=True)


    index_VMT = [yr for yr in range(2000,2051)] # i.e., 2000-2050

    df_aux_gasoline = df[df["Fuel"] == "Gasoline"]
    df_aux_diesel = df[df["Fuel"] == "Diesel"]
    df_aux_hybrid = df[df["Fuel"] == "Plug-in Hybrid"]
    df_aux_natgas = df[df["Fuel"] == "Natural Gas"]
    df_aux_hydrogen = df[df["Fuel"] == "Hydrogen"]
    df_aux_elec = df[df["Fuel"] == "Electricity"]

    VMT_scenario_data = {
        "VMT_gasoline_BAU_[mi/yr/vehicle]":         [df_aux_gasoline.loc[(df_aux_gasoline["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]] + np.linspace(df_VMT_scenarios.loc["Gasoline","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Gasoline","VMT_BAU_2045_[mi/yr/vehicle]"], num=27).tolist() + [df_VMT_scenarios.loc["Gasoline","VMT_BAU_2045_[mi/yr/vehicle]"]]*5,
        "VMT_gasoline_ScopingPlan_[mi/yr/vehicle]": [df_aux_gasoline.loc[(df_aux_gasoline["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]] + np.linspace(df_VMT_scenarios.loc["Gasoline","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Gasoline","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], num=12).tolist() + np.linspace(df_VMT_scenarios.loc["Gasoline","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Gasoline","VMT_ScopingPlan_2045_[mi/yr/vehicle]"], num=16).tolist()[1:] + [df_VMT_scenarios.loc["Gasoline","VMT_ScopingPlan_2045_[mi/yr/vehicle]"]]*5,
        "VMT_diesel_BAU_[mi/yr/vehicle]":           [df_aux_diesel.loc[(df_aux_diesel["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Diesel","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Diesel","VMT_BAU_2045_[mi/yr/vehicle]"], num=27).tolist() + [df_VMT_scenarios.loc["Diesel","VMT_BAU_2045_[mi/yr/vehicle]"]]*5,
        "VMT_diesel_ScopingPlan_[mi/yr/vehicle]":   [df_aux_diesel.loc[(df_aux_diesel["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Diesel","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Diesel","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], num=12).tolist() + np.linspace(df_VMT_scenarios.loc["Diesel","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Diesel","VMT_ScopingPlan_2045_[mi/yr/vehicle]"], num=16).tolist()[1:] + [df_VMT_scenarios.loc["Diesel","VMT_ScopingPlan_2045_[mi/yr/vehicle]"]]*5,
        "VMT_hybrid_BAU_[mi/yr/vehicle]":           [df_aux_hybrid.loc[(df_aux_hybrid["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Plug-in Hybrid","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Plug-in Hybrid","VMT_BAU_2045_[mi/yr/vehicle]"], num=27).tolist() + [df_VMT_scenarios.loc["Plug-in Hybrid","VMT_BAU_2045_[mi/yr/vehicle]"]]*5,
        "VMT_hybrid_ScopingPlan_[mi/yr/vehicle]":   [df_aux_hybrid.loc[(df_aux_hybrid["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Plug-in Hybrid","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Plug-in Hybrid","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], num=12).tolist() + np.linspace(df_VMT_scenarios.loc["Plug-in Hybrid","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Plug-in Hybrid","VMT_ScopingPlan_2045_[mi/yr/vehicle]"], num=16).tolist()[1:] + [df_VMT_scenarios.loc["Plug-in Hybrid","VMT_ScopingPlan_2045_[mi/yr/vehicle]"]]*5,
        "VMT_natgas_BAU_[mi/yr/vehicle]":           [df_aux_natgas.loc[(df_aux_natgas["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Natural Gas","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Natural Gas","VMT_BAU_2045_[mi/yr/vehicle]"], num=27).tolist() + [df_VMT_scenarios.loc["Natural Gas","VMT_BAU_2045_[mi/yr/vehicle]"]]*5,
        "VMT_natgas_ScopingPlan_[mi/yr/vehicle]":   [df_aux_natgas.loc[(df_aux_natgas["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Natural Gas","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Natural Gas","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], num=12).tolist() + np.linspace(df_VMT_scenarios.loc["Natural Gas","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Natural Gas","VMT_ScopingPlan_2045_[mi/yr/vehicle]"], num=16).tolist()[1:] + [df_VMT_scenarios.loc["Natural Gas","VMT_ScopingPlan_2045_[mi/yr/vehicle]"]]*5,
        "VMT_hydrogen_BAU_[mi/yr/vehicle]":         [df_aux_hydrogen.loc[(df_aux_hydrogen["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Hydrogen","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Hydrogen","VMT_BAU_2045_[mi/yr/vehicle]"], num=27).tolist() + [df_VMT_scenarios.loc["Hydrogen","VMT_BAU_2045_[mi/yr/vehicle]"]]*5,
        "VMT_hydrogen_ScopingPlan_[mi/yr/vehicle]": [df_aux_hydrogen.loc[(df_aux_hydrogen["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]     + np.linspace(df_VMT_scenarios.loc["Hydrogen","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Hydrogen","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], num=12).tolist() + np.linspace(df_VMT_scenarios.loc["Hydrogen","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Hydrogen","VMT_ScopingPlan_2045_[mi/yr/vehicle]"], num=16).tolist()[1:] + [df_VMT_scenarios.loc["Hydrogen","VMT_ScopingPlan_2045_[mi/yr/vehicle]"]]*5,
        "VMT_elec_BAU_[mi/yr/vehicle]":             [df_aux_elec.loc[(df_aux_elec["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]         + np.linspace(df_VMT_scenarios.loc["Electricity","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Electricity","VMT_BAU_2045_[mi/yr/vehicle]"], num=27).tolist() + [df_VMT_scenarios.loc["Electricity","VMT_BAU_2045_[mi/yr/vehicle]"]]*5,
        "VMT_elec_ScopingPlan_[mi/yr/vehicle]":     [df_aux_elec.loc[(df_aux_elec["Calendar Year"]==year), "VMT_[mi/yr/vehicle]"].values[0] for year in index_VMT[:19]]         + np.linspace(df_VMT_scenarios.loc["Electricity","VMT_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Electricity","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], num=12).tolist() + np.linspace(df_VMT_scenarios.loc["Electricity","VMT_ScopingPlan_2030_[mi/yr/vehicle]"], df_VMT_scenarios.loc["Electricity","VMT_ScopingPlan_2045_[mi/yr/vehicle]"], num=16).tolist()[1:] + [df_VMT_scenarios.loc["Electricity","VMT_ScopingPlan_2045_[mi/yr/vehicle]"]]*5,
        }
    
    df_VMT_scenario_data = pd.DataFrame(VMT_scenario_data, index_VMT)

    #print(df)
    #print(df_LCFS_scenario_data)
    #print(df_VMT_scenario_data)




    # COMBINE SCENARIOS IN FILE FOR EXPORT TO CSV:
    columns = ["Data_type", "Scenario", "Node", "Converter_name"] + list(range(2000, 2051))
    df_transport_supplementary = pd.DataFrame(columns=columns)

    # VMT: EMFAC
    for converter in df["Fuel"].unique().tolist(): # Include all fuels that exist in "EMFAC Base Case"
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "VMT_[mi/yr/vehicle]", 
                "Scenario": "EMFAC", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(row["Calendar Year"], row["VMT_[mi/yr/vehicle]"]) for index, row in df[df["Fuel"]==converter].iterrows()]}
            meta_dict.update(data_dict) # concat meta and data dictionary
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)

    # VMT: Scoping_Plan_BAU
    for converter, column_name in [("Diesel", "VMT_diesel_BAU_[mi/yr/vehicle]"), ("Gasoline", "VMT_gasoline_BAU_[mi/yr/vehicle]"), ("Plug-in Hybrid", "VMT_hybrid_BAU_[mi/yr/vehicle]"), ("Natural Gas", "VMT_natgas_BAU_[mi/yr/vehicle]"), ("Hydrogen", "VMT_hydrogen_BAU_[mi/yr/vehicle]"), ("Electricity", "VMT_elec_BAU_[mi/yr/vehicle]")]: # Include all fuels that exist in Scoping Plan Scenario
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "VMT_[mi/yr/vehicle]", 
                "Scenario": "Scoping_Plan_BAU", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(index, row[column_name]) for index, row in df_VMT_scenario_data.iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)
    
    # VMT: Scoping_Plan_Proposed
    for converter, column_name in [("Diesel", "VMT_diesel_ScopingPlan_[mi/yr/vehicle]"), ("Gasoline", "VMT_gasoline_ScopingPlan_[mi/yr/vehicle]"), ("Plug-in Hybrid", "VMT_hybrid_ScopingPlan_[mi/yr/vehicle]"), ("Natural Gas", "VMT_natgas_ScopingPlan_[mi/yr/vehicle]"), ("Hydrogen", "VMT_hydrogen_ScopingPlan_[mi/yr/vehicle]"), ("Electricity", "VMT_elec_ScopingPlan_[mi/yr/vehicle]")]: # Include all fuels that exist in Scoping Plan Scenario
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "VMT_[mi/yr/vehicle]", 
                "Scenario": "Scoping_Plan_Proposed", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(index, row[column_name]) for index, row in df_VMT_scenario_data.iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)

    # Carbon Intensity: EMFAC
    for converter in df["Fuel"].unique().tolist(): # Include all fuels that exist in EMFAC
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "CarbonIntensity_TransportFuel_[kgCO2e/MJ]", 
                "Scenario": "EMFAC", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(row["Calendar Year"], row["CarbonIntensity_TransportFuel_[kgCO2e/MJ]"]) for index, row in df[df["Fuel"]==converter].iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)

    # Carbon Intensity: LCFS_current
    for converter, column_name in [("Diesel", "LCFS_current_diesel_[kgCO2e/MJ]"), ("Gasoline", "LCFS_current_gasoline_[kgCO2e/MJ]"), ("Electricity", "LCFS_current_elec_[kgCO2e/MJ]"), ("Plug-in Hybrid", "LCFS_current_hybrid_[kgCO2e/MJ]"), ("Natural Gas", "LCFS_current_natgas_[kgCO2e/MJ]"), ("Hydrogen", "LCFS_current_hydrogen_[kgCO2e/MJ]")]:
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "CarbonIntensity_TransportFuel_[kgCO2e/MJ]", 
                "Scenario": "LCFS_current", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(index, row[column_name]) for index, row in df_LCFS_scenario_data.iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)
    
    # Carbon Intensity: LCFS_proposed
    for converter, column_name in [("Diesel", "LCFS_proposed_diesel_[kgCO2e/MJ]"), ("Gasoline", "LCFS_proposed_gasoline_[kgCO2e/MJ]"), ("Electricity", "LCFS_proposed_elec_[kgCO2e/MJ]"), ("Plug-in Hybrid", "LCFS_proposed_hybrid_[kgCO2e/MJ]"), ("Natural Gas", "LCFS_proposed_natgas_[kgCO2e/MJ]"), ("Hydrogen", "LCFS_proposed_hydrogen_[kgCO2e/MJ]")]: # Include all fuels that exist in LCFS regulation
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "CarbonIntensity_TransportFuel_[kgCO2e/MJ]", 
                "Scenario": "LCFS_proposed", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(index, row[column_name]) for index, row in df_LCFS_scenario_data.iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)
    
    # Carbon Intensity: LCFS_alternative1
    for converter, column_name in [("Diesel", "LCFS_alternative1_diesel_[kgCO2e/MJ]"), ("Gasoline", "LCFS_alternative1_gasoline_[kgCO2e/MJ]"), ("Electricity", "LCFS_alternative1_elec_[kgCO2e/MJ]"), ("Plug-in Hybrid", "LCFS_alternative1_hybrid_[kgCO2e/MJ]"), ("Natural Gas", "LCFS_alternative1_natgas_[kgCO2e/MJ]"), ("Hydrogen", "LCFS_alternative1_hydrogen_[kgCO2e/MJ]")]: # Include all fuels that exist in LCFS regulation
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "CarbonIntensity_TransportFuel_[kgCO2e/MJ]", 
                "Scenario": "LCFS_alternative1", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(index, row[column_name]) for index, row in df_LCFS_scenario_data.iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)

    # Carbon Intensity: LCFS_alternative2
    for converter, column_name in [("Diesel", "LCFS_alternative2_diesel_[kgCO2e/MJ]"), ("Gasoline", "LCFS_alternative2_gasoline_[kgCO2e/MJ]"), ("Electricity", "LCFS_alternative2_elec_[kgCO2e/MJ]"), ("Plug-in Hybrid", "LCFS_alternative2_hybrid_[kgCO2e/MJ]"), ("Natural Gas", "LCFS_alternative2_natgas_[kgCO2e/MJ]"), ("Hydrogen", "LCFS_alternative2_hydrogen_[kgCO2e/MJ]")]: # Include all fuels that exist in LCFS regulation
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "CarbonIntensity_TransportFuel_[kgCO2e/MJ]", 
                "Scenario": "LCFS_alternative2", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(index, row[column_name]) for index, row in df_LCFS_scenario_data.iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)

    # Fuel economy: EMFAC
    for converter in df["Fuel"].unique().tolist(): # Include all fuels that exist in EMFAC
        for node in range(1,number_of_nodes+1):
            meta_dict = {
                "Data_type": "VehicleFuelEconomy_[MJ/vehicle/mile]", 
                "Scenario": "EMFAC", 
                "Node": node, 
                "Converter_name": converter,
            }
            data_dict = {key: value for key, value in [(row["Calendar Year"], row["VehicleFuelEconomy_[MJ/vehicle/mile]"]) for index, row in df[df["Fuel"]==converter].iterrows()]}
            meta_dict.update(data_dict)
            df_transport_supplementary = pd.concat([df_transport_supplementary, pd.DataFrame([meta_dict])], ignore_index=True)

    #print(df_transport_supplementary)
    df_transport_supplementary.to_csv(output_filepath, index=False)
    # But still has some empty fields and missing values (i.e., values of CI in LCFS scenarios before 2011).


    return df, df_LCFS_scenario_data, df_VMT_scenario_data

    """
    # df = df[df['Calendar Year'].isin([2050])]
    # print(df.iloc[30,26:32])
    # print(df["Vehicle Category"].unique())

    number_of_nodes = 16
    considered_converter_types = {
        'Veh_Road_LDA_Gasoline':   {'VMT_[mi/yr/vehicle]': 18836, 'VehicleFuelEconomy_[MJ/mile]': 0, 'CarbonIntensity_TransportFuel_[kgCO2/MJ]': 0.257},
        'Veh_Road_LDA_Diesel':     {'VMT [mi/yr/vehicle]': 18836, 'VehicleFuelEconomy_[MJ/mile]': 0, 'CarbonIntensity_TransportFuel_[kgCO2/MJ]': 0.206},
        'Veh_Road_LDA_HybridElec': {'VMT [mi/yr/vehicle]': 18836, 'VehicleFuelEconomy_[MJ/mile]': 0, 'CarbonIntensity_TransportFuel_[kgCO2/MJ]': 0.134},
        'Veh_Road_LDA_NaturalGas': {'VMT [mi/yr/vehicle]': 18836, 'VehicleFuelEconomy_[MJ/mile]': 0, 'CarbonIntensity_TransportFuel_[kgCO2/MJ]': 0.390},
        'Veh_Road_LDA_HydrogenFC': {'VMT [mi/yr/vehicle]': 18836, 'VehicleFuelEconomy_[MJ/mile]': 0, 'CarbonIntensity_TransportFuel_[kgCO2/MJ]': 0.255}, #0 # Not in EMFAC database
        'Veh_Road_LDA_Elec':       {'VMT [mi/yr/vehicle]': 18836, 'VehicleFuelEconomy_[MJ/mile]': 0, 'CarbonIntensity_TransportFuel_[kgCO2/MJ]': 0}
        }

    columns = ['Data_Type', 'Scenario', 'Node', 'Converter_Name'] + [str(year) for year in range(2020, 2051)]
    df = pd.DataFrame(columns=columns)

    # ADD VMT PROJECTIONS


    # ADD VECHICLE FUEL ECONOMY PROJECTIONS

    # ADD FUEL CARBON INTENSITY PROJECTIONS

    """


def plot_curves(df, df_LCFS_scenario_data, df_VMT_scenario_data):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np

    mpl.rcParams['lines.linewidth'] = 1.25


    df.set_index('Calendar Year', inplace=True)
    print(df)

    y1_gasoline = df[df['Fuel'].isin(["Gasoline"])].sort_index()["VehicleFuelEconomy_[MJ/vehicle/mile]"]
    y2_gasoline = df[df['Fuel'].isin(["Gasoline"])].sort_index()["VMT_[mi/yr/vehicle]"]
    y2_gasoline_BAU =              df_VMT_scenario_data["VMT_gasoline_BAU_[mi/yr/vehicle]"]
    y2_gasoline_ScopingPlan =      df_VMT_scenario_data["VMT_gasoline_ScopingPlan_[mi/yr/vehicle]"]
    y3_gasoline = df[df['Fuel'].isin(["Gasoline"])].sort_index()["CarbonIntensity_TransportFuel_[kgCO2e/MJ]"]
    y3_gasoline_LCFScurrent =      df_LCFS_scenario_data['LCFS_current_gasoline_[kgCO2e/MJ]']
    y3_gasoline_LCFSproposed =     df_LCFS_scenario_data['LCFS_proposed_gasoline_[kgCO2e/MJ]']
    y3_gasoline_LCFSalternative1 = df_LCFS_scenario_data['LCFS_alternative1_gasoline_[kgCO2e/MJ]']
    y3_gasoline_LCFSalternative2 = df_LCFS_scenario_data['LCFS_alternative2_gasoline_[kgCO2e/MJ]']
    y4_gasoline = y1_gasoline.mul(y2_gasoline).mul(y3_gasoline)

    y1_diesel = df[df['Fuel'].isin(["Diesel"])].sort_index()["VehicleFuelEconomy_[MJ/vehicle/mile]"]
    y2_diesel = df[df['Fuel'].isin(["Diesel"])].sort_index()["VMT_[mi/yr/vehicle]"]
    y2_diesel_BAU =              df_VMT_scenario_data["VMT_diesel_BAU_[mi/yr/vehicle]"]
    y2_diesel_ScopingPlan =      df_VMT_scenario_data["VMT_diesel_ScopingPlan_[mi/yr/vehicle]"] 
    y3_diesel = df[df['Fuel'].isin(["Diesel"])].sort_index()["CarbonIntensity_TransportFuel_[kgCO2e/MJ]"]
    y3_diesel_LCFScurrent =      df_LCFS_scenario_data['LCFS_current_diesel_[kgCO2e/MJ]']
    y3_diesel_LCFSproposed =     df_LCFS_scenario_data['LCFS_proposed_diesel_[kgCO2e/MJ]']
    y3_diesel_LCFSalternative1 = df_LCFS_scenario_data['LCFS_alternative1_diesel_[kgCO2e/MJ]']
    y3_diesel_LCFSalternative2 = df_LCFS_scenario_data['LCFS_alternative2_diesel_[kgCO2e/MJ]']
    y4_diesel = y1_diesel.mul(y2_diesel).mul(y3_diesel)

    y1_hybrid = df[df['Fuel'].isin(["Plug-in Hybrid"])].sort_index()["VehicleFuelEconomy_[MJ/vehicle/mile]"]
    y2_hybrid = df[df['Fuel'].isin(["Plug-in Hybrid"])].sort_index()["VMT_[mi/yr/vehicle]"]
    y2_hybrid_BAU =              df_VMT_scenario_data["VMT_hybrid_BAU_[mi/yr/vehicle]"]
    y2_hybrid_ScopingPlan =      df_VMT_scenario_data["VMT_hybrid_ScopingPlan_[mi/yr/vehicle]"] 
    y3_hybrid = df[df['Fuel'].isin(["Plug-in Hybrid"])].sort_index()["CarbonIntensity_TransportFuel_[kgCO2e/MJ]"]
    y3_hybrid_LCFScurrent =      df_LCFS_scenario_data['LCFS_current_hybrid_[kgCO2e/MJ]']
    y3_hybrid_LCFSproposed =     df_LCFS_scenario_data['LCFS_proposed_hybrid_[kgCO2e/MJ]']
    y3_hybrid_LCFSalternative1 = df_LCFS_scenario_data['LCFS_alternative1_hybrid_[kgCO2e/MJ]']
    y3_hybrid_LCFSalternative2 = df_LCFS_scenario_data['LCFS_alternative2_hybrid_[kgCO2e/MJ]']
    y4_hybrid = y1_hybrid.mul(y2_hybrid).mul(y3_hybrid)

    y1_natgas = df[df['Fuel'].isin(["Natural Gas"])].sort_index()["VehicleFuelEconomy_[MJ/vehicle/mile]"]
    y2_natgas = df[df['Fuel'].isin(["Natural Gas"])].sort_index()["VMT_[mi/yr/vehicle]"]
    y2_natgas_BAU =              df_VMT_scenario_data["VMT_natgas_BAU_[mi/yr/vehicle]"]
    y2_natgas_ScopingPlan =      df_VMT_scenario_data["VMT_natgas_ScopingPlan_[mi/yr/vehicle]"] 
    y3_natgas = df[df['Fuel'].isin(["Natural Gas"])].sort_index()["CarbonIntensity_TransportFuel_[kgCO2e/MJ]"]
    y3_natgas_LCFScurrent =      df_LCFS_scenario_data['LCFS_current_natgas_[kgCO2e/MJ]']
    y3_natgas_LCFSproposed =     df_LCFS_scenario_data['LCFS_proposed_natgas_[kgCO2e/MJ]']
    y3_natgas_LCFSalternative1 = df_LCFS_scenario_data['LCFS_alternative1_natgas_[kgCO2e/MJ]']
    y3_natgas_LCFSalternative2 = df_LCFS_scenario_data['LCFS_alternative2_natgas_[kgCO2e/MJ]']
    y4_natgas = y1_natgas.mul(y2_natgas).mul(y3_natgas)

    y1_hydrogen = df[df['Fuel'].isin(["Hydrogen"])].sort_index()["VehicleFuelEconomy_[MJ/vehicle/mile]"]
    y2_hydrogen = df[df['Fuel'].isin(["Hydrogen"])].sort_index()["VMT_[mi/yr/vehicle]"]
    y2_hydrogen_BAU =              df_VMT_scenario_data["VMT_hydrogen_BAU_[mi/yr/vehicle]"]
    y2_hydrogen_ScopingPlan =      df_VMT_scenario_data["VMT_hydrogen_ScopingPlan_[mi/yr/vehicle]"] 
    y3_hydrogen = df[df['Fuel'].isin(["Hydrogen"])].sort_index()["CarbonIntensity_TransportFuel_[kgCO2e/MJ]"]
    y3_hydrogen_LCFScurrent =      df_LCFS_scenario_data['LCFS_current_hydrogen_[kgCO2e/MJ]']
    y3_hydrogen_LCFSproposed =     df_LCFS_scenario_data['LCFS_proposed_hydrogen_[kgCO2e/MJ]']
    y3_hydrogen_LCFSalternative1 = df_LCFS_scenario_data['LCFS_alternative1_hydrogen_[kgCO2e/MJ]']
    y3_hydrogen_LCFSalternative2 = df_LCFS_scenario_data['LCFS_alternative2_hydrogen_[kgCO2e/MJ]']
    y4_hydrogen = y1_hydrogen.mul(y2_hydrogen).mul(y3_hydrogen)

    y1_elec = df[df['Fuel'].isin(["Electricity"])].sort_index()["VehicleFuelEconomy_[MJ/vehicle/mile]"]
    y2_elec = df[df['Fuel'].isin(["Electricity"])].sort_index()["VMT_[mi/yr/vehicle]"]
    y2_elec_BAU =              df_VMT_scenario_data["VMT_elec_BAU_[mi/yr/vehicle]"]
    y2_elec_ScopingPlan =      df_VMT_scenario_data["VMT_elec_ScopingPlan_[mi/yr/vehicle]"] 
    y3_elec = df[df['Fuel'].isin(["Electricity"])].sort_index()["CarbonIntensity_TransportFuel_[kgCO2e/MJ]"]
    y3_elec_LCFScurrent =      df_LCFS_scenario_data['LCFS_current_elec_[kgCO2e/MJ]']
    y3_elec_LCFSproposed =     df_LCFS_scenario_data['LCFS_proposed_elec_[kgCO2e/MJ]']
    y3_elec_LCFSalternative1 = df_LCFS_scenario_data['LCFS_alternative1_elec_[kgCO2e/MJ]']
    y3_elec_LCFSalternative2 = df_LCFS_scenario_data['LCFS_alternative2_elec_[kgCO2e/MJ]']
    y4_elec = y1_elec.mul(y2_elec).mul(y3_elec)       


    # Create a figure and set up the axes
    fig, ax1 = plt.subplots(figsize=(10, 6))

    
    # Create first y-axis on the left side: Fuel economy
    
    curve1, = ax1.plot(y1_gasoline.index.tolist(), y1_gasoline.tolist(), color='goldenrod', label='Gasoline') 
    #curve1, = ax1.plot(y1_diesel.index.tolist(), y1_diesel.tolist(), color='goldenrod', linestyle='--', label='Diesel')
    #curve1, = ax1.plot(y1_hybrid.index.tolist(), y1_hybrid.tolist(), color='goldenrod', linestyle=':', label='Plug-In Hybrid')
    #curve1, = ax1.plot(y1_elec.index.tolist(), y1_elec.tolist(), color='goldenrod', linestyle='-.', label='Electric')
    #curve1, = ax1.plot(y1_natgas.index.tolist(), y1_natgas.tolist(), color='goldenrod', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas')
    #curve1, = ax1.plot(y1_hydrogen.index.tolist(), y1_hydrogen.tolist(), color='goldenrod', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen')
    
    ax1.set_ylim(0, 7.5) # (0,7)
    #plt.legend()
 
    

    # Create a second y-axis on the left side: VMT
    ax2 = ax1.twinx()
    
    #curve2, = ax2.plot(y2_gasoline.index.tolist(), y2_gasoline.tolist(), color='dodgerblue', label='Gasoline_EMFAC')
    curve2, = ax2.plot(y2_gasoline_BAU.index.tolist(), y2_gasoline_BAU.tolist(), color='lightskyblue', label='Gasoline_ScoPlan_BAU') # For poster: label='Scoping Plan: BAU'
    curve2, = ax2.plot(y2_gasoline_ScopingPlan.index.tolist(), y2_gasoline_ScopingPlan.tolist(), color='mediumblue', label='Gasoline_ScoPlan') # For poster: label='Scoping Plan: Proposed'
    #curve2, = ax2.plot(y2_diesel.index.tolist(), y2_diesel.tolist(), color='dodgerblue', linestyle='--', label='Diesel_EMFAC')
    #curve2, = ax2.plot(y2_diesel_BAU.index.tolist(), y2_diesel_BAU.tolist(), color='lightskyblue', linestyle='--', label='Diesel_ScoPlan_BAU')
    #curve2, = ax2.plot(y2_diesel_ScopingPlan.index.tolist(), y2_diesel_ScopingPlan.tolist(), color='mediumblue', linestyle='--', label='Diesel_ScoPlan')
    #curve2, = ax2.plot(y2_hybrid.index.tolist(), y2_hybrid.tolist(), color='dodgerblue', linestyle=':', label='Plug-In Hybrid_EMFAC')
    #curve2, = ax2.plot(y2_hybrid_BAU.index.tolist(), y2_hybrid_BAU.tolist(), color='lightskyblue', linestyle=':', label='Plug-In Hybrid_ScoPlan_BAU')
    #curve2, = ax2.plot(y2_hybrid_ScopingPlan.index.tolist(), y2_hybrid_ScopingPlan.tolist(), color='mediumblue', linestyle=':', label='Plug-In Hybrid_ScoPlan')
    #curve2, = ax2.plot(y2_elec.index.tolist(), y2_elec.tolist(), color='dodgerblue', linestyle='-.', label='Electricity_EMFAC')
    #curve2, = ax2.plot(y2_elec_BAU.index.tolist(), y2_elec_BAU.tolist(), color='lightskyblue', linestyle='-.', label='Electricity_ScoPlan_BAU')
    #curve2, = ax2.plot(y2_elec_ScopingPlan.index.tolist(), y2_elec_ScopingPlan.tolist(), color='mediumblue', linestyle='-.', label='Electricity_ScoPlan')
    #curve2, = ax2.plot(y2_natgas.index.tolist(), y2_natgas.tolist(), color='dodgerblue', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas_EMFAC')
    #curve2, = ax2.plot(y2_natgas_BAU.index.tolist(), y2_natgas_BAU.tolist(), color='lightskyblue', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas_ScoPlan_BAU')
    #curve2, = ax2.plot(y2_natgas_ScopingPlan.index.tolist(), y2_natgas_ScopingPlan.tolist(), color='mediumblue', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas_ScoPlan')
    #curve2, = ax2.plot(y2_hydrogen.index.tolist(), y2_hydrogen.tolist(), color='dodgerblue', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen_EMFAC')
    #curve2, = ax2.plot(y2_hydrogen_BAU.index.tolist(), y2_hydrogen_BAU.tolist(), color='lightskyblue', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen_ScoPlan_BAU')
    #curve2, = ax2.plot(y2_hydrogen_ScopingPlan.index.tolist(), y2_hydrogen_ScopingPlan.tolist(), color='mediumblue', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen_ScoPlan')
    
    ax2.spines['left'].set_position(('outward', 60))
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position('left')
    ax2.set_ylim(0, 15000) # (10000, 30000), (5000,20000)
    #plt.legend()
    plt.legend(loc="lower left")

    # Create a first y-axis on the right side: Fuel CI
    ax3 = ax1.twinx()
    
    #curve3, = ax3.plot(y3_gasoline.index.tolist(), y3_gasoline.tolist(), color='olivedrab', label='Gasoline_EMFAC')
    curve3, = ax3.plot(y3_gasoline_LCFScurrent.index.tolist(), y3_gasoline_LCFScurrent.tolist(), color='olive',  label='Gasoline_LCFS_current') # For poster: 'color='palegreen', label='LCFS: current'
    curve3, = ax3.plot(y3_gasoline_LCFSproposed.index.tolist(), y3_gasoline_LCFSproposed.tolist(), color='darkolivegreen', label='Gasoline_LCFS_proposed') # For poster: label='LCFS: proposed'
    #curve3, = ax3.plot(y3_gasoline_LCFSalternative1.index.tolist(), y3_gasoline_LCFSalternative1.tolist(), color='greenyellow', label='Gasoline_LCFS_alternative_1')
    #curve3, = ax3.plot(y3_gasoline_LCFSalternative2.index.tolist(), y3_gasoline_LCFSalternative2.tolist(), color='palegreen', label='Gasoline_LCFS_alternative_2')
    #curve3, = ax3.plot(y3_diesel.index.tolist(), y3_diesel.tolist(), color='olivedrab', linestyle='--', label='Diesel_EMFAC')
    #curve3, = ax3.plot(y3_diesel_LCFScurrent.index.tolist(), y3_diesel_LCFScurrent.tolist(), color='olive', linestyle='--', label='Diesel_LCFS_current')
    #curve3, = ax3.plot(y3_diesel_LCFSproposed.index.tolist(), y3_diesel_LCFSproposed.tolist(), color='darkolivegreen', linestyle='--', label='Diesel_LCFS_proposed')
    #curve3, = ax3.plot(y3_diesel_LCFSalternative1.index.tolist(), y3_diesel_LCFSalternative1.tolist(), color='greenyellow', linestyle='--', label='Diesel_LCFS_alternative_1')
    #curve3, = ax3.plot(y3_diesel_LCFSalternative2.index.tolist(), y3_diesel_LCFSalternative2.tolist(), color='palegreen', linestyle='--', label='Diesel_LCFS_alternative_2')
    #curve3, = ax3.plot(y3_hybrid.index.tolist(), y3_hybrid.tolist(), color='olivedrab', linestyle=':', label='Plug-In Hybrid')
    #curve3, = ax3.plot(y3_hybrid_LCFScurrent.index.tolist(), y3_hybrid_LCFScurrent.tolist(), color='olive', linestyle=':', label='Plug-In Hybrid_LCFS_current')
    #curve3, = ax3.plot(y3_hybrid_LCFSproposed.index.tolist(), y3_hybrid_LCFSproposed.tolist(), color='darkolivegreen', linestyle=':', label='Plug-In Hybrid_LCFS_proposed')
    #curve3, = ax3.plot(y3_hybrid_LCFSalternative1.index.tolist(), y3_hybrid_LCFSalternative1.tolist(), color='greenyellow', linestyle=':', label='Plug-In Hybrid_LCFS_alternative_1')
    #curve3, = ax3.plot(y3_hybrid_LCFSalternative2.index.tolist(), y3_hybrid_LCFSalternative2.tolist(), color='palegreen', linestyle=':', label='Plug-In Hybrid_LCFS_alternative_2')
    #curve3, = ax3.plot(y3_elec.index.tolist(), y3_elec.tolist(), color='olivedrab', linestyle='-.', label='Electricity')
    #curve3, = ax3.plot(y3_elec_LCFScurrent.index.tolist(), y3_elec_LCFScurrent.tolist(), color='olive', linestyle='-.', label='Electricity_LCFS_current')
    #curve3, = ax3.plot(y3_elec_LCFSproposed.index.tolist(), y3_elec_LCFSproposed.tolist(), color='darkolivegreen', linestyle='-.', label='Electricity_LCFS_proposed')
    #curve3, = ax3.plot(y3_elec_LCFSalternative1.index.tolist(), y3_elec_LCFSalternative1.tolist(), color='greenyellow', linestyle='-.', label='Electricity_LCFS_alternative_1')
    #curve3, = ax3.plot(y3_elec_LCFSalternative2.index.tolist(), y3_elec_LCFSalternative2.tolist(), color='palegreen', linestyle='-.', label='Electricity_LCFS_alternative_2')
    #curve3, = ax3.plot(y3_natgas.index.tolist(), y3_natgas.tolist(), color='olivedrab', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas')
    #curve3, = ax3.plot(y3_natgas_LCFScurrent.index.tolist(), y3_natgas_LCFScurrent.tolist(), color='olive', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas_LCFS_current')
    #curve3, = ax3.plot(y3_natgas_LCFSproposed.index.tolist(), y3_natgas_LCFSproposed.tolist(), color='darkolivegreen', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas_LCFS_proposed')
    #curve3, = ax3.plot(y3_natgas_LCFSalternative1.index.tolist(), y3_natgas_LCFSalternative1.tolist(), color='greenyellow', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas_LCFS_alternative_1')
    #curve3, = ax3.plot(y3_natgas_LCFSalternative2.index.tolist(), y3_natgas_LCFSalternative2.tolist(), color='palegreen', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas_LCFS_alternative_2')
    #curve3, = ax3.plot(y3_hydrogen.index.tolist(), y3_hydrogen.tolist(), color='olivedrab', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen')
    #curve3, = ax3.plot(y3_hydrogen_LCFScurrent.index.tolist(), y3_hydrogen_LCFScurrent.tolist(), color='olive', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen_LCFS_current')
    #curve3, = ax3.plot(y3_hydrogen_LCFSproposed.index.tolist(), y3_hydrogen_LCFSproposed.tolist(), color='darkolivegreen', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen_LCFS_proposed')
    #curve3, = ax3.plot(y3_hydrogen_LCFSalternative1.index.tolist(), y3_hydrogen_LCFSalternative1.tolist(), color='greenyellow', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen_LCFS_alternative_1')
    #curve3, = ax3.plot(y3_hydrogen_LCFSalternative2.index.tolist(), y3_hydrogen_LCFSalternative2.tolist(), color='palegreen', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen_LCFS_alternative_2')
    
    ax3.set_ylim(-0.025, 0.25) # (-0.01, 0.1), (-0.01, 0.15)
    #plt.legend()
    plt.legend(loc="lower center")


    # Create second y-axis on the right side: Total emission
    ax4 = ax1.twinx()
    
    #curve4, = ax4.plot(y4_gasoline.index.tolist(), y4_gasoline.tolist(), color='purple', label='Gasoline')
    #curve4, = ax4.plot(y4_diesel.index.tolist(), y4_diesel.tolist(), color='purple', linestyle='--', label='Diesel')
    #curve4, = ax4.plot(y4_hybrid.index.tolist(), y4_hybrid.tolist(), color='purple', linestyle=':', label='Plug-In Hybrid')
    #curve4, = ax4.plot(y4_elec.index.tolist(), y4_elec.tolist(), color='purple', linestyle='-.', label='Electric')
    #curve4, = ax4.plot(y4_natgas.index.tolist(), y4_natgas.tolist(), color='purple', linestyle=(0, (3, 1, 1, 1)), label='Natural Gas')
    #curve4, = ax4.plot(y4_hydrogen.index.tolist(), y4_hydrogen.tolist(), color='purple', linestyle=(0, (4, 1, 2, 1)), label='Hydrogen')
    
    ax4.spines['right'].set_position(('outward', 60))
    ax4.set_ylim(-2000, 7000) # (-1000, 8000)


    # Set labels and colors for the y-axes
    ax1.set_ylabel('Fuel economy  [MJ / mi]', color='goldenrod')
    ax2.set_ylabel('VMT  [mi / veh / yr]', color='dodgerblue')
    ax3.set_ylabel('Fuel carbon intensity  [kgCO2e / MJ]', color='olivedrab')
    ax4.set_ylabel('CO2 emissions  [kgCO2e / veh / yr]', color='purple')


    # Set labels and colors for the x-axis
    ax1.set_xlabel('Year')
    ax1.set_xlim(1998, 2052)

    # Adjust layout
    plt.tight_layout()
    plt.legend()

    # Show plot
    plt.show()
    


if __name__ == "__main__":
    df, df_LCFS_scenario_data, df_VMT_scenario_data = generate_transport_supplementary(
        #emfac_filepath = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\\NonDownloadableData\\transport\EMFAC2021-EI-202xClass-Statewide-All_CalYrs-Annual-20240408182311.csv',
        #output_filepath = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\Transport_Supplementary.csv'
        emfac_filepath=snakemake.input.emfac_dataset,
        output_filepath=snakemake.output[0]
    )
    #plot_curves(df, df_LCFS_scenario_data, df_VMT_scenario_data)