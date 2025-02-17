import pandas as pd

#path_speech_day_profile_weekday = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\HighHome_100p_NoTimers_weekday_CA_20240105.csv'
#path_speech_day_profile_weekend = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\HighHome_100p_NoTimers_weekend_CA_20240105.csv' 
#path_transport_network = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\TransportNetwork.csv' 
#output_filepath = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\TransportProfiles_ELECNetworkCold.csv'

def generate_transport_profile(profile_weekday_notimer_UH: str,
                               profile_weekend_notimer_UH: str,
                               profile_weekday_notimer_HH: str,
                               profile_weekend_notimer_HH: str,
                               profile_weekday_notimer_LHHW: str,
                               profile_weekend_notimer_LHHW: str,
                               profile_weekday_notimer_LHLW: str,
                               profile_weekend_notimer_LHLW: str,
                               profile_weekday_withtimer_UH: str,
                               profile_weekend_withtimer_UH: str,
                               profile_weekday_withtimer_HH: str,
                               profile_weekend_withtimer_HH: str,
                               profile_weekday_withtimer_LHHW: str,
                               profile_weekend_withtimer_LHHW: str,
                               profile_weekday_withtimer_LHLW: str,
                               profile_weekend_withtimer_LHLW: str,
                               path_transport_network: str,
                               transport_scenario_chargingtimer: str,
                               transport_scenario_chargeraccess: str,
                               output_filepath: str
                               ) -> None:

    """ 
    Creates the transport profile input file for BRIDGES.
    
    Parameters:
        path_speech_day_profile_weekday (str):
            For every minute of the (week)day lists the total EV electricity demand for CA by charging location.
        path_speech_day_profile_weekend (str):
            For every minute of the (weekend)day lists the total EV electricity demand for CA by charging location.
        path_transport_network (str):
            Lists all transport technology (EV, Diesel, ...) at each node the initial population, lifetime, etc.
        output_filepath (str):
            Lists for each transport technology at each node the electricity demand per car for every hour of the year. 
    """

    # Specifications
    # Vehicle numbers are added manually. Note that there is a problem with the vehicle numbers for weekenddays produced by the jupyter script (they are doubled).
    num_of_evs_weekday = 24234832
    num_of_evs_weekend = 24234832  # As numbers seem to be same, might be able to simplify duplicate code in this script.
    if transport_scenario_chargingtimer == "no_timer":
        if transport_scenario_chargeraccess == "UniversalHome":
            path_speech_day_profile_weekday = profile_weekday_notimer_UH
            path_speech_day_profile_weekend = profile_weekend_notimer_UH
        elif transport_scenario_chargeraccess == "HighHome":
            path_speech_day_profile_weekday = profile_weekday_notimer_HH
            path_speech_day_profile_weekend = profile_weekend_notimer_HH
        elif transport_scenario_chargeraccess == "LowHome_HighWork":
            path_speech_day_profile_weekday = profile_weekday_notimer_LHHW
            path_speech_day_profile_weekend = profile_weekend_notimer_LHHW
        elif transport_scenario_chargeraccess == "LowHome_LowWork":
            path_speech_day_profile_weekday = profile_weekday_notimer_LHLW
            path_speech_day_profile_weekend = profile_weekend_notimer_LHLW
    elif transport_scenario_chargingtimer == "with_timer":
        if transport_scenario_chargeraccess == "UniversalHome":
            path_speech_day_profile_weekday = profile_weekday_withtimer_UH
            path_speech_day_profile_weekend = profile_weekend_withtimer_UH
        elif transport_scenario_chargeraccess == "HighHome":
            path_speech_day_profile_weekday = profile_weekday_withtimer_HH
            path_speech_day_profile_weekend = profile_weekend_withtimer_HH
        elif transport_scenario_chargeraccess == "LowHome_HighWork":
            path_speech_day_profile_weekday = profile_weekday_withtimer_LHHW
            path_speech_day_profile_weekend = profile_weekend_withtimer_LHHW
        elif transport_scenario_chargeraccess == "LowHome_LowWork":
            path_speech_day_profile_weekday = profile_weekday_withtimer_LHLW
            path_speech_day_profile_weekend = profile_weekend_withtimer_LHLW

    number_of_nodes = 16
    considered_transport_types = {
        'Veh_Road_LDA_Gasoline': {'ElectricityDemand': 0}, # Electricity demand relative to full battery electric vehicles. Demand profiles are scaled with this value.
        'Veh_Road_LDA_Diesel': {'ElectricityDemand': 0},
        'Veh_Road_LDA_HybridElec': {'ElectricityDemand': 0.46},
        'Veh_Road_LDA_NaturalGas': {'ElectricityDemand': 0},
        'Veh_Road_LDA_HydrogenFC': {'ElectricityDemand': 0},
        'Veh_Road_LDA_Elec': {'ElectricityDemand': 1},
        }

    day_profile_weekday = pd.read_csv(path_speech_day_profile_weekday, index_col=0)
    day_profile_weekend = pd.read_csv(path_speech_day_profile_weekend, index_col=0)
    transport_network = pd.read_csv(path_transport_network)

    #print(transport_network)

    # Convert data from kW (speech) to MW (BRIDGES)
    day_profile_weekday = day_profile_weekday/1000
    day_profile_weekend = day_profile_weekend/1000

    # Convert from minute resolution to hourly resolution by averaging demand over every block of 60 rows
    day_profile_weekday = day_profile_weekday.groupby(day_profile_weekday.index // 60).agg('mean')
    day_profile_weekend = day_profile_weekend.groupby(day_profile_weekend.index // 60).agg('mean')
    day_profile_weekday.reset_index(drop=True, inplace=True)
    day_profile_weekend.reset_index(drop=True, inplace=True)

    # Aggregate the demand over all charging types ("segments": Residential L1, Residential L2, MUD L2, Workplace L2, Public L2, Public DCFC) by summing over columns
    day_profile_weekday['Sum [MW]'] = day_profile_weekday.sum(axis=1)
    day_profile_weekend['Sum [MW]'] = day_profile_weekend.sum(axis=1)

    # Get average demand per vehicle, by dividing by number of cars in modeled area
    day_profile_weekday = day_profile_weekday/num_of_evs_weekday
    day_profile_weekend = day_profile_weekend/num_of_evs_weekend

    # Concatenate the weekday and weekend profiles into a week profile
    week_profile = pd.concat([day_profile_weekday] * 5 + [day_profile_weekend] * 2, ignore_index=True)

    # Concatenate the week profile into an annual profile (concatenate 53 times and trim to 8760 hours)
    year_profile = pd.concat([week_profile] * 53, ignore_index=True).iloc[:8760]

    # Create transport profile data input file for BRIDGES. Use same per-car profile for EVs in all regions. Assumes BRIDGES input format with x1, x2, denoting the regions and transport types.
    TransportProfiles = pd.DataFrame(index=range(8760))
    for index, row in transport_network.iterrows():
        column_header = 'x' + str(index+1) # +1 since BRIDGES notation starts with x1
        TransportProfiles[column_header] = year_profile['Sum [MW]'].values * considered_transport_types[row['Converter Name']]['ElectricityDemand'] # multiply demand profile from speech with factor from specifications

    # Print to csv
    TransportProfiles.to_csv(output_filepath, index=False)


if __name__ == "__main__":
    generate_transport_profile(
        profile_weekday_notimer_UH   =snakemake.input.profile_weekday_notimer_UH,
        profile_weekend_notimer_UH   =snakemake.input.profile_weekend_notimer_UH,
        profile_weekday_notimer_HH   =snakemake.input.profile_weekday_notimer_HH,
        profile_weekend_notimer_HH   =snakemake.input.profile_weekend_notimer_HH,
        profile_weekday_notimer_LHHW =snakemake.input.profile_weekday_notimer_LHHW,
        profile_weekend_notimer_LHHW =snakemake.input.profile_weekend_notimer_LHHW,
        profile_weekday_notimer_LHLW =snakemake.input.profile_weekday_notimer_LHLW,
        profile_weekend_notimer_LHLW =snakemake.input.profile_weekend_notimer_LHLW,
        profile_weekday_withtimer_UH   =snakemake.input.profile_weekday_withtimer_UH,
        profile_weekend_withtimer_UH   =snakemake.input.profile_weekend_withtimer_UH,
        profile_weekday_withtimer_HH   =snakemake.input.profile_weekday_withtimer_HH,
        profile_weekend_withtimer_HH   =snakemake.input.profile_weekend_withtimer_HH,
        profile_weekday_withtimer_LHHW =snakemake.input.profile_weekday_withtimer_LHHW,
        profile_weekend_withtimer_LHHW =snakemake.input.profile_weekend_withtimer_LHHW,
        profile_weekday_withtimer_LHLW =snakemake.input.profile_weekday_withtimer_LHLW,
        profile_weekend_withtimer_LHLW =snakemake.input.profile_weekend_withtimer_LHLW,
        path_transport_network=snakemake.input.transport_network,
        transport_scenario_chargingtimer=snakemake.params.transport_scenario_chargingtimer,
        transport_scenario_chargeraccess=snakemake.params.transport_scenario_chargeraccess,
        output_filepath=snakemake.output[0]
    )