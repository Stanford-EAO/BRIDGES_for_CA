import pandas as pd

#path_speech_day_profile_weekday = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\HighHome_100p_NoTimers_weekday_CA_20240105.csv'
#path_speech_day_profile_weekend = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\HighHome_100p_NoTimers_weekend_CA_20240105.csv' 
#path_transport_network = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\TransportNetwork.csv' 
#output_filepath = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\TransportProfiles_ELECNetworkCold.csv'

def generate_transport_profile(path_speech_day_profile_weekday: str,
                               path_speech_day_profile_weekend: str,
                               path_transport_network: str,
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
    num_of_evs_weekday = 24234832
    num_of_evs_weekend = 48469664
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
        column_header = 'x' + str(index+1) # +1 since BRIDGES notation starts with x0
        TransportProfiles[column_header] = year_profile['Sum [MW]'].values * considered_transport_types[row['Converter Name']]['ElectricityDemand'] # multiply demand profile from speech with factor from specifications

    # Print to csv
    TransportProfiles.to_csv(output_filepath, index=False)


if __name__ == "__main__":
    generate_transport_profile(
        path_speech_day_profile_weekday=snakemake.input.weekday_profile,
        path_speech_day_profile_weekend=snakemake.input.weekend_profile,
        path_transport_network=snakemake.input.transport_network,
        output_filepath=snakemake.output[0]
    )