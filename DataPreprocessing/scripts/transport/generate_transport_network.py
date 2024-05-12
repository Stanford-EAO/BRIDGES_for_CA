
import pandas as pd

#vehicle_population_filepath = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\Vehicle_Population_Last_updated_04-28-2023_ada.xlsx'
#mapping_zip_climatezone_filepath = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\BuildingClimateZonesByZIPCode_ada.xlsx'
#output_filepath = 'C:\\Users\mheyer\Code_Stanford\BRIDGES_for_CA\Data\RawData\Transport\TransportNetwork.csv'

def generate_transport_network(path_vehicle_population: str, 
                               path_mapping_zip_climatezone: str, 
                               output_filepath: str) -> None:

    """ 
    Creates and populates the transport network input file for BRIDGES.
    
    Parameters:
        path_vehicle_population (str):
            Lists in "ZIP"-sheet number of registered vehicles by year, fuel type, and zip code.
        path_mapping_zip_climatezone (str):
            Lists for each California zip code the CEC climate zone.
        output_filepath (str):
            Location for generated transport network csv file. 
    """

    df_vehicle_population = pd.read_excel(path_vehicle_population, sheet_name='ZIP')
    df = df_vehicle_population
    df_mapping_zip_climatezone = pd.read_excel(path_mapping_zip_climatezone)

    # Add a column to df where climate zone is matched to the ZIP code
    df_mapping_zip_climatezone.rename(columns={'Zip Code': 'ZIP'}, inplace=True) # rename column in mapping data frame to enable matching
    df = pd.merge(df, df_mapping_zip_climatezone, on="ZIP", how='left')

    # Drop rows with ZIP codes not belonging to any climate zone or unknown ZIP codes
    df.dropna(subset=['Building CZ'], inplace=True)

    # Drop all data points for years other than 2022 and fuel types other than "Battery Electric (BEV)" and "Gasoline"
    df = df[(df['Data Year'] == 2022)]
    # & ((df['Fuel Type'] == 'Battery Electric (BEV)') | (df['Fuel Type'] == 'Gasoline'))]

    # Sum vehicle counts for electric, gasoline, diesel, natural gas, hydrogen, hybrid for each climate zone
    df_EVs = df[(df['Fuel Type'] == 'Battery Electric (BEV)')]
    df_EVs_summedByClimateZone = df_EVs.groupby('Building CZ')['Number of Vehicles'].sum().reset_index()
    df_EVs_summedByClimateZone.set_index('Building CZ', inplace=True)
    df_Gas = df[(df['Fuel Type'] == 'Gasoline')]
    df_Gas_summedByClimateZone = df_Gas.groupby('Building CZ')['Number of Vehicles'].sum().reset_index()
    df_Gas_summedByClimateZone.set_index('Building CZ', inplace=True)

    df_Diesel = df[(df['Fuel Type'] == 'Diesel')]
    df_Diesel_summedByClimateZone = df_Diesel.groupby('Building CZ')['Number of Vehicles'].sum().reset_index()
    df_Diesel_summedByClimateZone.set_index('Building CZ', inplace=True)

    df_NG = df[(df['Fuel Type'] == 'Natural Gas')]
    df_NG_summedByClimateZone = df_NG.groupby('Building CZ')['Number of Vehicles'].sum().reset_index()
    df_NG_summedByClimateZone.set_index('Building CZ', inplace=True)

    df_H2 = df[(df['Fuel Type'] == 'Fuel Cell (FCEV)')]
    df_H2_summedByClimateZone = df_H2.groupby('Building CZ')['Number of Vehicles'].sum().reset_index()
    df_H2_summedByClimateZone.set_index('Building CZ', inplace=True)

    df_Hybrid = df[(df['Fuel Type'] == 'Gasoline Hybrid') | (df['Fuel Type'] == 'Plug-in Hybrid (PHEV)')]
    df_Hybrid_summedByClimateZone = df_Hybrid.groupby('Building CZ')['Number of Vehicles'].sum().reset_index()
    df_Hybrid_summedByClimateZone.set_index('Building CZ', inplace=True)

    # Create empty data frame for TransportNetwork.csv
    columns = 'Node_ELEC','Node_GAS','DistID_ELEC','DistID_GAS','Service Name','Converter Name','Converter Count [no.]','Lifetime [years]','is_hybrid [bin.]','Upgrade Cost (Low) [$]','Upgrade Cost (High) [$]', 'VMT [mi/yr/vehicle]', 'Carbon Intensity [kgCO2/mile]'
    df_transport_network = pd.DataFrame(columns=columns)

    # Populate the transport network at each node (i.e., climate zone) with the characteristics of each transportation type/converter technology including vehicle count
    considered_transport_types = {
        'Veh_Road_LDA_Gasoline': {'Service Name': 'Transp_PassMiles', 'Lifetime [years]': 14, 'is_hybrid [bin.]': 0, 'Upgrade Cost (Low) [$]': 0, 'Upgrade Cost (High) [$]': 0, 'VMT [mi/yr/vehicle]': 18836, 'Carbon Intensity [kgCO2/mile]': 0.257},
        'Veh_Road_LDA_Diesel': {'Service Name': 'Transp_PassMiles', 'Lifetime [years]': 14, 'is_hybrid [bin.]': 0, 'Upgrade Cost (Low) [$]': 0, 'Upgrade Cost (High) [$]': 0, 'VMT [mi/yr/vehicle]': 18836, 'Carbon Intensity [kgCO2/mile]': 0.206},
        'Veh_Road_LDA_HybridElec': {'Service Name': 'Transp_PassMiles', 'Lifetime [years]': 14, 'is_hybrid [bin.]': 0, 'Upgrade Cost (Low) [$]': 200, 'Upgrade Cost (High) [$]': 200, 'VMT [mi/yr/vehicle]': 18836, 'Carbon Intensity [kgCO2/mile]': 0.134},
        'Veh_Road_LDA_NaturalGas': {'Service Name': 'Transp_PassMiles', 'Lifetime [years]': 14, 'is_hybrid [bin.]': 0, 'Upgrade Cost (Low) [$]': 0, 'Upgrade Cost (High) [$]': 0, 'VMT [mi/yr/vehicle]': 18836, 'Carbon Intensity [kgCO2/mile]': 0.390},
        'Veh_Road_LDA_HydrogenFC': {'Service Name': 'Transp_PassMiles', 'Lifetime [years]': 14, 'is_hybrid [bin.]': 0, 'Upgrade Cost (Low) [$]': 0, 'Upgrade Cost (High) [$]': 0, 'VMT [mi/yr/vehicle]': 18836, 'Carbon Intensity [kgCO2/mile]': 0.255}, #0
        'Veh_Road_LDA_Elec': {'Service Name': 'Transp_PassMiles', 'Lifetime [years]': 14, 'is_hybrid [bin.]': 0, 'Upgrade Cost (Low) [$]': 500, 'Upgrade Cost (High) [$]': 500, 'VMT [mi/yr/vehicle]': 18836, 'Carbon Intensity [kgCO2/mile]': 0}
        }
    for converter in considered_transport_types:
        for climate_zone in range(1,17):
            considered_transport_types[converter]['Node_ELEC'] = climate_zone
            considered_transport_types[converter]['Node_GAS'] = climate_zone
            considered_transport_types[converter]['DistID_ELEC'] = climate_zone
            considered_transport_types[converter]['DistID_GAS'] = climate_zone
            considered_transport_types[converter]['Converter Name'] = converter
            if converter == 'Veh_Road_LDA_Gasoline':
                considered_transport_types[converter]['Converter Count [no.]'] = df_Gas_summedByClimateZone.loc[climate_zone, 'Number of Vehicles']
            elif converter == 'Veh_Road_LDA_Elec':
                considered_transport_types[converter]['Converter Count [no.]'] = df_EVs_summedByClimateZone.loc[climate_zone, 'Number of Vehicles']
            elif converter == 'Veh_Road_LDA_Diesel':
                considered_transport_types[converter]['Converter Count [no.]'] = df_Diesel_summedByClimateZone.loc[climate_zone, 'Number of Vehicles']
            elif converter == 'Veh_Road_LDA_HybridElec':
                considered_transport_types[converter]['Converter Count [no.]'] = df_Hybrid_summedByClimateZone.loc[climate_zone, 'Number of Vehicles']
            elif converter == 'Veh_Road_LDA_NaturalGas':
                considered_transport_types[converter]['Converter Count [no.]'] = df_NG_summedByClimateZone.loc[climate_zone, 'Number of Vehicles']
            elif converter == 'Veh_Road_LDA_HydrogenFC':
                considered_transport_types[converter]['Converter Count [no.]'] = df_H2_summedByClimateZone.loc[climate_zone, 'Number of Vehicles']

            df_transport_network.loc[len(df_transport_network)] = considered_transport_types[converter]

    # Print to csv
    df_transport_network.to_csv(output_filepath, index=False)


if __name__ == "__main__":
    generate_transport_network(
        path_vehicle_population=snakemake.input.vehicle_population, # alternatively: snakemake.input[0]
        path_mapping_zip_climatezone=snakemake.input.zip_climate_zone_matching, # alternatively: snakemake.input[1]
        output_filepath=snakemake.output[0] # Note that the [0] is required
    )