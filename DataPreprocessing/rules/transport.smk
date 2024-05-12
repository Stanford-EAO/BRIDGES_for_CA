""" Rules to create input data for modeling transport in BRIDGES """

rule download_vehicle_population:
    message: "Download vehicle population data from CEC."
    params: url = config["data_preprocessing"]["transport"]["data-sets"]["vehicle-population"]
    output: "".join([BRIDGES_path, "Data/RawData/transport/vehicle_population.xlsx"])
    shell: "curl -sLo {output} {params.url}"

rule download_zip_climate_zone_matching:
    message: "Download matching between climate zones and zip codes from CEC."
    params: url = config["data_preprocessing"]["transport"]["data-sets"]["zip-cz-matching"]
    output: "".join([BRIDGES_path, "Data/RawData/transport/zip-cz-mapping.xlsx"])
    shell: "curl -sLo {output} {params.url}"

rule generate_transport_network:
    message: "Generate transport network input file for BRIDGES."
    input:
        vehicle_population = rules.download_vehicle_population.output[0],
        zip_climate_zone_matching = rules.download_zip_climate_zone_matching.output[0]
    output: "".join([BRIDGES_path, "Data/TransportNetwork.csv"])
    conda: "".join([BRIDGES_path, "DataPreprocessing/environments/transport.yaml"])
    script: "".join([BRIDGES_path, "DataPreprocessing/scripts/transport/generate_transport_network.py"])

rule generate_transport_profile:
    message: "Generate transport profile input file for BRIDGES."
    input: 
        weekday_profile = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/HighHome_100p_NoTimers_weekday_CA_20240105.csv"]),
        weekend_profile = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/HighHome_100p_NoTimers_weekend_CA_20240105.csv"]),
        transport_network = rules.generate_transport_network.output[0]
    output: "".join([BRIDGES_path, "Data/TransportProfiles_ELECNetworkCold.csv"])
    conda: "".join([BRIDGES_path, "DataPreprocessing/environments/transport.yaml"])
    script: "".join([BRIDGES_path, "DataPreprocessing/scripts/transport/generate_transport_profile.py"])