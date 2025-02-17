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
        profile_weekday_notimer_UH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/UniversalHome_100p_NoTimers_weekday_CA_20240428.csv"]),
        profile_weekend_notimer_UH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/UniversalHome_100p_NoTimers_weekend_CA_20240428.csv"]),
        profile_weekday_notimer_HH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/HighHome_100p_NoTimers_weekday_CA_20240428.csv"]),
        profile_weekend_notimer_HH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/HighHome_100p_NoTimers_weekend_CA_20240428.csv"]),
        profile_weekday_notimer_LHHW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/LowHome_HighWork_100p_NoTimers_weekday_CA_20240428.csv"]),
        profile_weekend_notimer_LHHW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/LowHome_HighWork_100p_NoTimers_weekend_CA_20240428.csv"]),
        profile_weekday_notimer_LHLW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/LowHome_LowWork_100p_NoTimers_weekday_CA_20240428.csv"]),
        profile_weekend_notimer_LHLW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/NoTimer/LowHome_LowWork_100p_NoTimers_weekend_CA_20240428.csv"]),
        profile_weekday_withtimer_UH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/UniversalHome_100p_ALLtimers_weekday_CA_20240428.csv"]),
        profile_weekend_withtimer_UH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/UniversalHome_100p_ALLtimers_weekend_CA_20240428.csv"]),
        profile_weekday_withtimer_HH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/HighHome_100p_ALLtimers_weekday_CA_20240428.csv"]),
        profile_weekend_withtimer_HH   = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/HighHome_100p_ALLtimers_weekend_CA_20240428.csv"]),
        profile_weekday_withtimer_LHHW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/LowHome_HighWork_100p_ALLtimers_weekday_CA_20240428.csv"]),
        profile_weekend_withtimer_LHHW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/LowHome_HighWork_100p_ALLtimers_weekend_CA_20240428.csv"]),
        profile_weekday_withtimer_LHLW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/LowHome_LowWork_100p_ALLtimers_weekday_CA_20240428.csv"]),
        profile_weekend_withtimer_LHLW = "".join([BRIDGES_path, "Data/NonDownloadableData/transport/speech_model/WithTimer/LowHome_LowWork_100p_ALLtimers_weekend_CA_20240428.csv"]),        
        transport_network = rules.generate_transport_network.output[0]
    params:
        transport_scenario_chargingtimer = config["data_preprocessing"]["transport"]["scenarios"]["transport_scenario_chargingtimer"],
        transport_scenario_chargeraccess = config["data_preprocessing"]["transport"]["scenarios"]["transport_scenario_chargeraccess"]
    output: "".join([BRIDGES_path, "Data/TransportProfiles_ELECNetworkCold.csv"])
    conda: "".join([BRIDGES_path, "DataPreprocessing/environments/transport.yaml"])
    script: "".join([BRIDGES_path, "DataPreprocessing/scripts/transport/generate_transport_profile.py"])

rule generate_transport_supplementary:
    message: "Generate transport supplementary input file for BRIDGES."
    input:
        emfac_dataset = "".join([BRIDGES_path, "Data\\NonDownloadableData\\transport\EMFAC2021-EI-202xClass-Statewide-All_CalYrs-Annual-20240408182311.csv"])
    output: "".join([BRIDGES_path, "Data/Transport_Supplementary.csv"])
    conda: "".join([BRIDGES_path, "DataPreprocessing/environments/transport.yaml"])
    script: "".join([BRIDGES_path, "DataPreprocessing/scripts/transport/generate_transport_supplementary.py"])
