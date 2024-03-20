""" Rules to create input data for modeling appliances in BRIDGES """

rule dummy_generate_appliance_network:
    message: "Dummy rule that copies the appliance profile from NonDownloadableData to Data"
    input: "".join([BRIDGES_path, "Data/NonDownloadableData/appliances/ApplianceProfiles_ELECNetworkCold.csv"])  # Specifies that file(s) that must exist before the rule starts
    params: base_year = config['params']['BaseYear'] # some parameter needed in the script and pulled from parameter yaml file
    output: "".join([BRIDGES_path, "Data/ApplianceProfiles_ELECNetworkCold.csv"])  # Specifies the file(s) that must exist after executing the script
    script: "".join([BRIDGES_path, "DataPreprocessing/scripts/appliances/dummy_generate_appliance_network.jl"])  # Script that is executed by rule. Should use the file specified in "input", should generate the file specified in "output".
