""" Rules to create input data for modeling storage in BRIDGES """

rule dummy_generate_storage_network:
    message: "Dummy rule that copies the storage network from NonDownloadableData to Data"
    input: "".join([BRIDGES_path, "Data/NonDownloadableData/storage/Storage_ELECNetwork.csv"])  # Specifies that file(s) that must exist before the rule starts
    params: base_year = config['params']['BaseYear'] # some parameter needed in the script and pulled from parameter yaml file
    output: "".join([BRIDGES_path, "Data/Storage_ELECNetwork.csv"])  # Specifies the file(s) that must exist after executing the script
    conda: "".join([BRIDGES_path, "DataPreprocessing/environments/transport.yaml"])
    script: "".join([BRIDGES_path, "DataPreprocessing/scripts/storage/dummy_generate_storage_network.py"])  # Script that is executed by rule. Should use the file specified in "input", should generate the file specified in "output".
