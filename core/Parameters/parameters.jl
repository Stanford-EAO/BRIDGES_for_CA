# This julia file imports the parameters from the yaml configuration file and 
# turns them into parameters that are accessible in the BRIDGES core scripts.
# In a few cases first calculations are performed on the parameters.

config = YAML.load_file("core/Parameters/"*CONFIG_FILE_NAME*".yaml")

foldername = config["params"]["foldername"]

################################################################################
########### TOGGLES ###########
# Turns on/off different model relaxations and functionality

# Indicate whether you want to fix gas flow in the direction specified in the GasTransmission and ElecTransmission csv files (=0) or introduce binary variables to allow bi-directional flow directions (=1)
GASFLOW_DIRECTIONS = config["params"]["GASFLOW_DIRECTIONS"]

# Indiciate whether gas distribution retirement should be contemplated (with the associated savings)
gasdistretirement_allowed = config["params"]["gasdistretirement_allowed"]
gasdistretirement_forced = config["params"]["gasdistretirement_forced"]
region_retire = config["params"]["region_retire"]

# Indicate whether transmission expansion/retirement should be explicitly considered
# 0 = No, 1 = Binary Yes, 2 = Continuous Yes Expansion Only (supports for elec-side expansion only)
TRANSMISSION_EXPANSION_GAS = config["params"]["TRANSMISSION_EXPANSION_GAS"]
TRANSMISSION_EXPANSION_ELEC = config["params"]["TRANSMISSION_EXPANSION_ELEC"]

# Indicate whether to include constraints that link representative time periods
# for tracking storage state of charge (if = 0, periodicity constraints are imposed for each rep. period)
LINKED_PERIODS_STORAGE = config["params"]["LINKED_PERIODS_STORAGE"]
# for generator operations such as min up/down times and ramp rates (if = 0, constraints only apply within each rep. period)
LINKED_PERIODS_GENOPS = config["params"]["LINKED_PERIODS_GENOPS"]

# Indicate whether steady-state physics should be simulated for the electric and gas systems
# If = 0, the flows of power and gas will be governed by simple transport models
STEADYSTATE_ELEC = config["params"]["STEADYSTATE_ELEC"]
STEADYSTATE_GAS = config["params"]["STEADYSTATE_GAS"]

appliance_decisions = config["params"]["appliance_decisions"] # whether retired appliances may be replaced by appliances that satisfy same demand but with different fuel
transport_decisions = config["params"]["transport_decisions"] # whether retired vehicles may be replaced by vehicles that satisfy same demand but with different fuel
hybrids_allowed = config["params"]["hybrids_allowed"]
bounding_steady_states = config["params"]["bounding_steady_states"] # default 0
toggle_variableNatGasPrice = config["params"]["toggle_variableNatGasPrice"]
force_retire_gasApps = config["params"]["force_retire_gasApps"]
allsector_emissions_constraint = config["params"]["allsector_emissions_constraint"]
consider_refrigerants = config["params"]["consider_refrigerants"]

################################################################################
#### CLUSTERING PARAMETERS ####

T_inv = config["params"]["T_inv"]                           # Number of investment time periods modeled
N_Periods = config["params"]["N_Periods"]                   # Number of representative operational time slices modeled for each investment period
HOURS_PER_PERIOD = config["params"]["HOURS_PER_PERIOD"]     # Number of hourly time steps in each rep. op. time slice

T_ops = N_Periods                                           # Number of operational periods simulated for each investment year
t_ops = HOURS_PER_PERIOD                                    # Number of hours simulated for each operational period
HOURS_PER_YEAR = config["params"]["HOURS_PER_YEAR"]         # hours/year
Periods_Per_Year = Int(HOURS_PER_YEAR/HOURS_PER_PERIOD)     # Number of rep. operational periods per calendar year

# Clustering technique to use for generating representative days
# Options include:
# (a) "average",
# (b) "ward",
# (c) "kmeans"
clustering_case = config["params"]["clustering_case"]
consider_extremedays = config["params"]["consider_extremedays"]
if clustering_case == "kmeans"  
    seed_no = 1234
    Random.seed!(seed_no)
end

BaseYear = config["params"]["BaseYear"]                     # Initial year
Years = config["params"]["Years"]                           # Modeled investment years
AllYears = config["params"]["AllYears"]                     # List of all calendar years

################################################################################
#### CONVERSION CONSTANTS ####

SEC_PER_HOUR = config["params"]["SEC_PER_HOUR"]             # sec/hour
MJ_PER_MWh = config["params"]["MJ_PER_MWh"]                 # MJ/MWh
MWh_PER_MMBTU = config["params"]["MWh_PER_MMBTU"]           # MWh/MMBtu
EF_NG = config["params"]["EF_NG"]                           # tCO2/MWh NG  This is "53.1/MWh_PER_MMBTU/1000". Per EPA emissions inventory (only CO2, no CH4 leakage)
HHV_H2 = config["params"]["HHV_H2"]                         # MJ/standard m3
HHV_CH4 = config["params"]["HHV_CH4"]                       # MJ/standard m3
LHV_H2 = config["params"]["LHV_H2"]                         # MJ/standard m3
LHV_CH4 = config["params"]["LHV_CH4"]                       # MJ/standard m3

## Gas/Power flow parameters
SLACK_BUS = config["params"]["SLACK_BUS"]
BASEMVA = config["params"]["BASEMVA"]
SLACK_NODE = config["params"]["SLACK_NODE"]
PRESSURE_MIN = config["params"]["PRESSURE_MIN"] 
PRESSURE_MAX = config["params"]["PRESSURE_MAX"] 
SLACK_NODE_Pressure = PRESSURE_MAX


################################################################################
######## CASE NAMING & ########
###### MODEL PARAMETERS #######

# For file import, and to permit many sensitivity scenarios, each configuration can be specified by a system and a region which corresponds to the necessary import files
system = config["params"]["system"]
region = config["params"]["region"]
num = config["params"]["num"]

### Biomethane
biomethane = config["params"]["biomethane"]
max_biomethane_share = config["params"]["max_biomethane_share_scenarios"]["Mid"]            # Annual system-wide limitation on biomethane production (as a share of initial core gas demands)
if biomethane == "High"
    global max_biomethane_share = config["params"]["max_biomethane_share_scenarios"]["High"]
end
if biomethane == "Low"
    global max_biomethane_share = config["params"]["max_biomethane_share_scenarios"]["Low"]
end

### Industrials
industrials = config["params"]["industrials"]            
baselinegasdemand_multiplier = config["params"]["baselinegasdemand_multiplier_scenarios"]["None"] # as % of baseline electricity demands
if industrials == "No"
    global baselinegasdemand_multiplier = config["params"]["baselinegasdemand_multiplier_scenarios"]["No"]
end
if industrials == "Mid"
    global baselinegasdemand_multiplier = config["params"]["baselinegasdemand_multiplier_scenarios"]["Mid"]
end
if industrials =="High"
    global baselinegasdemand_multiplier = config["params"]["baselinegasdemand_multiplier_scenarios"]["High"]
end

### Energy Storage Cost
# Li-ion
costScenario_LiIon = config["params"]["costScenario_LiIon"]
#
cost_LiIon_multiplier = config["params"]["cost_LiIon_multiplier_scenarios"]["Mid"]
if costScenario_LiIon == "Low"
    global cost_LiIon_multiplier = config["params"]["cost_LiIon_multiplier_scenarios"]["Low"]
end
if costScenario_LiIon == "High"
    global cost_LiIon_multiplier = config["params"]["cost_LiIon_multiplier_scenarios"]["High"]
end
        
# Fe-Air
costScenario_FeAir = config["params"]["costScenario_FeAir"]
#
cost_FeAir_multiplier = config["params"]["cost_FeAir_multiplier_scenarios"]["Mid"]
if costScenario_FeAir == "Low"
    global cost_FeAir_multiplier = config["params"]["cost_FeAir_multiplier_scenarios"]["Low"]
end
if costScenario_FeAir == "High"
    global cost_FeAir_multiplier = config["params"]["cost_FeAir_multiplier_scenarios"]["High"]
end

# Hydrogen
costScenario_HydrogenStorage = config["params"]["costScenario_HydrogenStorage"]
#
cost_HydrogenStorage_multiplier = config["params"]["cost_HydrogenStorage_multiplier_scenarios"]["Mid"]
if costScenario_HydrogenStorage == "Low"
    global cost_HydrogenStorage_multiplier = config["params"]["cost_HydrogenStorage_multiplier_scenarios"]["Low"]
end
if costScenario_HydrogenStorage == "High"
    global cost_HydrogenStorage_multiplier = config["params"]["cost_HydrogenStorage_multiplier_scenarios"]["High"]
end


### Transport
transport_scenario_zevmandate = config["params"]["transport_scenario_zevmandate"]
transport_scenario_stockshare = config["params"]["transport_scenario_stockshare"]
transport_scenario_vmt = config["params"]["transport_scenario_vmt"]
transport_scenario_fuelcarbonintensity = config["params"]["transport_scenario_fuelcarbonintensity"]
transport_scenario_fueleconomy = config["params"]["transport_scenario_fueleconomy"]

buildingretrofits = config["params"]["buildingretrofits"]                                   # "Low" or "High" cost of building retrofit
transportretrofits = config["params"]["transportretrofits"]                                 # "Low" or "High" cost of building retrofit

CleanCosts = config["params"]["CleanCosts"]
cost_case = config["params"]["cost_case"]


### Offsets
offsets_case = config["params"]["offsets_case"]                                             # where No offsets = 0, Unlimited Offsets = 1.0
maxOffsets_elec = config["params"]["maxOffsets_elec"] * ones(T_inv)                         # % of gross emissions
maxOffsets_gas = config["params"]["maxOffsets_gas"] * ones(T_inv)
maxOffsets = config["params"]["maxOffsets"] * ones(T_inv)

initialEmissions = config["params"]["initialEmissions"]                    # tCO2

offsets_Cost = config["params"]["offsets_Cost"]                                             # $/tCO2e  

GasQuality = config["params"]["GasQuality"]                                                 # "Nodal", "Annual", "No"

br = config["params"]["br"]                                                                 # build rate multiplier
transmission_multiplier = config["params"]["transmission_multiplier"]                       # electric transmission rating multiplier
forceretire_multiplier = config["params"]["forceretire_multiplier"]                         # multiplier for upper limit on appliance retirement (as share of natural retirement), min = 1.0


################################################################################
##### EMISSIONS INTENSITY# ####

EITrajectory = config["params"]["EITrajectory"]

# Specify emissions intensity targets for the electricity sector and gas sector with Slow and Fast sensitivity scenarios possible.
# For reference, a natural gas-fired generator will yield ~500kg/MWh elec.; coal-fired generators will yield ~1000kg/MWh elec.; fossil natural gas delivered for direct-use will release ~181kg/MWh thermal
EI_ElecSector = config["params"]["EI_ElecSector_scenarios"]["MidEI"]                        # kg/MWh electricity generated
EI_GasSector = config["params"]["EI_GasSector_scenarios"]["MidEI"]                          # kg/MWh gas delivered (to core customers)

TotalEmissions_Allowed = config["params"]["TotalEmissions_Allowed"]

EC_TransportSector = config["params"]["EC_TransportSector_scenarios"]["MidEI"]              # emission constraint for transport sector in t_CO2 emitted annually by non-electric passenger vehicles (2022 in CA: 100e6 t)
if EITrajectory == "SlowEI"
    global EI_ElecSector = config["params"]["EI_ElecSector_scenarios"]["SlowEI"]            # kg/MWh electricity generated
    global EI_GasSector = config["params"]["EI_GasSector_scenarios"]["SlowEI"]              # kg/MWh gas delivered (to core customers)
    global EC_TransportSector = config["params"]["EC_TransportSector_scenarios"]["SlowEI"]  # emission constraint for transport sector in t_CO2 emitted annually by non-electric passenger vehicles (2022 in CA: 100e6 t)
end
if EITrajectory == "FastEI"
    global EI_ElecSector = config["params"]["EI_ElecSector_scenarios"]["FastEI"]            # kg/MWh electricity generated
    global EI_GasSector = config["params"]["EI_GasSector_scenarios"]["FastEI"]              # kg/MWh gas delivered (to core customers)
    global EC_TransportSector = config["params"]["EC_TransportSector_scenarios"]["FastEI"]  # emission constraint for transport sector in t_CO2 emitted annually by non-electric passenger vehicles (2022 in CA: 100e6 t)
end

H2molfrac_max = config["params"]["H2molfrac_max"]
ADDITIONAL_SLACK_NODE = config["params"]["ADDITIONAL_SLACK_NODE"]
SLACK_GAS = config["params"]["SLACK_GAS"]                                                   # MW

################################################################################
#### DISCOUNT RATES & WACC ####

WACC = config["params"]["WACC"]                                                             # Weighted average cost of capital (WACC) applied to annualize utility-scale capital investments
WACC_APPLIANCES = config["params"]["WACC_APPLIANCES"]                                       # Weighted average cost of capital (WACC) applied to annualize customer-scale appliance investments
WACC_TRANSPORT = config["params"]["WACC_TRANSPORT"]                                         # Weighted average cost of capital (WACC) applied to annualize customer-scale vehicle investments
societal_discounting = config["params"]["societal_discounting"]                             # discount rate [%] applied to discount future costs to present value
LoadGrowthRate = config["params"]["LoadGrowthRate"]                                         # [%/year] of baseline electricity/gas demand growth 

## Compute the discounting factor for each investment period's annualized costs based on the number of years represented by each period
# Eq. 2.72 in Von Wald thesis.
EndOfCostHorizon = config["params"]["EndOfCostHorizon"]                                     # Specifies the horizon over which societal costs should be included in objective function
discountfactor = zeros(T_inv)
for i = 1:T_inv
    if i < T_inv
        for j = 1:Int(Years[i+1]-Years[i])
            discountfactor[i] = discountfactor[i] + 1/((1+societal_discounting)^(Years[i]-BaseYear+j-1))
        end
    end
    if i == T_inv
        for j = 1:Int(EndOfCostHorizon - Years[i] - 1)
            discountfactor[i] = discountfactor[i] + 1/((1+societal_discounting)^(Years[i]-BaseYear+j-1))
        end
    end
end

################################################################################
#### GAS DISTRIBUTION UTILITY FINANCIAL ASSUMPTIONS ####

# Using a simplified set of financial assumptions
equity = config["params"]["equity"]                                                         # return on equity afforded to utility shareholders [%]
debt = config["params"]["debt"]                                                             # interest rate on debt associated with securitization of the gas system  [%]
shareFOM = config["params"]["shareFOM"]                                                     # share of total annual revenue requirement that is fixed operating costs (as opposed to capital investment)
ReinvestmentRate = config["params"]["ReinvestmentRate"]                                     # % of reinvestment 
AvgDepreciation = config["params"]["AvgDepreciation"]                                       # average % of depreciation per year

# Cost of electric distribution infrastructure is based on (Fares, et. al, ) as a function of the peak electrical demand
Cost_DistributionInfrastructure = config["params"]["Cost_DistributionInfrastructure"]       # $/kW peak

# Cost values for transmission expansion/retirement modeling
ElecTransmissionCapitalCosts = config["params"]["ElecTransmissionCapitalCosts"]             # $/MW-m from "Cost of long-distance energy transmission by different carriers," DeSantis
ElecTransmissionOperatingCosts = config["params"]["ElecTransmissionOperatingCosts"]         # $/MW
GasTransmissionCapitalCosts = config["params"]["GasTransmissionCapitalCosts"]               # $/km
GasTransmissionOperatingCosts = config["params"]["GasTransmissionOperatingCosts"]           # $/km

################################################################################
#### FIRM GEN OPTIONS ####
# nuclear retirement year 
techScenario_Nuclear = config["params"]["techScenario_Nuclear"]
# unrestricted
nuclear_RetirementYear = config["params"]["nuclear_RetirementYear"]                   # normal retirement year, like all other generators
# retirement by 2030
if techScenario_Nuclear == "2030"
    global nuclear_RetirementYear = config["params"]["nuclear_RetirementYear_scenarios"][2030] 
end
# retirement by 2045
if techScenario_Nuclear == "2045"
    global nuclear_RetirementYear = config["params"]["nuclear_RetirementYear_scenarios"][2045] 
end

techScenario_OffshoreWind = config["params"]["techScenario_OffshoreWind"] # "No Offshore" activates restriction
techScenario_NGCC = config["params"]["techScenario_NGCC"] # if "No" restricts to no new build of NG CC,CT,CC-CCS

################################################################################
#### STORAGE OPTIONS ####

### Options: FormEnergy and PumpedHydroStorage and HydrogenStorage
FormEnergy_allowed = config["params"]["FormEnergy_allowed"]
PHS_allowed = config["params"]["PHS_allowed"]
H2Storage_allowed = config["params"]["H2Storage_allowed"]
# starting SOC
SOC_fraction = config["params"]["SOC_fraction"]

### Heat Storage and nonGasHeat_ON
nonGasHeat_ON = config["params"]["nonGasHeat_ON"]
fraction_electrifiableHeat = config["params"]["fraction_electrifiableHeat"]
simpleHeatElectrification_ON = config["params"]["simpleHeatElectrification_ON"]   # simple == without heat storage


GWP100_methane = config["params"]["GWP100_methane"]
EnergyContent_methane = config["params"]["EnergyContent_methane"]
methaneLeak_ON = config["params"]["methaneLeak_ON"]
if methaneLeak_ON == 1
    methane_leakage = config["params"]["methaneLeak_scenarios"]["On"]
else
    methane_leakage = config["params"]["methaneLeak_scenarios"]["Off"]
end

################################################################################
#### PRINT OUT CASE SCENARIOS ####

println("")
#
println("PARAMETERS")
println("T_inv: ", T_inv)
println("N_Periods: ", N_Periods)
println("")

println("Clustering case: ", clustering_case)
if clustering_case == "kmeans"
    println("Random seed: ", seed_no)
end
println("Consider extreme days: ", consider_extremedays)

println("Linked storage: ", LINKED_PERIODS_STORAGE)
println("Linked generation: ", LINKED_PERIODS_GENOPS)
println("Steady state elec: ", STEADYSTATE_ELEC)
println("Steady state gas: ", STEADYSTATE_GAS)

println("Constraint emissions of all sector: ", allsector_emissions_constraint)
if allsector_emissions_constraint == 1
    println("Total emissions allowed [MMTCO2e]: ", TotalEmissions_Allowed)
    println("Max offset as % initial emissions: ", maxOffsets)
    println("Initial emissions [MMTCO2e]: ", initialEmissions/1e6)
else
    println("EI trajectory electric: ", EI_ElecSector)
    println("EI trajectory gas: ", EI_GasSector)
    println("Max offset electric: ", maxOffsets_elec)
    println("Max offset gas: ", maxOffsets_gas)
end
println("Offset cost: ", offsets_Cost)

println("Gas quality: ", GasQuality)
println("Max H2 injection frac: ", H2molfrac_max)
println("Max appliance retirement multiplier: ", forceretire_multiplier)
println("Max biomethane share: ", max_biomethane_share)

println("Allow gas dist. retirement: ", gasdistretirement_allowed)
println("Force gas dist. retirement: ", gasdistretirement_forced)
if sum(gasdistretirement_forced) > 0
    println("Region forced retire: ", region_retire)
end

println("Electric transmission expansion: ", TRANSMISSION_EXPANSION_ELEC)
println("Electric transmission expansion cost: ", ElecTransmissionCapitalCosts, " \$/MW-m")
println("")

println("Nuclear Retirement Year: ", nuclear_RetirementYear)
println("Offshore Wind Build: ", techScenario_OffshoreWind)
println("NG New Build: ", techScenario_NGCC)
println("")

println("Appliance ban by 2045: ", force_retire_gasApps)
println("Consider refrigerants: ", consider_refrigerants)

println("Multi-day Storage: ", FormEnergy_allowed)
println("Pumped Hydro Storage: ", PHS_allowed)
println("Hydrogen Storage: ", H2Storage_allowed)
println("")

println("Li-ion Cost Multiplier: ", cost_LiIon_multiplier)
println("Fe-Air Cost Multiplier: ", cost_FeAir_multiplier)
println("H2 Storage Cost Multiplier: ", cost_HydrogenStorage_multiplier)
println("")

### RONDO EDIT
println("Non-gas Heat Allowed: ", nonGasHeat_ON)
println("Simple Heat Electrification: ", simpleHeatElectrification_ON)
println("Fraction of Simple Heat Electrification: ", fraction_electrifiableHeat)
println("Tracking Methane Leakage: ", if methaneLeak_ON == 1 "Yes with a $(methane_leakage*100) % leakage" else "No" end)
println("")
println("")
#
### RONDO EDIT

println("Starting SOC: ", SOC_fraction)
