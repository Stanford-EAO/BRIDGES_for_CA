foldername = "Data"

################################################################################
########### TOGGLES ###########
# Turns on/off different model relaxations and functionality

# Indicate whether you want to fix gas flow in the direction specified in the GasTransmission and ElecTransmission csv files (=0) or introduce binary variables to allow bi-directional flow directions (=1)
GASFLOW_DIRECTIONS = 0

# Indiciate whether gas distribution retirement should be contemplated (with the associated savings)
gasdistretirement_allowed = 0

# Indicate whether transmission expansion/retirement should be explicitly considered
# 0 = No, 1 = Binary Yes, 2 = Continuous Yes Expansion Only (supports for elec-side expansion only)
TRANSMISSION_EXPANSION_GAS = 0
TRANSMISSION_EXPANSION_ELEC = 0

# Indicate whether to include constraints that link representative time periods
# for tracking storage state of charge (if = 0, periodicity constraints are imposed for each rep. period)
LINKED_PERIODS_STORAGE = 1
# for generator operations such as min up/down times and ramp rates (if = 0, constraints only apply within each rep. period)
LINKED_PERIODS_GENOPS = 0

# Indicate whether steady-state physics should be simulated for the electric and gas systems
# If = 0, the flows of power and gas will be governed by simple transport models
STEADYSTATE_ELEC = 1
STEADYSTATE_GAS = 0

appliance_decisions = 1
hybrids_allowed = 0
bounding_steady_states = 0              # default 0
toggle_variableNatGasPrice = true

################################################################################
#### CLUSTERING PARAMETERS ####

T_inv = 1               # Number of investment time periods modeled
N_Periods = 6          # Number of representative operational time slices modeled for each investment period
HOURS_PER_PERIOD = 24   # Number of hourly time steps in each rep. op. time slice

T_ops = N_Periods                                           # Number of operational periods simulated for each investment year
t_ops = HOURS_PER_PERIOD                                    # Number of hours simulated for each operational period
HOURS_PER_YEAR = 8760   # hours/year
Periods_Per_Year = Int(HOURS_PER_YEAR/HOURS_PER_PERIOD)     # Number of rep. operational periods per calendar year

# Clustering technique to use for generating representative days
# Options include:
# (a) "average",
# (b) "ward",
# (c) "kmeans"
clustering_case = "ward"
# clustering_case = "kmeans"  
# seed_no = 1234
# Random.seed!(seed_no) 

BaseYear = 2021                             # Initial year
Years = [2025,2030,2035,2040,2045]          # Modeled investment years
AllYears = [2021,2025,2030,2035,2040,2045]  # List of all calendar years

################################################################################
#### CONVERSION CONSTANTS ####

SEC_PER_HOUR = 3600     # sec/hour
MJ_PER_MWh = 3600       # MJ/MWh
MWh_PER_MMBTU = 0.293   # MWh/MMBtu
EF_NG = 53.1/MWh_PER_MMBTU/1000    #tCO2/MWh NG     # Per EPA emissions inventory (only CO2, no CH4 leakage)
HHV_H2 = 12.7  # MJ/standard m3
HHV_CH4 = 37.7  # MJ/standard m3
LHV_H2 = 10.24 # MJ/standard m3
LHV_CH4 = 33.9  # MJ/standard m3

## Gas/Power flow parameters
SLACK_BUS = 20
BASEMVA = 100
SLACK_NODE = 20
PRESSURE_MIN = (3447380/10^6)^2 # 500 psi squared
PRESSURE_MAX = (10342136/10^6)^2 #(5515808/10^6)^2 # 800 psi squared
SLACK_NODE_Pressure = PRESSURE_MAX


################################################################################
######## CASE NAMING & ########
###### MODEL PARAMETERS #######

# For file import, and to permit many sensitivity scenarios, each configuration can be specified by a system and a region which corresponds to the necessary import files
system = "Network"
region = "Cold"
num = ""

### Biomethane
biomethane = "Mid"
max_biomethane_share = 0.1      # annual system-wide limitation on biomethane production (as a share of initial core gas demands)
if biomethane == "High"
    global max_biomethane_share = 0.50
end
if biomethane == "Low"
    global max_biomethane_share = 0.1
end

### Industrials
industrials = "None"            # as % of baseline electricity demands
baselinegasdemand_multiplier = 1
if industrials == "No"
    global baselinegasdemand_multiplier = 0.001
end
if industrials == "Mid"
    global baselinegasdemand_multiplier = 2.5
end
if industrials =="High"
    global baselinegasdemand_multiplier = 5
end

buildingretrofits = "Low"   # "Low" or "High" cost of building retrofit

CleanElecCosts = "Mid"
CleanGasCosts = "Mid"

retirements_case = "NoGasDistRet"
gasdistretirement_forced = [0,0,0,0,0]

### Offsets
offsets_case = "NoOffsets"  # where No offsets = 0, Unlimited Offsets = 1.0
maxOffsets_elec = 1.0*ones(T_inv)                  # % of gross emissions
maxOffsets_gas = 1.0*ones(T_inv)

offsets_Cost = [650, 550, 450, 350, 250]                        # $/tCO2e
# offsets_Cost = [1600, 1400, 1200, 1000, 800]   

GasQuality = "Nodal" # "Annual", "No"

br = 1.0                        # build rate multiplier
transmission_multiplier = 1.0   # electric transmission rating multiplier
forceretire_multiplier = 1.5    # multiplier for upper limit on appliance retirement (as share of natural retirement), min = 1.0


################################################################################
##### EMISSIONS INTENSITY# ####

EITrajectory = "MidEI"

# Specify emissions intensity targets for the electricity sector and gas sector with Slow and Fast sensitivity scenarios possible.
# For reference, a natural gas-fired generator will yield ~500kg/MWh elec.; coal-fired generators will yield ~1000kg/MWh elec.; fossil natural gas delivered for direct-use will release ~181kg/MWh thermal
EI_ElecSector = [400,200,100,50,0] #[500,250,75,50,0]  # kg/MWh electricity generated
EI_GasSector = [200,150,100,50,0] #[200,150,50,15,0]   # kg/MWh gas delivered (to core customers)
if EITrajectory == "SlowEI"
    global EI_ElecSector = [500,500,250,250,0]  # kg/MWh electricity generated
    global EI_GasSector = [200,200,90,90,0.0]   # kg/MWh gas delivered (to core customers)
end
if EITrajectory == "FastEI"
    global EI_ElecSector = [1000,500,100,10,0]  # kg/MWh electricity generated
    global EI_GasSector = [200,100,20,2,0.0]    # kg/MWh gas delivered (to core customers)
end

H2molfrac_max = 0.2
ADDITIONAL_SLACK_NODE = 18
SLACK_GAS = 100000000000 # MW

################################################################################
#### DISCOUNT RATES & WACC ####

WACC = 0.07                         # Weighted average cost of capital (WACC) applied to annualize utility-scale capital investments
WACC_APPLIANCES = 0.15              # Weighted average cost of capital (WACC) applied to annualize customer-scale appliance investments
societal_discounting = 0.05         # discount rate [%] applied to discount future costs to present value
LoadGrowthRate = 0.00               # [%/year] of baseline electricity/gas demand growth 

## Compute the discounting factor for each investment period's annualized costs based on the number of years represented by each period
# Eq. 2.72 in Von Wald thesis.
EndOfCostHorizon = 2050             # Specifies the horizon over which societal costs should be included in objective function
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

# Specify the number of residential and commercial customers on each distribution system
N_ResCust = 600000
N_CommCust = 1000
# Residential system costs are estimated at $350/customer-year
# Commercial system costs are estimated at $1200/cust.-year
Costs_GasDistSys = 350*N_ResCust + 1200*N_CommCust # Dollars per distsyst. per year


################################################################################
#### GAS DISTRIBUTION UTILITY FINANCIAL ASSUMPTIONS ####

# For gas distribution retirement evaluation, we need to estimate the potential avoided costs
# of gas system maintenance and reinvestment. To do this, we use the estimated revenue requirement
# and how it evolves across the planning time horizon using a simplified set of assumptions.
RR_est = Costs_GasDistSys     # Each distribution system has an associated total revenue requirement [$/year]
# Here, we assess the potential annual costs of gas system maintenance, depreciation, and reinvestment
# for two cases: business as usual (BAU) and Accelerated Depreciation (AccDep).
BAUGasSyst_FixedCosts = zeros(T_inv)
AccDepGasSyst_FixedCosts = zeros(T_inv,T_inv)
# Using a simplified set of financial assumptions
equity = 0.10               # return on equity afforded to utility shareholders [%]
debt = 0.04                 # interest rate on debt associated with securitization of the gas system  [%]
shareFOM = 0.1              # share of total annual revenue requirement that is fixed operating costs (as opposed to capital investment)
ReinvestmentRate = 0.025    # % of reinvestment 
AvgDepreciation = 0.03      # average % of depreciation per year
println("Estimated Revenue Requirement = $(RR_est) per year")
RB_est = ((1-shareFOM)*RR_est/(equity + AvgDepreciation))
println("Estimated Ratebase = $(RB_est)")
for i = 2:T_inv
    depTimeHorizon = Years[i] - Years[1]
    syd = depTimeHorizon*(depTimeHorizon+1)/2
    nb = RB_est
    depTimeHorizon_remaining = Years[i] - Years[1]
    for j = 1:T_inv-1
        cost = 0
        if j < i
            for y = 1:(Years[j+1]-Years[j])
                cost = cost + depTimeHorizon_remaining/syd*RB_est + nb*debt
                nb = nb - depTimeHorizon_remaining/syd*RB_est
                depTimeHorizon_remaining = depTimeHorizon_remaining - 1
            end
            AccDepGasSyst_FixedCosts[i,j] = cost/(Years[j+1]-Years[j])
        end
        if j >= i
            AccDepGasSyst_FixedCosts[i,j] = 0
        end
    end
    BAUGasSyst_FixedCosts[i] = sum(RR_est*shareFOM + (equity+AvgDepreciation)*((1-AvgDepreciation+ReinvestmentRate)^y)*((1-shareFOM)*RR_est/(equity + AvgDepreciation)) for y = (AllYears[i]-BaseYear+1):(AllYears[i+1]-BaseYear))/(AllYears[i+1]-AllYears[i])
    println("BAU Gas Costs = $(BAUGasSyst_FixedCosts[i])")
    println("ShutDown Gas Costs = $(AccDepGasSyst_FixedCosts[i,:])")
end

BAUGasSyst_FixedCosts[1] = RR_est
# If you retire the gas system in investment period 1, then you must pay off the entire rate base
# in this year
AccDepGasSyst_FixedCosts[1,1] = RB_est

# Cost of electric distribution infrastructure is based on (Fares, et. al, )
# as a function of the peak electrical demand
Cost_DistributionInfrastructure = 73        # $/kW peak

# Cost values for transmission expansion/retirement modeling
ElecTransmissionCapitalCosts = 0.93 # $/MW-m from "Cost of long-distance energy transmission by different carriers," DeSantis
ElecTransmissionOperatingCosts = 0 # $/MW
GasTransmissionCapitalCosts = 0 # $/km
GasTransmissionOperatingCosts = 0 # $/km


################################################################################
#### STORAGE OPTIONS ####

### Options: FormEnergy and PumpedHydroStorage
FormEnergy_allowed = 0
PHS_allowed = 0


################################################################################
#### PRINT OUT CASE SCENARIOS ####

println("PARAMETERS")
println("T_inv: ", T_inv)
println("N_Periods: ", N_Periods)

println("Clustering case: ", clustering_case)
if clustering_case == "kmeans"
    println("Random seed: ", seed_no)
end

println("Max offset electric: ", maxOffsets_elec)
println("Max offset gas: ", maxOffsets_gas)
println("Offset cost: ", offsets_Cost)

println("Gas quality: ", GasQuality)
println("Max H2 injection frac: ", H2molfrac_max)
println("Max appliance retirement multiplier: ", forceretire_multiplier)
println("Max biomethane share: ", max_biomethane_share)
println("EI trajectory electric: ", EI_ElecSector)
println("EI trajectory gas: ", EI_GasSector)

println("Electric transmission expansion: ", TRANSMISSION_EXPANSION_ELEC)
println("Electric transmission expansion cost: ", ElecTransmissionCapitalCosts, " \$/MW-m")

println("Multi-day Storage: ", FormEnergy_allowed)
println("Pumped Hydro Storage: ", PHS_allowed)