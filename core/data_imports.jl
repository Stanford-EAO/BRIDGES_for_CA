################################################################################
# Baseline Energy Demands
################################################################################
# Import the baseline electrical demands [MWh/hr] across all system nodes
D_Elec2 = CSV.read("$(foldername)/BaselineElectricDemands$(system)$(region).csv",DataFrame)

NODES_ELEC = length(D_Elec2[1,:])       # Number of electrical nodes is specified based on the number of columns in D_Elec2
# Set up a new array to hold electrical demand info.
D_Elec = zeros(8760,NODES_ELEC)         
for n = 1:NODES_ELEC
    D_Elec[:,n] = D_Elec2[:,n]
end

# Import the baseline gas demands [MWh/hr] across all system nodes
D_Gas2 = CSV.read("$(foldername)/BaselineGasDemands$(system)$(region).csv",DataFrame)

NODES_GAS = length(D_Gas2[1,:])         # Number of gas nodes is specified based on the number of columns in D_Gas2
# Set up a new array to hold gas demand info.
D_Gas = zeros(8760,NODES_GAS)
for n = 1:NODES_GAS
    D_Gas[:,n] = baselinegasdemand_multiplier*D_Gas2[:,n]
end

# Names/Numbers of nodes to facilitate look-ups
REGIONS_ELEC = String.(names(D_Elec2))
REGIONS_GAS = String.(names(D_Gas2))

# Maximum fossil gas supply at the boundary/slack node
MAXSLACK = zeros(1,NODES_GAS)   # Set slack supply to zero everywhere except for the boundary node
MAXSLACK[SLACK_NODE] = SLACK_GAS    # [MW] (Arbitrarily large, but not so large as to trigger numerical issues in optimization)
MAXSLACK[ADDITIONAL_SLACK_NODE] = SLACK_GAS    # additional gas supply at non slack node


################################################################################
# End-use appliance demands
################################################################################
EndUseAppliances = CSV.read("$(foldername)/EndUseAppliances$(system)$(num).csv",DataFrame)
APPLIANCES = length(EndUseAppliances[:, :1])        # Number of appliance classes modeled
ApplianceServices = EndUseAppliances[:,5]           # End-use services satisfied by each appliance class
PrimeMover_APPLIANCES = EndUseAppliances[:,6]       # Technology type for each appliance class
InitialAppliancePopulation = EndUseAppliances[:,7]  # Initial appliance population [no. units]
ApplianceLifetime = EndUseAppliances[:,8]           # Expected appliance lifetime [years]
IS_HYBRID = EndUseAppliances[:,9]                   # Indicator for whether the appliance is hybrid gas-electric
upgrade_cost = EndUseAppliances[:,10]               # Building infrastructure upgrade costs associated with transitioning to this appliance [$]
APP_avg_charge = EndUseAppliances[:,12]             # [lbs refrigerant]
APP_eol_charge = EndUseAppliances[:,13]             # [lbs refrigerant]
APP_annual_leak = EndUseAppliances[:,14]            # [%/year]
APP_eol_loss = EndUseAppliances[:,15]               # [% at EOL]
CRF_APPLIANCES = (WACC_APPLIANCES.*(1+WACC_APPLIANCES).^ApplianceLifetime)./((1+WACC_APPLIANCES).^ApplianceLifetime .- 1)   # Capital recovery factor [yr^-1] for annualizing appliance investments

## Create a matrix that maps each appliance to the energy service that it satisfies
################################################################################
SERVICES = length(unique(EndUseAppliances[:,5]))
ServiceList = unique(EndUseAppliances[:,5]) # List of all energy services (i.e., residential space heating, residential water heating, commercial space heating, commercial water heating, etc.)
AppliancesToServices = zeros(APPLIANCES,SERVICES)
# For each appliance, a, put a 1 in the column corresponding to its energy service s
for a = 1:APPLIANCES
    AppliancesToServices[a,findfirst(occursin.([(ApplianceServices[a])],ServiceList))] = 1
end

# Pre-compute the cumulative failure fraction for each appliance in each investment period
# See Eq. 2.9 in Von Wald thesis
################################################################################
cumulativefailurefrac = zeros(APPLIANCES,T_inv,T_inv)
failureProb = zeros(APPLIANCES,150)
failureArchive = CSV.read("$(foldername)/failureProb.csv",DataFrame)
# First, calculate failure probabilities for each appliance class in each year of its lifetime from 1 to 50.
for a = 1:APPLIANCES
    for i = 1:50
#       Using Poisson probability distribution to assess failure fractions
#       failureProb[a,i] = exp(-ApplianceLifetime[a])*(ApplianceLifetime[a]^(i))/factorial(i) 
#       However, the factorial function in Julia won't go over 20! which limits our ability to model long-lived equipment
#       Instead, we use an exogenous file generated using python's factorial function.
        failureProb[a,i] = failureArchive[Int(ApplianceLifetime[a]),i]
    end
    # Ensures that the sum across each row equals 1 (i.e., no appliance lasts longer than 50 years)
    failureProb[a,50] = 1 - sum(failureProb[a,1:49])
end
# Second, compute the cumulative failure fraction for each appliance type, in each investment year
# i.e., cumulativefailurefrac[a,v,t] corresponds to the cumulative failure fraction of appliances of type a
# that were installed in investment period v, that will fail by investment period t.
for a = 1:APPLIANCES
    for v = 1:T_inv
        for t = 1:T_inv
            cumulativefailurefrac[a,v,t] = round(sum(failureProb[a,1:max(Years[t]-Years[v],1)]),digits = 4)  # rounding to avoid numerical issues in the optimization program due to small coefficients
            if t == v
                cumulativefailurefrac[a,v,t] = 0.0
            end
        end
    end
end

## Appliance level energy demand profiles (hourly)
# In MWh/hr per unit, for each hour in a typical year; then must be clustered down
################################################################################
ApplianceProfilesGAS2 = CSV.read("$(foldername)/ApplianceProfiles_GAS$(system)$(region).csv",DataFrame)
ApplianceProfilesELEC2 = CSV.read("$(foldername)/ApplianceProfiles_ELEC$(system)$(region).csv",DataFrame)

ApplianceProfilesGAS = zeros(8760,length(ApplianceProfilesGAS2[1,:]))
ApplianceProfilesELEC = zeros(8760,length(ApplianceProfilesELEC2[1,:]))
# Rounded to avoid introducing numerical issues
for i = 1:length(ApplianceProfilesGAS2[1,:])
    ApplianceProfilesGAS[:,i] = round.(ApplianceProfilesGAS2[:,i], digits = 8)
    ApplianceProfilesELEC[:,i] = round.(ApplianceProfilesELEC2[:,i], digits = 8)
end

## Growth rates used for forecasting and back-casting appliance sales
# Set all growth rates to zero
################################################################################
ServicesGrowthRate = zeros(SERVICES,1)      # %\year
HistoricalGrowthRate = zeros(APPLIANCES,1)  # %\year  
ForecastGrowthRate = zeros(APPLIANCES,1)    # %\year    

# Distribution systems are set up to potentially exist at the sub-transmission nodal level
# i.e., multiple distribution systems may exist and operate independently at the
# same transmission node.
################################################################################
DISTSYS_ELEC = (unique(EndUseAppliances[:,3]))
DISTSYS_GAS  = (unique(EndUseAppliances[:,4]))
DIST_ELEC = length(DISTSYS_ELEC)
DIST_GAS = length(DISTSYS_GAS)

# APP_DistSystemLoc_GAS to tie appliances to gas distribution systems
# APP_DistSystemLoc_ELEC to tie appliances to electric distribution systems
APP_DistSystemLoc_ELEC = zeros(DIST_ELEC, APPLIANCES)
APP_DistSystemLoc_GAS = zeros(DIST_GAS, APPLIANCES)

Loc_ELEC = EndUseAppliances[:,3]
Loc_GAS = EndUseAppliances[:,4]
for a = 1:APPLIANCES
    APP_DistSystemLoc_ELEC[findfirst(occursin.([string(Loc_ELEC[a])],string.(DISTSYS_ELEC))),a] = 1
    APP_DistSystemLoc_GAS[findfirst(occursin.([string(Loc_GAS[a])],string.(DISTSYS_GAS))),a] = 1
end

APPLIANCES_NodalLoc_ELEC = zeros(NODES_ELEC, APPLIANCES)
APPLIANCES_NodalLoc_GAS = zeros(NODES_GAS, APPLIANCES)

Loc_ELEC = EndUseAppliances[:,1]
Loc_GAS = EndUseAppliances[:,2]
for a = 1:APPLIANCES
    APPLIANCES_NodalLoc_ELEC[findfirst(occursin.([string(Loc_ELEC[a])],REGIONS_ELEC)),a] = 1
    APPLIANCES_NodalLoc_GAS[findfirst(occursin.([string(Loc_GAS[a])],REGIONS_GAS)),a] = 1
end

### RONDO EDIT
maxBuild_mult = 1 #6e3
maxBuild_mult_tot = maxBuild_mult
### RONDO EDIT


################################################################################
# Transmission interchanges
################################################################################
TransmissionLinks_ELEC = CSV.read("$(foldername)/ElecTransmission$(system).csv",DataFrame)
EDGES_ELEC = length(TransmissionLinks_ELEC[:,1])
MAXFLOW_ELEC = TransmissionLinks_ELEC[:,3].*transmission_multiplier
Line_Rating = TransmissionLinks_ELEC[:,4].*transmission_multiplier
Line_Reactance = TransmissionLinks_ELEC[:,5]
ExistingUnits_ElecTrans = TransmissionLinks_ELEC[:,6]
### RONDO EDIT
MaxNewUnits_ElecTrans = TransmissionLinks_ELEC[:,7] * maxBuild_mult
### RONDO EDIT
Length_ElecTrans = TransmissionLinks_ELEC[:,8]                          # m
EconomicLifetime_ELECTrans = 50                                         # years
CAPEX_ELECTrans = ElecTransmissionCapitalCosts.*Length_ElecTrans        # $/MW-m * m => $/MW
CRF_ELECTrans = (WACC.*(1+WACC).^EconomicLifetime_ELECTrans)./((1+WACC).^EconomicLifetime_ELECTrans .- 1)
AMMORTIZED_ELECTrans = CAPEX_ELECTrans.*Line_Rating.*(CRF_ELECTrans+ElecTransmissionOperatingCosts)

TransmissionLinks_GAS = CSV.read("$(foldername)/GasTransmission$(system).csv",DataFrame)
EDGES_GAS = length(TransmissionLinks_GAS[:,1])
MAXFLOW_GAS = TransmissionLinks_GAS[:,3]./10
Diameter_Pipes = TransmissionLinks_GAS[:,4]  #[m]
Length_Pipes = TransmissionLinks_GAS[:,5]    #[m]
FrictionFactor_Pipes = TransmissionLinks_GAS[:,6]
ExistingUnits_GasTrans = TransmissionLinks_GAS[:,7]
MaxNewUnits_GasTrans = TransmissionLinks_GAS[:,8]
CompressionRatio_MAX_Branch = TransmissionLinks_GAS[:,9]
# CAPEX_GASTrans = GasTransmissionCapitalCosts.*Length_Pipes./1000
# FOM_GASTrans = GasTransmissionOperatingCosts.*Length_Pipes./1000

## Parameters for gas pipeline flow simulation
################################################################################
Temp_GAS = 300 # [K]
Temp_N = 298.15  # [K]
Pressure_N = 101325   # [Pa]
pi = 3.14
SpecGravity = 0.64
Compressibility = 0.96
# Pressure is in Pascals, which puts the actual pressure variables in Pa^2
# Compressibility and specific gravity will vary depending on actual injection
# of alternative fuels, but for the purposes of flow evauluation we assume them
# constant to avoid a fully nonlinear problem
K = zeros(size(Diameter_Pipes))
K1 = zeros(size(Diameter_Pipes))
K2 = zeros(size(Diameter_Pipes))
V = zeros(size(Diameter_Pipes))
C = zeros(size(Diameter_Pipes))

# Here, we present three different approaches to the gas flow equation:
# (1) General flow equation
for e = 1:EDGES_GAS
    K[e] = 1/((13.2986*Temp_N/Pressure_N)^2*Diameter_Pipes[e]^5/(Length_Pipes[e]*SpecGravity*Temp_GAS*Compressibility*FrictionFactor_Pipes[e]))
    V[e] = pi/4*Diameter_Pipes[e]^2*Length_Pipes[e]       # m3
    C[e] = V[e]*Temp_N/Pressure_N/Compressibility/Temp_GAS
end
# (2) Weymouth equation
for e = 1:EDGES_GAS
    K[e] = 1/((137.2364*Temp_N/Pressure_N)^2*Diameter_Pipes[e]^5.33/(Length_Pipes[e]*SpecGravity*Temp_GAS*Compressibility))
end
M_CH4 = 16/1000                         # kg/mol
UnivGasConstant = 8.314                 # J/mol-K
UnivDensity_GAS = 1/0.024465            # moles/m3
GasConstant = 8.314/M_CH4               # J/kg-K
Density_GAS = 1/0.024465*M_CH4          # kg/m3 (converted from moles/m3)
# (3) Per Correa-Posada, Carlos M., and Pedro Sanchez-Martin. "Integrated power and natural gas model for energy adequacy in short-term operation." IEEE Transactions on Power Systems 30.6 (2014): 3347-3355.
for e = 1:EDGES_GAS
    K1[e] = (pi/4)*Diameter_Pipes[e]^2/GasConstant/Temp_GAS/Compressibility/Density_GAS
    K2[e] = (pi/4)^2*Diameter_Pipes[e]^5/FrictionFactor_Pipes[e]/GasConstant/Temp_N/Compressibility/Density_GAS^2
end

# Correcting all constants to bring pressures up to MPa
K1 = K1.*10^6
K2 = K2.*10^12
K = K./10^12
C = C.*10^6

################################################################################
### Import set of energy supply/storage/demand units
################################################################################
Generators = CSV.read("$(foldername)/Generators$(system).csv",DataFrame)

HourlyVRE2 = CSV.read("$(foldername)/HourlyVRE$(system)$(region).csv",DataFrame)
HourlyVRE = zeros(8760,length(HourlyVRE2[1,:]))
for i = 1:length(HourlyVRE2[1,:])
    # Capacity factors must be greater than 0
    HourlyVRE[:,i] = max.(HourlyVRE2[:,i],0)
end

GEN = length(Generators[:, :1])
PrimeMover_GEN = Generators[:,4]
Fuel_GEN = Generators[:,5]
## add geothermal capacity anyway?
# increase geothermal
idx_geothermal = Generators[!, "Prime Mover"] .== fill("Geothermal EGS", size(Generators[!, 8],1), size(Generators[!, 8],2))    
idx_geothermal = [all(row) for row in eachrow(idx_geothermal)]
# previously 4/20
Generators[idx_geothermal, 8] = vec( fill(4 * 1.0, size(Generators[idx_geothermal, 8],1), size(Generators[idx_geothermal, 8],2)) )
Generators[idx_geothermal, 9] = vec( fill(20 * 1.0, size(Generators[idx_geothermal, 9],1), size(Generators[idx_geothermal, 9],2)) )
#

# retirement by 2030 or 2045: force them to build no new nuclear
if techScenario_Nuclear == "2030" || techScenario_Nuclear == "2045"
    idx_nuclear = Generators[!, "Prime Mover"] .== fill("Nuclear", size(Generators[!, 8],1), size(Generators[!, 8],2))    
    idx_nuclear = [all(row) for row in eachrow(idx_nuclear)]
    #
    Generators[idx_nuclear, 8] = vec( fill(0, size(Generators[idx_nuclear, 8],1), size(Generators[idx_nuclear, 8],2)) )
    Generators[idx_nuclear, 9] = vec( fill(0, size(Generators[idx_nuclear, 9],1), size(Generators[idx_nuclear, 9],2)) )
    # NEW FOR Retirement FOR 1 inv period *****
    # Generators[idx_nuclear, 6] = vec( fill(0, size(Generators[idx_nuclear, 6],1), size(Generators[idx_nuclear, 6],2)) )
end
if techScenario_NGCC == "No"
    ng_primemovers = ["Natural Gas CC-CCS","Natural Gas CC","Natural Gas CT"]
    # ng_primemovers = ["Natural Gas CC-CCS"]
    for i = 1:length(ng_primemovers)
        idx_ngccs = Generators[!, "Prime Mover"] .== fill(ng_primemovers[i], size(Generators[!, 8],1), size(Generators[!, 8],2))    
        idx_ngccs = [all(row) for row in eachrow(idx_ngccs)]
        #
        Generators[idx_ngccs, 8] = vec( fill(0, size(Generators[idx_ngccs, 8],1), size(Generators[idx_ngccs, 8],2)) )
        Generators[idx_ngccs, 9] = vec( fill(0, size(Generators[idx_ngccs, 9],1), size(Generators[idx_ngccs, 9],2)) )
    end
end
if techScenario_OffshoreWind == "Limited Offshore"
    idx_offshorewind = Generators[!, "Prime Mover"] .== fill("OffshoreWind", size(Generators[!, 8],1), size(Generators[!, 8],2))    
    idx_offshorewind = [all(row) for row in eachrow(idx_offshorewind)]
    #
    Generators[idx_offshorewind, 8] = vec( fill(2, size(Generators[idx_offshorewind, 8],1), size(Generators[idx_offshorewind, 8],2)) )
    Generators[idx_offshorewind, 9] = vec( fill(5, size(Generators[idx_offshorewind, 9],1), size(Generators[idx_offshorewind, 9],2)) )
elseif techScenario_OffshoreWind == "No Offshore"
    idx_offshorewind = Generators[!, "Prime Mover"] .== fill("OffshoreWind", size(Generators[!, 8],1), size(Generators[!, 8],2))    
    idx_offshorewind = [all(row) for row in eachrow(idx_offshorewind)]
    #
    Generators[idx_offshorewind, 8] = vec( fill(0, size(Generators[idx_offshorewind, 8],1), size(Generators[idx_offshorewind, 8],2)) )
    Generators[idx_offshorewind, 9] = vec( fill(0, size(Generators[idx_offshorewind, 9],1), size(Generators[idx_offshorewind, 9],2)) )
end
#
#
NumUnits_GEN = Generators[:,6]                  # [units]
UnitSize_GEN = Generators[:,7]                  # [MW]
MaxNewUnitsAnnual_GEN = Generators[:,8].*br * maxBuild_mult     # [units/year]
MaxNewUnitsTotal_GEN = Generators[:,9].*br * maxBuild_mult_tot      # [units]
Pmin_GEN = Generators[:,10]                     # [p.u.]
Pmax_GEN = Generators[:,11]                     # [p.u.]
RampDownRate_GEN = Generators[:,12]             # [p.u.]
RampUpRate_GEN = Generators[:,13]               # [p.u.]
MinUpTime_GEN = Generators[:,14]                # [hours]
MinDownTime_GEN = Generators[:,15]              # [hours]
IS_RENEWABLE = Generators[:,16]                 # [bin.]
HeatRate = Generators[:,17]                     # [MMBtu fuel/MWh elec.]
NG_fueled = Generators[:,18]                    # [bin.]
emissions_factors = Generators[:,19]./1000      # [tCO2/MMBtu fuel]
StartUpCosts = Generators[:,20]                 # [$/start]
EconomicLifetime_GEN = Generators[:,21]         # [years]
Lifetime_GEN = Generators[:,22]                 # [years]
StartupFuel = Generators[:,23]                  # [MMBtu/start]
## scenarios for variable retirement year
# nuclear
idx_nuclear = Generators[!, "Prime Mover"] .== fill("Nuclear", size(Generators[!, "Prime Mover"],1), size(Generators[!, "Prime Mover"],2))    
idx_nuclear = [all(row) for row in eachrow(idx_nuclear)]
Generators[idx_nuclear, "Forced Retirement"] = vec( fill(nuclear_RetirementYear, size(Generators[idx_nuclear, "Forced Retirement"],1), size(Generators[idx_nuclear, "Forced Retirement"],2)) )
#
RetirementYear_GEN = min.(Generators[:,24]+Lifetime_GEN,Generators[:,25])
CRF_GEN = (WACC.*(1+WACC).^EconomicLifetime_GEN)./((1+WACC).^EconomicLifetime_GEN .- 1)


### P2G
PowerToGas = CSV.read("$(foldername)/PowerToGas$(system).csv",DataFrame)
#
P2G = length(PowerToGas[:, :1])
#
PrimeMover_P2G = PowerToGas[:,4]
NumUnits_P2G = PowerToGas[:,5]                  # [units]
UnitSize_P2G = PowerToGas[:,6]                  # [MW]
MaxNewUnitsAnnual_P2G = PowerToGas[:,7].*br * maxBuild_mult * 1.0    # [units/year]
MaxNewUnitsTotal_P2G = PowerToGas[:,8].*br * maxBuild_mult_tot * 1.0     # [units]
Pmin_P2G = PowerToGas[:,9]                      # [p.u.]
Pmax_P2G = PowerToGas[:,10]                     # [p.u.]
RampDownRate_P2G = PowerToGas[:,11]             # [p.u.]
RampUpRate_P2G = PowerToGas[:,12]               # [p.u.]
MinUpTime_P2G = PowerToGas[:,13]                # [hours]
MinDownTime_P2G = PowerToGas[:,14]              # [hours]
eta_P2G = PowerToGas[:,15]                      # [MJ gas/MJ elec.]
eta_P2L = PowerToGas[:,16]                      # [MJ LPG/MJ elec.]
EconomicLifetime_P2G = PowerToGas[:,17]         # [years]
Lifetime_P2G = PowerToGas[:,18]                 # [years]
ISBIOMETHANE = PowerToGas[:,19]                 # [bin.]
ISBIOMASS = PowerToGas[:,20]                    # [bin.]
MoleFracs_P2G = Matrix(PowerToGas[:,23:24])             # [%]
CRF_P2G = (WACC.*(1+WACC).^EconomicLifetime_P2G)./((1+WACC).^EconomicLifetime_P2G .- 1)
RetirementYear_P2G = min.(PowerToGas[:,21]+Lifetime_P2G, PowerToGas[:,22])


### RONDO EDIT
### P2H
PowerToHeat = CSV.read("$(foldername)/PowerToHeat$(system).csv",DataFrame)
#
P2H = length(PowerToHeat[:, :1])
#
PrimeMover_P2H = PowerToHeat[:,4]
NumUnits_P2H = PowerToHeat[:,5]                  # [units]
UnitSize_P2H = PowerToHeat[:,6]                  # [MW]
MaxNewUnitsAnnual_P2H = PowerToHeat[:,7].*br     # [units/year]
MaxNewUnitsTotal_P2H = PowerToHeat[:,8].*br      # [units]
Pmin_P2H = PowerToHeat[:,9]                      # [p.u.]
Pmax_P2H = PowerToHeat[:,10]                     # [p.u.]
RampDownRate_P2H = PowerToHeat[:,11]             # [p.u.]
RampUpRate_P2H = PowerToHeat[:,12]               # [p.u.]
MinUpTime_P2H = PowerToHeat[:,13]                # [hours]
MinDownTime_P2H = PowerToHeat[:,14]              # [hours]
eta_P2H = PowerToHeat[:,15]                      # [MJ heat/MJ elec.]
EconomicLifetime_P2H = PowerToHeat[:,16]         # [years]
Lifetime_P2H = PowerToHeat[:,17]                 # [years]
CRF_P2H = (WACC.*(1+WACC).^EconomicLifetime_P2H)./((1+WACC).^EconomicLifetime_P2H .- 1)
RetirementYear_P2H = min.(PowerToHeat[:,18]+Lifetime_P2H, PowerToHeat[:,19])
### RONDO EDIT




### Electrical Storage
ElectricalStorage = CSV.read("$(foldername)/Storage_ELEC$(system).csv",DataFrame)
### choose storage options
# formEnergy
if FormEnergy_allowed == 0
    idx_allowed = ElectricalStorage[!, "Prime Mover"] .!= fill("Multi-day storage", size(ElectricalStorage[!, "Prime Mover"],1), size(ElectricalStorage[!, "Prime Mover"],2))
    idx_allowed = [all(row) for row in eachrow(idx_allowed)]
    ElectricalStorage = ElectricalStorage[idx_allowed,:]
end
# pumped hydro storage
if PHS_allowed == 0
    idx_allowed = ElectricalStorage[!, "Prime Mover"] .!= fill("Pumped hydro storage", size(ElectricalStorage[!, "Prime Mover"],1), size(ElectricalStorage[!, "Prime Mover"],2))
    idx_allowed = [all(row) for row in eachrow(idx_allowed)]
    ElectricalStorage = ElectricalStorage[idx_allowed,:]
end
# hydrogen storage
if H2Storage_allowed == 0
    idx_allowed = ElectricalStorage[!, "Prime Mover"] .!= fill("Long-duration storage", size(ElectricalStorage[!, "Prime Mover"],1), size(ElectricalStorage[!, "Prime Mover"],2))
    idx_allowed = [all(row) for row in eachrow(idx_allowed)]
    ElectricalStorage = ElectricalStorage[idx_allowed,:]
end
###
STORAGE_ELEC = length(ElectricalStorage[:, :1])
PrimeMover_STORAGE_ELEC = ElectricalStorage[:,4]
NumUnits_STORAGE_ELEC = ElectricalStorage[:,5]                  # [units]
UnitSize_STORAGE_ELEC = ElectricalStorage[:,6]                  # [MW]
MaxNewUnitsAnnual_STORAGE_ELEC = ElectricalStorage[:,7].*br * maxBuild_mult     # [units/year]
MaxNewUnitsTotal_STORAGE_ELEC = ElectricalStorage[:,8].*br  * maxBuild_mult_tot     # [units]
duration_ELEC = ElectricalStorage[:,9]                          # [hours]
eta_charging_ELEC = ElectricalStorage[:,10]                     # [%]
eta_discharging_ELEC = ElectricalStorage[:,11]                  # [%]
eta_loss_ELEC = ElectricalStorage[:,12]                         # [%]
EconomicLifetime_STORAGE_ELEC = ElectricalStorage[:,13]         # [years]
Lifetime_STORAGE_ELEC = ElectricalStorage[:,14]                 # [years]
CRF_STORAGE_ELEC = (WACC.*(1+WACC).^EconomicLifetime_STORAGE_ELEC)./((1+WACC).^EconomicLifetime_STORAGE_ELEC .- 1)
RetirementYear_STORAGE_ELEC = min.(ElectricalStorage[:,15]+Lifetime_STORAGE_ELEC,ElectricalStorage[:,16])

### RONDO EDIT
#
HeatStorage = CSV.read("$(foldername)/Storage_HEAT$(system).csv",DataFrame)
#
H2Heating_allowed = 0
# hydrogen for heating
if H2Heating_allowed == 0
    idx_allowed = HeatStorage[!, "Prime Mover"] .!= fill("Hydrogen-to-Heat", size(HeatStorage[!, "Prime Mover"],1), size(HeatStorage[!, "Prime Mover"],2))
    idx_allowed = [all(row) for row in eachrow(idx_allowed)]
    HeatStorage = HeatStorage[idx_allowed,:]
end
#
STORAGE_HEAT = length(HeatStorage[:, :1])
PrimeMover_STORAGE_HEAT = HeatStorage[:,4]
NumUnits_STORAGE_HEAT = HeatStorage[:,5]                  # [units]
UnitSize_STORAGE_HEAT = HeatStorage[:,6]                  # [MWh]   *** energy capacity
MaxNewUnitsAnnual_STORAGE_HEAT = HeatStorage[:,7].*br * maxBuild_mult     # [units/year]
MaxNewUnitsTotal_STORAGE_HEAT = HeatStorage[:,8].*br  * maxBuild_mult_tot     # [units]
#
maxCharge_HEAT    = HeatStorage[:,9]                          # [MW]
maxDischarge_HEAT = HeatStorage[:,10]                          # [MW]
#
eta_charging_HEAT = HeatStorage[:,11]                     # [%]
eta_discharging_HEAT = HeatStorage[:,12]                  # [%]
eta_loss_HEAT = HeatStorage[:,13]                         # [%]
EconomicLifetime_STORAGE_HEAT = HeatStorage[:,14]         # [years]
Lifetime_STORAGE_HEAT = HeatStorage[:,15]                 # [years]
CRF_STORAGE_HEAT = (WACC.*(1+WACC).^EconomicLifetime_STORAGE_HEAT)./((1+WACC).^EconomicLifetime_STORAGE_HEAT .- 1)
RetirementYear_STORAGE_HEAT = min.(HeatStorage[:,16]+Lifetime_STORAGE_HEAT,HeatStorage[:,17])
### RONDO EDIT



#
GasStorage = CSV.read("$(foldername)/Storage_GAS$(system).csv",DataFrame)
#
STORAGE_GAS = length(GasStorage[:, :1])
PrimeMover_STORAGE_GAS = GasStorage[:,4]
NumUnits_STORAGE_GAS = GasStorage[:,5]                          # [units]
UnitSize_STORAGE_GAS = GasStorage[:,6]                          # [MW]
MaxNewUnitsAnnual_STORAGE_GAS = GasStorage[:,7] * maxBuild_mult                 # [units/year]
MaxNewUnitsTotal_STORAGE_GAS = GasStorage[:,8]  * maxBuild_mult_tot                 # [units]
duration_GAS = GasStorage[:,9]                                  # [hours]
eta_charging_GAS = GasStorage[:,10]                             # [%]
eta_discharging_GAS = GasStorage[:,11]                          # [%]
eta_loss_GAS = GasStorage[:,12]                                 # [%]
EconomicLifetime_STORAGE_GAS = GasStorage[:,13]                 # [years]
Lifetime_STORAGE_GAS = GasStorage[:,14]                         # [years]
CRF_STORAGE_GAS = (WACC.*(1+WACC).^EconomicLifetime_STORAGE_GAS)./((1+WACC).^EconomicLifetime_STORAGE_GAS .- 1)
RetirementYear_STORAGE_GAS = min.(GasStorage[:,15]+Lifetime_STORAGE_GAS,GasStorage[:,16])
MoleFracs_STORAGE = Matrix(GasStorage[:,17:18])
initialStorage_GAS = GasStorage[:,19]                           # [MWh]

# Gas storage facilities are assumed to be maintained regardless of decisions made in optimization
CAPEX_STORAGE_GAS = 0*ones(T_inv,STORAGE_GAS)                   # [$]
FOM_STORAGE_GAS = 0*ones(T_inv,STORAGE_GAS)                     # [$]

# define cost of storage as:
costOfGasStorage = 0.5 / MWh_PER_MMBTU                          # [$/MWh]



################################################################################
### Gas quality tracking information
################################################################################
GAS_COMPONENTS = 2                     # Currently set up for CH4, H2
V_m = 40.87                            # moles/standard m3
MolarMass = [16, 2]                    # kg/kmol
LHV = [50, 120]                        # MJ/kg
MoleFrac_MAX = [1.0, H2molfrac_max]              # kmol/kmol gas
HV_MIN = 40                            # MJ/kg
HV_MAX = 120                           # MJ/kg
MoleFracs_SLACK = zeros(NODES_GAS,GAS_COMPONENTS)
MoleFracs_SLACK[:,1] .= 1.0

## For each source of gas, calculate the molar mass [kg/kmol] and LHV [MJ/kg] of gas provided
MolarMass_SLACK = sum(MoleFracs_SLACK.*transpose(MolarMass), dims = 2)         # [kg/kmol gas]
MolarMass_STORAGE = sum(MoleFracs_STORAGE.*transpose(MolarMass), dims = 2)     # [kg/kmol gas]
MolarMass_P2G = sum(MoleFracs_P2G.*transpose(MolarMass), dims = 2)             # [kg/kmol gas]

LHV_SLACK = sum(MoleFracs_SLACK.*transpose(MolarMass.*LHV), dims = 2)./MolarMass_SLACK        # [MJ/kg gas]
LHV_STORAGE = sum(MoleFracs_STORAGE.*transpose(MolarMass.*LHV), dims = 2)./MolarMass_STORAGE    # [MJ/kg gas]
LHV_P2G = sum(MoleFracs_P2G.*transpose(MolarMass.*LHV), dims = 2)./MolarMass_P2G            # [MJ/kg gas]

################################################################################
### CAPEX, FOM, VOM, and fuel costs
################################################################################
CAPEXLookup = CSV.read("$(foldername)/CAPEXLookup.csv",DataFrame)
FOMLookup = CSV.read("$(foldername)/FOMLookup.csv",DataFrame)
VOMLookup = CSV.read("$(foldername)/VOMLookup.csv",DataFrame)
FuelCostLookup = CSV.read("$(foldername)/FuelCostLookUp.csv",DataFrame)

CAPEX_GEN = zeros(T_inv,GEN)
FOM_GEN = zeros(T_inv,GEN)
VOM_GEN = zeros(T_inv,GEN)
FuelCosts = zeros(T_inv,GEN)
CAPEX_P2G = zeros(T_inv,P2G)
FOM_P2G = zeros(T_inv,P2G)
VOM_P2G = zeros(T_inv,P2G)
CAPEX_STORAGE_ELEC = zeros(T_inv,STORAGE_ELEC)
FOM_STORAGE_ELEC = zeros(T_inv,STORAGE_ELEC)
CAPEX_APPLIANCES = zeros(T_inv, APPLIANCES)
FOM_APPLIANCES = zeros(T_inv, APPLIANCES)
### RONDO EDIT
CAPEX_STORAGE_HEAT = zeros(T_inv,STORAGE_HEAT)
FOM_STORAGE_HEAT = zeros(T_inv,STORAGE_HEAT)
#
CAPEX_P2H = zeros(T_inv,P2H)
FOM_P2H = zeros(T_inv,P2H)
VOM_P2H = zeros(T_inv,P2H)
### RONDO EDIT


## Apply cost multipliers for different energy storage cost scenarios
################################################################################
# define function for scaling
function scaleCAPEX(techType, df, multiplier)
    
    # Get the index of the columns corresponding to the range you want (2020 to 2050)
    start_col = findfirst(names(CAPEXLookup) .== "2020")
    end_col = findlast(names(CAPEXLookup) .== "2050")
    #
    # Filter the DataFrame for rows with "Li-ion battery" in the "Technology" column
    filtered_data = CAPEXLookup[CAPEXLookup[!, "Technology"] .== techType, :]

    # Select columns from '2020' to '2030' for the filtered data
    selected_columns = filtered_data[:, start_col:end_col]

    # Multiply values in selected columns by multiplier
    selected_columns = selected_columns .* multiplier

    CAPEXLookup[CAPEXLookup[!, "Technology"] .== techType, start_col:end_col] .= selected_columns
    
    return df
end

## run it on all costs; if multiplier is 1 then it won't change anything
# Li-ion
techType = "Li-ion battery"
CAPEXLookup = scaleCAPEX(techType, CAPEXLookup, cost_LiIon_multiplier)
println("")
println("Li-ion 2020 cost is: ", CAPEXLookup[CAPEXLookup[!, "Technology"] .== techType, "2020"][1], " \$/kW")
# Fe-Air
techType = "Multi-day storage"
CAPEXLookup = scaleCAPEX(techType, CAPEXLookup, cost_FeAir_multiplier)
println("Fe-Air 2020 cost is: ", CAPEXLookup[CAPEXLookup[!, "Technology"] .== techType, "2020"][1], " \$/kW")
# Hydrogen
techType = "Long-duration storage"
CAPEXLookup = scaleCAPEX(techType, CAPEXLookup, cost_HydrogenStorage_multiplier)
println("Hydrogen 2020 cost is: ", CAPEXLookup[CAPEXLookup[!, "Technology"] .== techType, "2020"][1], " \$/kW")
println("")



## Assign the appropriate cost scenario based on CleanElecCosts and CleanGasCosts
################################################################################
CostScenarios = CSV.read("$(foldername)/CostScenarios.csv",DataFrame)

if CleanCosts == "Low"
    global CostScenarios = CSV.read("$(foldername)/CostScenariosLow$(cost_case).csv",DataFrame)
elseif CleanCosts == "High"
    global CostScenarios = CSV.read("$(foldername)/CostScenariosHigh$(cost_case).csv",DataFrame)
end


# Look up each technology, the associated calendar year in the data tables and assign
# it a cost value
################################################################################
for i = 1:T_inv
    for g = 1:GEN
        subset = findall(in([PrimeMover_GEN[g]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_GEN[i,g] = CAPEXLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_GEN[g]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_GEN[i,g] = FOMLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_GEN[g]]),VOMLookup.Technology)
        scen = findall(in([PrimeMover_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),VOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        VOM_GEN[i,g] = VOMLookup[index, Int(Years[i]-2019)]
        subset = findall(in([Fuel_GEN[g]]),FuelCostLookup.Fuel)
        scen = findall(in([Fuel_GEN[g]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FuelCostLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FuelCosts[i,g] = FuelCostLookup[index, Int(Years[i]-2019)]
    end
    for d = 1:P2G
        subset = findall(in([PrimeMover_P2G[d]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_P2G[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_P2G[i,d] = CAPEXLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_P2G[d]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_P2G[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_P2G[i,d] = FOMLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_P2G[d]]),VOMLookup.Technology)
        scen = findall(in([PrimeMover_P2G[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),VOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        VOM_P2G[i,d] = VOMLookup[index, Int(Years[i]-2019)]
    end
    for s = 1:STORAGE_ELEC
        subset = findall(in([PrimeMover_STORAGE_ELEC[s]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_STORAGE_ELEC[s]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_STORAGE_ELEC[i,s] = CAPEXLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_STORAGE_ELEC[s]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_STORAGE_ELEC[s]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_STORAGE_ELEC[i,s] = FOMLookup[index, Int(Years[i]-2019)]
    end
    ### RONDO EDIT
    for s = 1:STORAGE_HEAT
        subset = findall(in([PrimeMover_STORAGE_HEAT[s]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_STORAGE_HEAT[s]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_STORAGE_HEAT[i,s] = CAPEXLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_STORAGE_HEAT[s]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_STORAGE_HEAT[s]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_STORAGE_HEAT[i,s] = FOMLookup[index, Int(Years[i]-2019)]
    end
    for d = 1:P2H
        subset = findall(in([PrimeMover_P2H[d]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_P2H[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_P2H[i,d] = CAPEXLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_P2H[d]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_P2H[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_P2H[i,d] = FOMLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_P2H[d]]),VOMLookup.Technology)
        scen = findall(in([PrimeMover_P2H[d]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),VOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        VOM_P2H[i,d] = VOMLookup[index, Int(Years[i]-2019)]
    end
    ### RONDO EDIT
    for a = 1:APPLIANCES
        subset = findall(in([PrimeMover_APPLIANCES[a]]),CAPEXLookup.Technology)
        scen = findall(in([PrimeMover_APPLIANCES[a]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),CAPEXLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        CAPEX_APPLIANCES[i,a] = CAPEXLookup[index, Int(Years[i]-2019)]
        subset = findall(in([PrimeMover_APPLIANCES[a]]),FOMLookup.Technology)
        scen = findall(in([PrimeMover_APPLIANCES[a]]),CostScenarios.Technology)
        scen = findall(in([CostScenarios.Cost[scen[1]]]),FOMLookup.Cost)
        index = findall(in(subset),scen)
        index = scen[index[1]]
        FOM_APPLIANCES[i,a] = FOMLookup[index, Int(Years[i]-2019)]
    end
end
