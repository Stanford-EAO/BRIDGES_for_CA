###############################################################################
### Expansion and retirement of energy supply/demand units
# See Eq. 2.5a/2.5b/2.5c in Von Wald thesis
################################################################################
@variable(m, 0 <= unitsbuilt_GEN[I = 1:T_inv, g = 1:GEN] <= MaxNewUnitsAnnual_GEN[g])                                      # [units]
@variable(m, 0 <= unitsretired_GEN[I = 1:T_inv, g = 1:GEN])                                                                 # [units]
@constraint(m, [g = 1:GEN], sum(unitsbuilt_GEN[i,g] for i = 1:T_inv) <= MaxNewUnitsTotal_GEN[g])                            # [units]
# new constraint to check if we can limit building
# @constraint(m, [I = 1:T_inv, g = 1:GEN], unitsbuilt_GEN[I, g] <= MaxNewUnitsAnnual_GEN[g] * min( max(RetirementYear_GEN[g] - Years[I],0), 1) )                                      # [units]

@variable(m, 0 <= unitsbuilt_P2G[I = 1:T_inv, d = 1:P2G] <= MaxNewUnitsAnnual_P2G[d])                                      # [units]
@variable(m, 0 <= unitsretired_P2G[I = 1:T_inv, d = 1:P2G])                                                                 # [units]
@constraint(m, [d = 1:P2G], sum(unitsbuilt_P2G[i,d] for i = 1:T_inv) <= MaxNewUnitsTotal_P2G[d])                            # [units]

@variable(m, 0 <= unitsbuilt_STORAGE_ELEC[I = 1:T_inv, s = 1:STORAGE_ELEC] <= MaxNewUnitsAnnual_STORAGE_ELEC[s])           # [units]
@variable(m, 0 <= unitsretired_STORAGE_ELEC[I = 1:T_inv, s = 1:STORAGE_ELEC])                                               # [units]
@constraint(m, [s = 1:STORAGE_ELEC], sum(unitsbuilt_STORAGE_ELEC[i,s] for i = 1:T_inv) <= MaxNewUnitsTotal_STORAGE_ELEC[s]) # [units]

### RONDO EDIT
@variable(m, 0 <= unitsbuilt_STORAGE_HEAT[I = 1:T_inv, s = 1:STORAGE_HEAT] <= MaxNewUnitsAnnual_STORAGE_HEAT[s])           # [units]
@variable(m, 0 <= unitsretired_STORAGE_HEAT[I = 1:T_inv, s = 1:STORAGE_HEAT])                                               # [units]
@constraint(m, [s = 1:STORAGE_HEAT], sum(unitsbuilt_STORAGE_HEAT[i,s] for i = 1:T_inv) <= MaxNewUnitsTotal_STORAGE_HEAT[s]) # [units]
# #
@variable(m, 0 <= unitsbuilt_P2H[I = 1:T_inv, d = 1:P2H] <= MaxNewUnitsAnnual_P2H[d])                                      # [units]
@variable(m, 0 <= unitsretired_P2H[I = 1:T_inv, d = 1:P2H])                                                                 # [units]
@constraint(m, [d = 1:P2H], sum(unitsbuilt_P2H[i,d] for i = 1:T_inv) <= MaxNewUnitsTotal_P2H[d])                            # [units]
### RONDO EDIT

@variable(m, 0 <= unitsbuilt_STORAGE_GAS[I = 1:T_inv, s = 1:STORAGE_GAS] <= MaxNewUnitsAnnual_STORAGE_GAS[s])              # [units]
@variable(m, 0 <= unitsretired_STORAGE_GAS[I = 1:T_inv, s = 1:STORAGE_GAS])                                                 # [units]
@constraint(m, [s = 1:STORAGE_GAS], sum(unitsbuilt_STORAGE_GAS[i,s] for i = 1:T_inv) <= MaxNewUnitsTotal_STORAGE_GAS[s])    # [units]

@variable(m, unitsbuilt_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)       # [thousands of units]
@variable(m, unitsremaining_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)   # [thousands of units]
@variable(m, unitsretired_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)     # [thousands of units]


###############################################################################
### Retirement functions for generators, p2g, p2h, and storage units
# See Eq. 2.5d/2.5e in Von Wald thesis 
################################################################################                       
@constraint(m, [I = 1, g = 1:GEN], unitsretired_GEN[I,g] <= NumUnits_GEN[g] + unitsbuilt_GEN[I,g])                          
@constraint(m, [I = 1, g = 1:GEN], unitsretired_GEN[I,g] >= NumUnits_GEN[g]*max(min(Years[I] - RetirementYear_GEN[g],1),0)) 
if T_inv > 1
    @constraint(m, [I = 2:T_inv, g = 1:GEN], unitsretired_GEN[I,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i0,g] - unitsretired_GEN[i0,g]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, g = 1:GEN], unitsretired_GEN[I,g] >= NumUnits_GEN[g]*max(min(Years[I] - RetirementYear_GEN[g],1),0) + sum(unitsbuilt_GEN[i0,g]*max(min(Years[I] - (Years[i0] + Lifetime_GEN[g]),1),0) for i0 = 1:I-1) - sum(unitsretired_GEN[i0,g] for i0 = 1:I-1))
end

@constraint(m, [I = 1, d = 1:P2G], unitsretired_P2G[I,d] <= NumUnits_P2G[d] + unitsbuilt_P2G[I,d])
@constraint(m, [I = 1, d = 1:P2G], unitsretired_P2G[I,d] >= NumUnits_P2G[d]*max(min(Years[I] - RetirementYear_P2G[d],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, d = 1:P2G], unitsretired_P2G[I,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d] - unitsretired_P2G[i0,d]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, d = 1:P2G], unitsretired_P2G[I,d] >= NumUnits_P2G[d]*max(min(Years[I] - RetirementYear_P2G[d],1),0) + sum(unitsbuilt_P2G[i0,d]*max(min(Years[I] - (Years[i0] + Lifetime_P2G[d]),1),0) for i0 = 1:I-1) - sum(unitsretired_P2G[i0,d] for i0 = 1:I-1))
end

@constraint(m, [I = 1, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] <= NumUnits_STORAGE_ELEC[s] + unitsbuilt_STORAGE_ELEC[I,s])
@constraint(m, [I = 1, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] >= NumUnits_STORAGE_ELEC[s]*max(min(Years[I] - RetirementYear_STORAGE_ELEC[s],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] <= NumUnits_STORAGE_ELEC[s] + sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_ELEC], unitsretired_STORAGE_ELEC[I,s] >= NumUnits_STORAGE_ELEC[s]*max(min(Years[I] - RetirementYear_STORAGE_ELEC[s],1),0) + sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min(Years[I] - (Years[i0] + Lifetime_STORAGE_ELEC[s]),1),0) for i0 = 1:I-1)  -  sum(unitsretired_STORAGE_ELEC[i0,s]  for i0 = 1:I-1))
end

### RONDO EDIT
@constraint(m, [I = 1, s = 1:STORAGE_HEAT], unitsretired_STORAGE_HEAT[I,s] <= NumUnits_STORAGE_HEAT[s] + unitsbuilt_STORAGE_HEAT[I,s])
@constraint(m, [I = 1, s = 1:STORAGE_HEAT], unitsretired_STORAGE_HEAT[I,s] >= NumUnits_STORAGE_HEAT[s]*max(min(Years[I] - RetirementYear_STORAGE_HEAT[s],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_HEAT], unitsretired_STORAGE_HEAT[I,s] <= NumUnits_STORAGE_HEAT[s] + sum(unitsbuilt_STORAGE_HEAT[i0,s] - unitsretired_STORAGE_HEAT[i0,s]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_HEAT], unitsretired_STORAGE_HEAT[I,s] >= NumUnits_STORAGE_HEAT[s]*max(min(Years[I] - RetirementYear_STORAGE_HEAT[s],1),0) + sum(unitsbuilt_STORAGE_HEAT[i0,s]*max(min(Years[I] - (Years[i0] + Lifetime_STORAGE_HEAT[s]),1),0) for i0 = 1:I-1)  -  sum(unitsretired_STORAGE_HEAT[i0,s]  for i0 = 1:I-1))
end
#
@constraint(m, [I = 1, d = 1:P2H], unitsretired_P2H[I,d] <= NumUnits_P2H[d] + unitsbuilt_P2H[I,d])
@constraint(m, [I = 1, d = 1:P2H], unitsretired_P2H[I,d] >= NumUnits_P2H[d]*max(min(Years[I] - RetirementYear_P2H[d],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, d = 1:P2H], unitsretired_P2H[I,d] <= NumUnits_P2H[d] + sum(unitsbuilt_P2H[i0,d] - unitsretired_P2H[i0,d]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, d = 1:P2H], unitsretired_P2H[I,d] >= NumUnits_P2H[d]*max(min(Years[I] - RetirementYear_P2H[d],1),0) + sum(unitsbuilt_P2H[i0,d]*max(min(Years[I] - (Years[i0] + Lifetime_P2H[d]),1),0) for i0 = 1:I-1) - sum(unitsretired_P2H[i0,d] for i0 = 1:I-1))
end
### RONDO EDIT

@constraint(m, [I = 1, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] <= NumUnits_STORAGE_GAS[s] + unitsbuilt_STORAGE_GAS[I,s])
@constraint(m, [I = 1, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] >= NumUnits_STORAGE_GAS[s]*max(min(Years[I] - RetirementYear_STORAGE_GAS[s],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] <= NumUnits_STORAGE_GAS[s] + sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, s = 1:STORAGE_GAS], unitsretired_STORAGE_GAS[I,s] >= NumUnits_STORAGE_GAS[s]*max(min(Years[I] - RetirementYear_STORAGE_GAS[s],1),0) + sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min(Years[I] - (Years[i0] + Lifetime_STORAGE_GAS[s]),1),0) for i0 = 1:I-1)  -  sum(unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:I-1))
end


###############################################################################
### Retirement and replacement functions for end-use appliances
# See Eq. 2.10-2.13 in Von Wald thesis
###############################################################################
BaseYear_Sales = zeros(APPLIANCES)
unitsremaining_APPS_historical = zeros(T_inv, APPLIANCES)
for a = 1:APPLIANCES
    BaseYear_Sales[a] = (InitialAppliancePopulation[a]/1000) / sum( max(1 - round(sum(failureProb[a,1:k]),digits = 4),0)/(1+HistoricalGrowthRate[a])^k   for k = 1:50)
    for I = 1:T_inv
        unitsremaining_APPS_historical[I,a] = BaseYear_Sales[a] * sum((1+HistoricalGrowthRate[a])^(-1*k) * (1 - round(sum(failureProb[a,1:k+Int(Years[I] - BaseYear)]),digits = 4))   for k = 1:50)
    end
end

gasApps = cat(collect(1:16), collect(33:48), dims =(1))
gasApps = cat(gasApps, collect(65:80), dims =(1))
gasApps = cat(gasApps, collect(65+32:80+32), dims =(1))
gasApps = cat(gasApps, collect(65+32*2:80+32*2), dims =(1))
GASAPPS = length(gasApps)
unitsremaining_APPS_historical[T_inv,gasApps] .= 0

# See Eq. 2.5 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, a = 1:APPLIANCES], unitsremaining_APPS[I,a] == unitsremaining_APPS_historical[I,a] + sum(unitsbuilt_APPS[i0,a] - unitsretired_APPS[i0,a] for i0 = 1:I))
@constraint(m, [I = 1, a = 1:APPLIANCES], unitsretired_APPS[I,a] == 0)
# @constraint(m, [I = 1:T_inv], sum(unitsbuilt_APPS[I,a] for a = 1:APPLIANCES) <= sum(InitialAppliancePopulation[a]/1000 - unitsremaining_APPS_historical[I,a] for a = 1:APPLIANCES))
@constraint(m, [I = 1:T_inv], sum(unitsremaining_APPS[I,a] for a = 1:APPLIANCES) <= 1.0005*sum(InitialAppliancePopulation[a]/1000 for a = 1:APPLIANCES))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, a = 1:APPLIANCES], unitsretired_APPS[I,a] <= forceretire_multiplier* (sum(round(cumulativefailurefrac[a,v,I],digits = 4)*unitsbuilt_APPS[v,a] for v = 1:I-1) - sum(unitsretired_APPS[i0,a] for i0 = 1:I-1)))
    @constraint(m, [I = 2:T_inv, a = 1:APPLIANCES], unitsretired_APPS[I,a] >= sum(round(cumulativefailurefrac[a,v,I],digits = 4)*unitsbuilt_APPS[v,a] for v = 1:I-1) - sum(unitsretired_APPS[i0,a] for i0 = 1:I-1))
end
if force_retire_gasApps == 1
    if T_inv == 5
        @constraint(m, [I = T_inv, a = gasApps], unitsremaining_APPS[I,a] <= unitsremaining_APPS_historical[I,a])
    end
end

# See Eq. 2.14 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, s = 1:SERVICES, d = 1:DIST_GAS], sum(APP_DistSystemLoc_GAS[d,a]*AppliancesToServices[a,s]*(unitsremaining_APPS[I,a]) for a = 1:APPLIANCES) >= (1+ServicesGrowthRate[s])^(Years[I]-BaseYear)*sum(AppliancesToServices[a,s]*APP_DistSystemLoc_GAS[d,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)/1000)
@constraint(m, [I = 1:T_inv, s = 1:SERVICES, d = 1:DIST_ELEC], sum(APP_DistSystemLoc_ELEC[d,a]*AppliancesToServices[a,s]*(unitsremaining_APPS[I,a]) for a = 1:APPLIANCES) >= (1+ServicesGrowthRate[s])^(Years[I]-BaseYear)*sum(AppliancesToServices[a,s]*APP_DistSystemLoc_ELEC[d,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)/1000)

# See Eq. 2.15 in Von Wald thesis
###############################################################################
if appliance_decisions == 0
    @constraint(m, [I = 1:T_inv, a = 1:APPLIANCES], unitsremaining_APPS[I,a] >= ((1+ForecastGrowthRate[a])^(Years[I]-BaseYear))*InitialAppliancePopulation[a]/1000)
end


###############################################################################
### Retirement and replacement functions for transport classes
# See Eq. 2.5-2.16 in Von Wald thesis
###############################################################################
# COMMENT ON THIS FLEET MODEL:
#   Variables:
#      unitsremaining_Transport_historical[I,tr] ... Number of legacy/historic vehicles (i.e., vehicles present in the initial population) that are still on the road at the end of investment period I (i.e., in year 2025 for I = 1 (assuming 5 investment periods between 2020 and 2045))
#      unitsbuilt_TRANSPORT[I,tr] ... Number of new vehicles built in period I
#      units_retired_TRANSPORT[I,tr] ... Number of vehicles retiring in period I (i.e., retiring legacy vehicles + failing vehicles that have been built in a previous investment period + additional retiring vehicles)
#      sum(cumulativefailurefrac_Transport[tr,v,I]*unitsbuilt_TRANSPORT[v,tr] for v = 1:I-1) ... Number of vehicles build in periods 1 to I-1 that have failed by the end of period I-1.
#      unitsremaining_TRANSPORT ... Total vehicle stock at end of period I
#      unitsbuilt_TRANSPORT_total ... Total number of vehicles added in an investment period (= vehicle sales). Used for implementing ZEV mandate. Not in Van Wald thesis.
#      unitsbuilt_TRANSPORT_zev ... Total number of new ZEV vehicles added in an investment period (= vehicle sales). Used for implementing ZEV mandate. Not in Van Wald thesis.
#   Constraints:
#      Eq. 2.5a: Lower bound for new vehicles per period: 0
#      Eq. 2.5a: Lower bound for retirements per period: 0
#      Eq. 2.5a: Lower bound for total cars in fleet at end of period I: 0
#      Eq. 2.5b: Upper bound for new vehicles per period: currently unbounded
#      Eq. 2.5c: Upper limit for total new vehicles: currently unbounded
#      Eq. 2.5d: Upper bound for retirements in period I>1: Total vehicle stock at beginning of period I (consistent with generator stock modeling above): InitialTransportPopulation + Vehicles built in periods 1:I-1 - Vehicles retired in periods 1:I-1
#      Eq. 2.5d: Upper bound for retirements in period I=1: Total vehicle stock at beginning of period 1 (consistent with modeling of years I>1, not consistent with generator stock modeling which allows retiring of units built in the first period): InitialTransportPopulation
#      Eq. 2.5e: Total number of vehicles at end of period: Initial vehicle population + total vehicles built - total vehicles retired
#      Eq. 2.6:  Lower bound for retirements in period I>1 (consistent with generator stock modeling above): retiring legacy vehicles in period I + vehicles that were built in previous periods and fail in period I - total number of vehicles retired in periods 1:I-1 (the last term accounts conservatively for the fact that some of the previously built vehicles that were set to fail in period I might have already been retired as additional retirementsin in previous years)
#      Eq. 2.6:  Lower bound for retirements in period I=1 (consistent with modeling of years I>1, not consistent with generator stock modeling which allows retiring of units built in the first period): retiring legacy vehicles in period I
#      Added:    Compute total added cars across technologies (used for ZEV mandate scenario)
#      Added:    Compute total added ZEV cars (used for ZEV mandate scenario)
#      Added:    Constraints to implement different ZEV mandate scenarios.

# See Eq. 2.10-2.13 in Von Wald thesis
###############################################################################
BaseYear_Sales_Transport = zeros(TRANSPORTS)
unitsremaining_Transport_historical = zeros(T_inv, TRANSPORTS)
for tr = 1:TRANSPORTS
    BaseYear_Sales_Transport[tr] = (InitialTransportPopulation[tr]/1000) / sum( max(1 - round(sum(failureProb_Transport[tr,1:k]),digits = 4),0)/(1+HistoricalGrowthRate_Transport[tr])^k   for k = 1:50)
    for I = 1:T_inv
        unitsremaining_Transport_historical[I,tr] = BaseYear_Sales_Transport[tr] * sum((1+HistoricalGrowthRate_Transport[tr])^(-1*k) * (1 - round(sum(failureProb_Transport[tr,1:k+Int(Years[I] - BaseYear)]),digits = 4))   for k = 1:50)
    end
end

# See Eq. 2.5-2.6 in Von Wald thesis
###############################################################################
@variable(m, unitsbuilt_TRANSPORT[I = 1:T_inv, tr = 1:TRANSPORTS] >= 0)       # [thousands of units] Eq. 2.5a
@variable(m, unitsremaining_TRANSPORT[I = 1:T_inv, tr = 1:TRANSPORTS] >= 0)   # [thousands of units] Eq. 2.5a
@variable(m, unitsretired_TRANSPORT[I = 1:T_inv, tr = 1:TRANSPORTS] >= 0)     # [thousands of units] Eq. 2.5a

# Eq. 2.5d
@constraint(m, [I = 1, tr = 1:TRANSPORTS], unitsretired_TRANSPORT[I,tr] <= InitialTransportPopulation[tr]/1000)                          
# Eq. 2.6
@constraint(m, [I = 1, tr = 1:TRANSPORTS], unitsretired_TRANSPORT[I,tr] >= InitialTransportPopulation[tr]/1000 - unitsremaining_Transport_historical[I,tr]) 
if T_inv > 1
    # Eq. 2.5d
    @constraint(m, [I = 2:T_inv, tr = 1:TRANSPORTS], unitsretired_TRANSPORT[I,tr] <= InitialTransportPopulation[tr]/1000 + sum(unitsbuilt_TRANSPORT[i0,tr] - unitsretired_TRANSPORT[i0,tr]  for i0 = 1:I-1))
    # Eq. 2.6:                                       retiring vehicles            >= retiring legacy vehicles                                                                  + vehicles build in previous periods and failing in I                                                                                                                                                               - total number of vehicles retiring
    @constraint(m, [I = 2:T_inv, tr = 1:TRANSPORTS], unitsretired_TRANSPORT[I,tr] >= (unitsremaining_Transport_historical[I-1,tr] - unitsremaining_Transport_historical[I,tr]) + (sum(round(cumulativefailurefrac_Transport[tr,v,I],digits = 4)*unitsbuilt_TRANSPORT[v,tr] for v = 1:I) - sum(round(cumulativefailurefrac_Transport[tr,v,I],digits = 4)*unitsbuilt_TRANSPORT[v,tr] for v = 1:I-1)) - sum(unitsretired_TRANSPORT[i0,tr] for i0 = 1:I-1))
end

# Eq. 2.5e
@constraint(m, [I = 1:T_inv, tr = 1:TRANSPORTS], unitsremaining_TRANSPORT[I,tr] == InitialTransportPopulation[tr]/1000 + sum(unitsbuilt_TRANSPORT[i0,tr] - unitsretired_TRANSPORT[i0,tr] for i0 = 1:I))


# See Eq. 2.14 in Von Wald thesis
###############################################################################
# For each investment period, for each service, for each node, ensures that (growing) demand of service (PassMiles) is satisfied. Can be satisfied by any transport type that offers that service at that node.                                                                                                                                                                         
#@constraint(m, [I = 1:T_inv, s = 1:TRANSPORT_SERVICES, n = 1:NODES_GAS], sum(TRANSPORT_NodalLoc_GAS[n,tr]*TransportsToServices[tr,s]*(unitsremaining_TRANSPORT[I,tr]) for tr = 1:TRANSPORTS) >= (1+ServicesGrowthRate_Transport[s])^(Years[I]-BaseYear)*sum(TransportsToServices[tr,s]*TRANSPORT_NodalLoc_GAS[n,a]*InitialTransportPopulation[tr] for tr = 1:TRANSPORTS)/1000)
@constraint(m, [I = 1:T_inv, s = 1:TRANSPORT_SERVICES, n = 1:NODES_ELEC], sum(TRANSPORT_NodalLoc_ELEC[n,tr]*TransportsToServices[tr,s]*(unitsremaining_TRANSPORT[I,tr]) for tr = 1:TRANSPORTS) >= (1+ServicesGrowthRate_Transport[s])^(Years[I]-BaseYear)*sum(TransportsToServices[tr,s]*TRANSPORT_NodalLoc_ELEC[n,tr]*InitialTransportPopulation[tr] for tr = 1:TRANSPORTS)/1000)

# See Eq. 2.15 in Von Wald thesis
###############################################################################
if transport_decisions == 0
    # At the end of each period, the total number of vehicles (within each transport class (no switching)) has to fulill at least the projections
    @constraint(m, [I = 1:T_inv, tr = 1:TRANSPORTS], unitsremaining_TRANSPORT[I,tr] >= ((1+ForecastGrowthRate_Transport[tr])^(Years[I]-BaseYear))*InitialTransportPopulation[tr]/1000)
end

# Scenarios (not in Von Wald thesis)
###############################################################################
@variable(m, unitsbuilt_TRANSPORT_total[I = 1:T_inv] >= 0)                    # [thousands of units] not in thesis
@variable(m, unitsbuilt_TRANSPORT_zev[I = 1:T_inv] >= 0)                      # [thousands of units] not in thesis
@constraint(m, [I = 1:T_inv], unitsbuilt_TRANSPORT_total[I] == sum(unitsbuilt_TRANSPORT[I,tr] for tr = 1:TRANSPORTS)) # computes total number of new vehicles/sales per investment period
ZEV_row_indices = collect(findall(x -> x in ["Veh_Road_LDA_HybridElec", "Veh_Road_LDA_HydrogenFC", "Veh_Road_LDA_Elec"], Transport[:, 6])) # finds rows in TransportNetwork.csv (or the "Transport" DataFrame) that has ZEV technologies.
@constraint(m, [I = 1:T_inv], unitsbuilt_TRANSPORT_zev[I] == sum(unitsbuilt_TRANSPORT[I,tr] for tr = ZEV_row_indices)) # computes total number of new ZEV vehicles/sales per investment period

if transport_scenario_zevmandate == "unconstrained"
    nothing
elseif transport_scenario_zevmandate == "scoping_plan"
    inv_year_indexes_2026 = collect(findall(year -> year >= 2026, Years)) # Finds the "indexes of all investment periods" (T_inv) after a given year.
    println("inv_year_indexes_2026:")
    println(inv_year_indexes_2026)
    inv_year_indexes_2035 = collect(findall(year -> year >= 2035, Years)) # Finds the "indexes of all investment periods" (T_inv) after a given year.
    println("inv_year_indexes_2035:")
    println(inv_year_indexes_2035)
    @constraint(m, [I = inv_year_indexes_2026], unitsbuilt_TRANSPORT_zev[I] >= 0.4 * unitsbuilt_TRANSPORT_total[I])   # 40% ZEV sales by 2026
    @constraint(m, [I = inv_year_indexes_2035], unitsbuilt_TRANSPORT_zev[I] >= 1 * unitsbuilt_TRANSPORT_total[I])     # 100% ZEV sales by 2035
elseif transport_scenario_zevmandate == "scoping_plan_bau"
    inv_year_indexes_2030 = collect(findall(year -> year >= 2030, Years)) # Finds the "indexes of all investment periods" (T_inv) after a given year.
    @constraint(m, [I = inv_year_indexes_2030], unitsbuilt_TRANSPORT_zev[I] >= 0.4 * unitsbuilt_TRANSPORT_total[I])   # 40% ZEV sales by 2030.
elseif transport_scenario_zevmandate == "early"
    inv_year_indexes_2030 = collect(findall(year -> year >= 2030, Years)) # Finds the "indexes of all investment periods" (T_inv) after a given year.
    @constraint(m, [I = inv_year_indexes_2030], unitsbuilt_TRANSPORT_zev[I] >= 1 * unitsbuilt_TRANSPORT_total[I])     # 100% ZEV sales by 2030.
elseif transport_scenario_zevmandate == "late"
    inv_year_indexes_2040 = collect(findall(year -> year >= 2040, Years)) # Finds the "indexes of all investment periods" (T_inv) after a given year.
    @constraint(m, [I = inv_year_indexes_2040], unitsbuilt_TRANSPORT_zev[I] >= 1 * unitsbuilt_TRANSPORT_total[I])     # 100% ZEV sales by 2040.
else
    println("ATTENTION - Unrecognized transport scenario!")
end


if transport_scenario_stockshare == "unconstrained"
    nothing
elseif transport_scenario_stockshare == "early"
    inv_year_indexes_2030 = collect(findall(year -> year >= 2030, Years)) # Finds the "indexes of all investment periods" (T_inv) after a given year.
    @constraint(m, [I = inv_year_indexes_2030], sum(unitsremaining_TRANSPORT[I,tr] for tr = ZEV_row_indices) >= sum(unitsremaining_TRANSPORT[I,tr] for tr = 1:TRANSPORTS))   # 100% ZEV stock by 2030.
elseif transport_scenario_stockshare == "late"
    inv_year_indexes_2045 = collect(findall(year -> year >= 2045, Years)) # Finds the "indexes of all investment periods" (T_inv) after a given year.
    @constraint(m, [I = inv_year_indexes_2045], sum(unitsremaining_TRANSPORT[I,tr] for tr = ZEV_row_indices) <= 0.8 * sum(unitsremaining_TRANSPORT[I,tr] for tr = 1:TRANSPORTS))   # 20% Non-ZEV stock by 2045 (as projected by Madalsa). Unclear if offsets are currently allowed.
else
    println("ATTENTION - Unrecognized transport scenario!")
end

###############################################################################
### Addition of appliances' and transport's incremental demand to baselines
# See Eq. 2.16 in Von Wald thesis
###############################################################################
@variable(m, Demand_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC] >= 0)
@variable(m, Demand_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS] >= 0)
#@variable(m, Demand_LPG[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS] >= 0)

@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], Demand_ELEC[I,T,t,n] == BaselineDemand_ELEC[I,T,t,n] + 1000*sum(APPLIANCES_NodalLoc_ELEC[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_ELEC[T,t,a] for a = 1:APPLIANCES) + 1000*sum(TRANSPORT_NodalLoc_ELEC[n,tr]*(unitsremaining_TRANSPORT[I,tr])*TransportProfiles_ELEC[T,t,tr] for tr = 1:TRANSPORTS))
### RONDO EDIT
# add new variable for baseline demand met by gas: BaselineDemand_fromGas
@variable(m, 0 <= BaselineDemand_fromGAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS])
#
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], Demand_GAS[I,T,t,n] == BaselineDemand_fromGAS[I,T,t,n] + 1000*sum(APPLIANCES_NodalLoc_GAS[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES))

# create similar equations for baseline demand met by clean/direct heat: BaselineDemand_fromDirectHeat
# could add to it later something like hydrogen burning from pipeline or from existing P2G; rather than storage
# could also constrain this by BaselineDemand_HEAT
@variable(m, 0 <= BaselineDemand_fromDirectHeat[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS])
# sum of both BaselineDemands == BaselineDemand_HEAT; which is an input
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], BaselineDemand_HEAT[I,T,t,n] == BaselineDemand_fromGAS[I,T,t,n] + BaselineDemand_fromDirectHeat[I,T,t,n] )
#


# To permit sensitivity testing to disallowing hybrid appliance strategies
if hybrids_allowed == 0
    @constraint(m, [a = 1:APPLIANCES, I = 1:T_inv], unitsbuilt_APPS[I,a] <= (1-IS_HYBRID[a])*10^3)
end
