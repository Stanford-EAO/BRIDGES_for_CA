###############################################################################
### Expansion and retirement of energy supply/demand units
# See Eq. 2.5a/2.5b/2.5c in Von Wald thesis
################################################################################
@variable(m, 0 <= unitsbuilt_GEN[I = 1:T_inv, g = 1:GEN] <= MaxNewUnitsAnnual_GEN[g])                                      # [units]
@variable(m, 0 <= unitsretired_GEN[I = 1:T_inv, g = 1:GEN])                                                                 # [units]
@constraint(m, [g = 1:GEN], sum(unitsbuilt_GEN[i,g] for i = 1:T_inv) <= MaxNewUnitsTotal_GEN[g])                            # [units]

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

# CDR
@variable(m, 0 <= unitsbuilt_CDR[I = 1:T_inv, d = 1:CDR] <= MaxNewUnitsAnnual_CDR[d])                                       # [units]
@variable(m, 0 <= unitsretired_CDR[I = 1:T_inv, d = 1:CDR])                                                                 # [units]
@constraint(m, [d = 1:CDR], sum(unitsbuilt_CDR[i,d] for i = 1:T_inv) <= MaxNewUnitsTotal_CDR[d])                            # [units]


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

# CDR
@constraint(m, [I = 1, d = 1:CDR], unitsretired_CDR[I,d] <= NumUnits_CDR[d] + unitsbuilt_CDR[I,d])
@constraint(m, [I = 1, d = 1:CDR], unitsretired_CDR[I,d] >= NumUnits_CDR[d]*max(min(Years[I] - RetirementYear_CDR[d],1),0))
if T_inv > 1
    @constraint(m, [I = 2:T_inv, d = 1:CDR], unitsretired_CDR[I,d] <= NumUnits_CDR[d] + sum(unitsbuilt_CDR[i0,d] - unitsretired_CDR[i0,d]  for i0 = 1:I-1))
    @constraint(m, [I = 2:T_inv, d = 1:CDR], unitsretired_CDR[I,d] >= NumUnits_CDR[d]*max(min(Years[I] - RetirementYear_CDR[d],1),0) + sum(unitsbuilt_CDR[i0,d]*max(min(Years[I] - (Years[i0] + Lifetime_CDR[d]),1),0) for i0 = 1:I-1) - sum(unitsretired_CDR[i0,d] for i0 = 1:I-1))
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

# See Eq. 2.16 in Von Wald thesis
###############################################################################
@variable(m, Demand_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC] >= 0)
@variable(m, Demand_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS] >= 0)
#@variable(m, Demand_LPG[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS] >= 0)

@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], Demand_ELEC[I,T,t,n] == BaselineDemand_ELEC[I,T,t,n] + 1000*sum(APPLIANCES_NodalLoc_ELEC[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_ELEC[T,t,a] for a = 1:APPLIANCES))
### RONDO EDIT
# add new variable for baseline demand met by gas: BaselineDemand_fromGas
@variable(m, 0 <= BaselineDemand_fromGAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS])
# now for CDR
@variable(m, 0 <= CDR_Demand_fromGAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS])
#
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], Demand_GAS[I,T,t,n] == BaselineDemand_fromGAS[I,T,t,n] + CDR_Demand_fromGAS[I,T,t,n] + 1000*sum(APPLIANCES_NodalLoc_GAS[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES))



# To permit sensitivity testing to disallowing hybrid appliance strategies
if hybrids_allowed == 0
    @constraint(m, [a = 1:APPLIANCES, I = 1:T_inv], unitsbuilt_APPS[I,a] <= (1-IS_HYBRID[a])*10^3)
end