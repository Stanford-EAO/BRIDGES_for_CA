################################################################################
################################################################################
## Optimization program
################################################################################
################################################################################
m = Model(optimizer_with_attributes(Gurobi.Optimizer,"Threads" => 46,"BarHomogeneous" => 1,"ScaleFlag"=>2, "FeasibilityTol"=> 0.005, "OptimalityTol" => 0.001, "BarConvTol"=> 0.0001, "Method"=> 2, "Crossover"=> 0)) # "NumericFocus"=>2, "Presolve"=>2))

###############################################################################
### Expansion and retirement of energy supply/demand units
# See Eq. 2.5a/2.5b/2.5c in Von Wald thesis
################################################################################
@variable(m, 0 <= unitsbuilt_GEN[I = 1:T_inv, g = 1:GEN]  <= MaxNewUnitsAnnual_GEN[g])                                      # [units]
@variable(m, 0 <= unitsretired_GEN[I = 1:T_inv, g = 1:GEN])                                                                 # [units]
@constraint(m, [g = 1:GEN], sum(unitsbuilt_GEN[i,g] for i = 1:T_inv) <= MaxNewUnitsTotal_GEN[g])                            # [units]

@variable(m, 0 <= unitsbuilt_P2G[I = 1:T_inv, d = 1:P2G]  <= MaxNewUnitsAnnual_P2G[d])                                      # [units]
@variable(m, 0 <= unitsretired_P2G[I = 1:T_inv, d = 1:P2G])                                                                 # [units]
@constraint(m, [d = 1:P2G], sum(unitsbuilt_P2G[i,d] for i = 1:T_inv) <= MaxNewUnitsTotal_P2G[d])                            # [units]

@variable(m, 0 <= unitsbuilt_STORAGE_ELEC[I = 1:T_inv, s = 1:STORAGE_ELEC]  <= MaxNewUnitsAnnual_STORAGE_ELEC[s])           # [units]
@variable(m, 0 <= unitsretired_STORAGE_ELEC[I = 1:T_inv, s = 1:STORAGE_ELEC])                                               # [units]
@constraint(m, [s = 1:STORAGE_ELEC], sum(unitsbuilt_STORAGE_ELEC[i,s] for i = 1:T_inv) <= MaxNewUnitsTotal_STORAGE_ELEC[s]) # [units]

@variable(m, 0 <= unitsbuilt_STORAGE_GAS[I = 1:T_inv, s = 1:STORAGE_GAS]  <= MaxNewUnitsAnnual_STORAGE_GAS[s])              # [units]
@variable(m, 0 <= unitsretired_STORAGE_GAS[I = 1:T_inv, s = 1:STORAGE_GAS])                                                 # [units]
@constraint(m, [s = 1:STORAGE_GAS], sum(unitsbuilt_STORAGE_GAS[i,s] for i = 1:T_inv) <= MaxNewUnitsTotal_STORAGE_GAS[s])    # [units]

@variable(m, unitsbuilt_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)       # [thousands of units]
@variable(m, unitsremaining_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)   # [thousands of units]
@variable(m, unitsretired_APPS[I = 1:T_inv, a = 1:APPLIANCES] >= 0)     # [thousands of units]


###############################################################################
### Retirement functions for generators, p2g, and storage units
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

# See Eq. 2.5 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, a = 1:APPLIANCES], unitsremaining_APPS[I,a] == unitsremaining_APPS_historical[I,a] + sum(unitsbuilt_APPS[i0,a] - unitsretired_APPS[i0,a] for i0 = 1:I))

@constraint(m, [I = 1, a = 1:APPLIANCES], unitsretired_APPS[I,a] <= forceretire_multiplier * (InitialAppliancePopulation[a]/1000 - unitsremaining_APPS_historical[I,a]))

if T_inv > 1
    @constraint(m, [I = 2:T_inv, a = 1:APPLIANCES], unitsretired_APPS[I,a] <= forceretire_multiplier* (sum(round(cumulativefailurefrac[a,v,I],digits = 4)*unitsbuilt_APPS[v,a] for v = 1:I-1) - sum(unitsretired_APPS[i0,a] for i0 = 1:I-1)))
    @constraint(m, [I = 2:T_inv, a = 1:APPLIANCES], unitsretired_APPS[I,a] >= sum(round(cumulativefailurefrac[a,v,I],digits = 4)*unitsbuilt_APPS[v,a] for v = 1:I-1) - sum(unitsretired_APPS[i0,a] for i0 = 1:I-1))
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
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], Demand_GAS[I,T,t,n] == BaselineDemand_GAS[I,T,t,n] + 1000*sum(APPLIANCES_NodalLoc_GAS[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES))
#@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], Demand_LPG[I,T,t,n] == 1000*sum(APPLIANCES_NodalLoc_GAS[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_LPG[T,t,a] for a = 1:APPLIANCES))

# To permit sensitivity testing to disallowing hybrid appliance strategies
if hybrids_allowed == 0
    @constraint(m, [a = 1:APPLIANCES, I = 1:T_inv], unitsbuilt_APPS[I,a] <= (1-IS_HYBRID[a])*10^3)
end


###############################################################################
### Dispatch operations of electric generators and gas production
# See Eq. 2.22a in Von Wald thesis
################################################################################
@variable(m, generation[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)         # [MWh]
@variable(m, commit_GEN[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)         # [no. units]
@variable(m, startup_GEN[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)        # [no. units]
@variable(m, shutdown_GEN[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)       # [no. units]

@variable(m, P2G_dispatch[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)       # [MWh]
@variable(m, commit_P2G[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)         # [no. units]
@variable(m, startup_P2G[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)        # [no. units]
@variable(m, shutdown_P2G[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G] >= 0)       # [no. units]


### Startup and shutdown events
# See Eq. 2.22a/2.22e in Von Wald thesis
################################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], commit_GEN[I,T,t,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], startup_GEN[I,T,t,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], shutdown_GEN[I,T,t,g] <= NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, g = 1:GEN], commit_GEN[I,T,t,g] == commit_GEN[I,T,t-1,g] + startup_GEN[I,T,t,g] - shutdown_GEN[I,T,t,g])

@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], commit_P2G[I,T,t,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i,d] - unitsretired_P2G[i,d] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], startup_P2G[I,T,t,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i,d] - unitsretired_P2G[i,d] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], shutdown_P2G[I,T,t,d] <= NumUnits_P2G[d] + sum(unitsbuilt_P2G[i,d] - unitsretired_P2G[i,d] for i = 1:I))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, d = 1:P2G], commit_P2G[I,T,t,d] == commit_P2G[I,T,t-1,d] + startup_P2G[I,T,t,d] - shutdown_P2G[I,T,t,d])

### Min and Max generation constraints
# See Eq. 2.22b/2.22c in Von Wald thesis
################################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], generation[I,T,t,g] >= Pmin_GEN[g]*UnitSize_GEN[g]*commit_GEN[I,T,t,g])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], generation[I,T,t,g] <= Pmax_GEN[g]*UnitSize_GEN[g]*commit_GEN[I,T,t,g])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], P2G_dispatch[I,T,t,d] >= Pmin_P2G[d]*UnitSize_P2G[d]*commit_P2G[I,T,t,d])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, d = 1:P2G], P2G_dispatch[I,T,t,d] <= Pmax_P2G[d]*UnitSize_P2G[d]*commit_P2G[I,T,t,d])

### Constraints on fixed profile generation resources
# See Eq. 2.22d in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], generation[I,T,t,g] <= HourlyVRE[T,t,g]*UnitSize_GEN[g]*(NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I))) # allows for curtailment
# Explicit calculation of renewable energy curtailments
@variable(m, curtailmentRE[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN] >= 0)
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, g = 1:GEN], curtailmentRE[I,T,t,g] == IS_RENEWABLE[g]*(HourlyVRE[T,t,g]*UnitSize_GEN[g]*(NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I)) - generation[I,T,t,g])) # allows for curtailment

### Ramping constraint
# See Eq. 2.23 in Von Wald thesis
###############################################################################
if t_ops > 1
    minimax = min.(Pmax_GEN, max.(Pmin_GEN,RampDownRate_GEN))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, g = 1:GEN], generation[I,T,t-1,g]-generation[I,T,t,g] <= RampDownRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,T,t,g]-startup_GEN[I,T,t,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*startup_GEN[I,T,t,g] + minimax[g]*UnitSize_GEN[g]*shutdown_GEN[I,T,t,g])
    minimax = min.(Pmax_GEN, max.(Pmin_GEN,RampUpRate_GEN))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 2:t_ops, g = 1:GEN], generation[I,T,t,g]-generation[I,T,t-1,g] <= RampUpRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,T,t,g]-startup_GEN[I,T,t,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*shutdown_GEN[I,T,t,g] + minimax[g]*UnitSize_GEN[g]*startup_GEN[I,T,t,g])
end

### Min up time/down time constraints
# See Eq. 2.25/2.26 in Von Wald thesis
###############################################################################
for g = 1:GEN
    if t_ops > MinUpTime_GEN[g]
        for t = Int(MinUpTime_GEN[g]+1):t_ops
            @constraint(m, [I = 1:T_inv, T = 1:T_ops], commit_GEN[I,T,t,g] >= sum(startup_GEN[I,T,t0,g] for t0 = (t-Int(MinUpTime_GEN[g])):t))
        end
    end
end
for g = 1:GEN
    if t_ops > MinDownTime_GEN[g]
        for t = Int(MinDownTime_GEN[g]+1):t_ops
            @constraint(m, [I = 1:T_inv, T = 1:T_ops], NumUnits_GEN[g] + sum(unitsbuilt_GEN[i,g] - unitsretired_GEN[i,g] for i = 1:I) - commit_GEN[I,T,t,g] >= sum(shutdown_GEN[I,T,t0,g] for t0 = (t-MinDownTime_GEN[g]):t))
        end
    end
end

### Generator operational constraints across linked time periods
# See Eq. 2.24/2.25/2.26 in Von Wald thesis
###############################################################################
if LINKED_PERIODS_GENOPS == 1
    for g = 1:GEN
        for i = 2:Int(Periods_Per_Year)
        # Ramping constraint
            minmax = min(Pmax_GEN[g], max(Pmin_GEN[g],RampDownRate_GEN[g]))
            @constraint(m, [I = 1:T_inv], generation[I,Int(RepDays[I,i-1]),t_ops,g]-generation[I,Int(RepDays[I,i]),1,g] <= RampDownRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,Int(RepDays[I,i]),1,g]-startup_GEN[I,Int(RepDays[I,i]),1,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*startup_GEN[I,Int(RepDays[I,i]),1,g] + minmax*UnitSize_GEN[g]*shutdown_GEN[I,Int(RepDays[I,i]),1,g])
            minmax = min(Pmax_GEN[g], max(Pmin_GEN[g],RampUpRate_GEN[g]))
            @constraint(m, [I = 1:T_inv], generation[I,Int(RepDays[I,i]),1,g]-generation[I,Int(RepDays[I,i-1]),t_ops,g] <= RampUpRate_GEN[g]*UnitSize_GEN[g]*(commit_GEN[I,Int(RepDays[I,i]),1,g]-startup_GEN[I,Int(RepDays[I,i]),1,g]) - Pmin_GEN[g]*UnitSize_GEN[g]*shutdown_GEN[I,Int(RepDays[I,i]),1,g] + minmax*UnitSize_GEN[g]*startup_GEN[I,Int(RepDays[I,i]),1,g])
        end
        # Min up time/down time constraints
        firstpd_up = Int(round(MinUpTime_GEN[g]/t_ops)+2)
        X = repeat(collect(range(t_ops,stop = 1,length = t_ops)), outer = [firstpd_up+3])
        Y = repeat(collect(range(0, stop = 24, length = 25)), inner = [t_ops])
        for i = firstpd_up:Int(Periods_Per_Year)
            for t = 1:t_ops
                @constraint(m, [I = 1:T_inv], commit_GEN[I,Int(RepDays[I,i]),t,g] >= sum(startup_GEN[I,Int(RepDays[I,Int(i-Y[Int(t0-t+t_ops)])]),Int(X[Int(t0+t_ops-t)]),g] for t0 = 1:MinUpTime_GEN[g]))
            end
        end
        firstpd_down = Int(round(MinDownTime_GEN[g]/t_ops)+2)
        X = repeat(collect(range(t_ops,stop = 1,length = t_ops)), outer = [firstpd_up+3])
        Y = repeat(collect(range(0, stop = 24, length = 25)), inner = [t_ops])
        for i = firstpd_down:Int(Periods_Per_Year)
            for t = 1:HOURS_PER_PERIOD
                @constraint(m, [I = 1:T_inv], NumUnits_GEN[g] + sum(unitsbuilt_GEN[i0,g] - unitsretired_GEN[i0,g] for i0 = 1:I) - commit_GEN[I,Int(RepDays[I,i]),t,g] >= sum(shutdown_GEN[I,Int(RepDays[I,Int(i-Y[Int(t0-t+t_ops)])]),Int(X[Int(t0+t_ops-t)]),g] for t0 = 1:MinDownTime_GEN[g]))
            end
        end
    end
end


###############################################################################
### Dispatch of storage
# See Eq. 2.51-2.55 in Von Wald thesis
###############################################################################
@variable(m, storedEnergy_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_ELEC] >= 0)       # [MWh]
@variable(m, charging_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC] >= 0)             # [MW]
@variable(m, discharging_ELEC[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC] >= 0)          # [MW]

@variable(m, storedEnergy_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_GAS] >= 0)         # [MWh]
@variable(m, charging_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS] >= 0)               # [MW]
@variable(m, discharging_GAS[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS] >= 0)            # [MW]


@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], storedEnergy_ELEC[I,T,t+1,s] == storedEnergy_ELEC[I,T,t,s] + eta_charging_ELEC[s]*charging_ELEC[I,T,t,s] - (1/eta_discharging_ELEC[s])*discharging_ELEC[I,T,t,s]-eta_loss_ELEC[s]*storedEnergy_ELEC[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_ELEC], storedEnergy_ELEC[I,T,t,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], charging_ELEC[I,T,t,s] <= (1/eta_charging_ELEC[s])*UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], charging_ELEC[I,T,t,s] <= duration_ELEC[s]*UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)) - storedEnergy_ELEC[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], discharging_ELEC[I,T,t,s] <= eta_discharging_ELEC[s]*UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], discharging_ELEC[I,T,t,s] <= storedEnergy_ELEC[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], discharging_ELEC[I,T,t,s]/eta_discharging_ELEC[s] + eta_charging_ELEC[s]*charging_ELEC[I,T,t,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I)))

@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,t+1,s] == storedEnergy_GAS[I,T,t,s] + eta_charging_GAS[s]*charging_GAS[I,T,t,s] - (1/eta_discharging_GAS[s])*discharging_GAS[I,T,t,s]-eta_loss_GAS[s]*storedEnergy_GAS[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops+1, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,t,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I))*duration_GAS[s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], charging_GAS[I,T,t,s] <= (1/eta_charging_GAS[s])*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], charging_GAS[I,T,t,s] <= duration_GAS[s]*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)) - storedEnergy_GAS[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s] <= eta_discharging_GAS[s]*UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s] <= storedEnergy_GAS[I,T,t,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s]/eta_discharging_GAS[s] + eta_charging_GAS[s]*charging_GAS[I,T,t,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s] for i = 1:I)))

## To limit flexible charge/discharge associated with gas storage
# See Eq. 2.52 in Von Wald thesis
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops-1, s = 1:STORAGE_GAS], charging_GAS[I,T,t,s] == charging_GAS[I,T,t+1,s])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops-1, s = 1:STORAGE_GAS], discharging_GAS[I,T,t,s] == discharging_GAS[I,T,t+1,s])


### Constraints for linking energy storage across representative periods
# See Eq. 2.56-2.59 in Von Wald thesis
################################################################################
if LINKED_PERIODS_STORAGE == 0
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC], storedEnergy_ELEC[I,T,1,s] == storedEnergy_ELEC[I,T,t_ops+1,s])       # [MWh]
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS], storedEnergy_GAS[I,T,1,s] == storedEnergy_GAS[I,T,t_ops+1,s])       # [MWh]
end
if LINKED_PERIODS_STORAGE == 1
    @variable(m, MinSOC_ELEC[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC] >= 0)
    @variable(m, MaxSOC_ELEC[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC] >= 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC], MinSOC_ELEC[I,T,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_ELEC], MaxSOC_ELEC[I,T,s] <= UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
    @variable(m, SOCTracked_ELEC[I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_ELEC] >= 0)
    @constraint(m, [I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,d,s] <=  UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i,s] - unitsretired_STORAGE_ELEC[i,s] for i = 1:I))*duration_ELEC[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], MinSOC_ELEC[I,T,s] <= storedEnergy_ELEC[I,T,t,s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_ELEC], MaxSOC_ELEC[I,T,s] >= storedEnergy_ELEC[I,T,t,s])
    for i = 1:Int(Periods_Per_Year)
        @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,i,s] == storedEnergy_ELEC[I,Int(RepDays[I,1]),1,s] + sum(storedEnergy_ELEC[I,Int(RepDays[I,t]),t_ops+1,s] - storedEnergy_ELEC[I,Int(RepDays[I,t]),1,s] for t = 1:i))
        if i < Periods_Per_Year
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,i,s] + (MaxSOC_ELEC[I,Int(RepDays[I,i+1]),s] - storedEnergy_ELEC[I,Int(RepDays[I,i+1]),1,s]) <=  UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:I))*duration_ELEC[s])
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,i,s] - (storedEnergy_ELEC[I,Int(RepDays[I,i+1]),1,s] - MinSOC_ELEC[I,Int(RepDays[I,i+1]),s]) >= 0)
        end
    end
    @constraint(m, [I = 1:T_inv, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,Int(Periods_Per_Year),s] == storedEnergy_ELEC[I,Int(RepDays[I,1]),1,s])
    
    @variable(m, MinSOC_GAS[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS] >= 0)
    @variable(m, MaxSOC_GAS[I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS] >= 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS], MinSOC_GAS[I,T,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I))*duration_GAS[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, s = 1:STORAGE_GAS], MaxSOC_GAS[I,T,s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I))*duration_GAS[s])
    @variable(m, SOCTracked_GAS[I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_GAS] >= 0)
    @constraint(m, [I = 1:T_inv, d = 1:Int(Periods_Per_Year), s = 1:STORAGE_GAS], SOCTracked_GAS[I, d, s] <= UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i,s] - unitsretired_STORAGE_GAS[i,s]  for i = 1:I))*duration_GAS[s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], MinSOC_GAS[I,T,s] <= storedEnergy_GAS[I,T,t,s])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, s = 1:STORAGE_GAS], MaxSOC_GAS[I,T,s] >= storedEnergy_GAS[I,T,t,s])
    
    for i = 1:Int(Periods_Per_Year)
        @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,i,s] == storedEnergy_GAS[I,Int(RepDays[I,1]),1,s] + sum(storedEnergy_GAS[I,Int(RepDays[I,t]),t_ops+1,s] - storedEnergy_GAS[I,Int(RepDays[I,t]),1,s] for t = 1:i))
        if i < Periods_Per_Year
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,i,s] + (MaxSOC_GAS[I,Int(RepDays[I,i+1]),s] - storedEnergy_GAS[I,Int(RepDays[I,i+1]),1,s]) <=  UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s] for i0 = 1:I))*duration_GAS[s])
            @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,i,s] - (storedEnergy_GAS[I,Int(RepDays[I,i+1]),1,s] - MinSOC_GAS[I,Int(RepDays[I,i+1]),s]) >= 0)
        end
    end
    @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,Int(Periods_Per_Year),s] == storedEnergy_GAS[I,Int(RepDays[I,1]),1,s])

    # Set initial condition for gas storage based on specified in input file
    # @constraint(m, [I = 1, t = 1, s = 1:STORAGE_GAS], storedEnergy_GAS[I,Int(RepDays[I,t]),t,s] == initialStorage_GAS[s])
    @constraint(m, [I = 1:T_inv, d = 1, s = 1:STORAGE_GAS], SOCTracked_GAS[I,d,s] == initialStorage_GAS[s])

end


###############################################################################
### Transmission of electric power
# See Eq. 2.18 in Von Wald thesis
###############################################################################
@variable(m, Flows_Elec[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC])

# Built to permit eventual inclusion of transmission expansion
if TRANSMISSION_EXPANSION_ELEC == 0
    unitsbuilt_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
    unitsretired_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
    addflow_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
elseif TRANSMISSION_EXPANSION_ELEC == 1
    @variable(m, unitsbuilt_TRANS_ELEC[I = 1:T_inv, e = 1:EDGES_ELEC], Bin)
    @variable(m, unitsretired_TRANS_ELEC[I = 1:T_inv, e = 1:EDGES_ELEC], Bin)
    @constraint(m, [e = 1:EDGES_ELEC], sum(unitsbuilt_TRANS_ELEC[I,e] for I = 1:T_inv) <=  MaxNewUnits_ElecTrans[e])
    @constraint(m, [I = 1:T_inv, e = 1:EDGES_ELEC], sum(unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) <=  ExistingUnits_ElecTrans[e] + sum(unitsbuilt_TRANS_ELEC[i0,e] for i0 = 1:I))
    addflow_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
elseif TRANSMISSION_EXPANSION_ELEC == 2
    unitsbuilt_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
    unitsretired_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
    @variable(m, addflow_TRANS_ELEC[I = 1:T_inv,e = 1:EDGES_ELEC] >= 0)        # MW
end

## If not simulating steady-state power flows, flows are only bound by the maximum power flow specified for the line (i.e., simple transport model)
if STEADYSTATE_ELEC == 0
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] <= (sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*MAXFLOW_ELEC[e] + sum(addflow_TRANS_ELEC[i0,e] for i0 = 1:I))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] >= -1*(sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*MAXFLOW_ELEC[e] - sum(addflow_TRANS_ELEC[i0,e] for i0 = 1:I))
end

## If simulating steady-state power flows, we introduce additional voltage angle variables and constraints governing the power flow
if STEADYSTATE_ELEC == 1
    SLACK_BUS = 1
    @variable(m, -2*pi <= BusAngle[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC] <= 2*pi)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops], BusAngle[I,T,t,SLACK_BUS] == 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] <= BASEMVA*-1/Line_Reactance[e]*sum(A_ELEC[n,e]*BusAngle[I,T,t,n] for n = 1:NODES_ELEC) + MAXFLOW_ELEC[e]*(1-sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e]) + sum(addflow_TRANS_ELEC[i0,e] for i0 = 1:I))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] >= BASEMVA*-1/Line_Reactance[e]*sum(A_ELEC[n,e]*BusAngle[I,T,t,n] for n = 1:NODES_ELEC) - MAXFLOW_ELEC[e]*(1-sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e]) - sum(addflow_TRANS_ELEC[i0,e] for i0 = 1:I))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] <= (sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*Line_Rating[e] + sum(addflow_TRANS_ELEC[i0,e] for i0 = 1:I))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, e = 1:EDGES_ELEC], Flows_Elec[I,T,t,e] >= -1*(sum(unitsbuilt_TRANS_ELEC[i0,e] - unitsretired_TRANS_ELEC[i0,e] for i0 = 1:I) + ExistingUnits_ElecTrans[e])*Line_Rating[e] - sum(addflow_TRANS_ELEC[i0,e] for i0 = 1:I))
end

###############################################################################
### Energy Balance for electricity grid
# See Eq. 2.19 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], sum(GEN_NodalLoc_ELEC[n,g]*generation[I,T,t,g] for g = 1:GEN) + sum(-1*A_ELEC[n,e]*Flows_Elec[I,T,t,e] for e = 1:EDGES_ELEC) - sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(charging_ELEC[I,T,t,s]-discharging_ELEC[I,T,t,s]) for s = 1:STORAGE_ELEC) - Demand_ELEC[I,T,t,n] - sum(P2G_NodalLoc_ELEC[n,d]*P2G_dispatch[I,T,t,d]*(1-ISBIOMETHANE[d]) for d = 1:P2G) >= 0)


###############################################################################
### Transmission of gaseous fuel
# See Eq. 2.X in Von Wald thesis
###############################################################################
# Nodal gaseous energy demand balance imposed on average across the day
@variable(m, Flows_Gas[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS])                      # [standard m3/sec]

if TRANSMISSION_EXPANSION_GAS == 0
    unitsbuilt_TRANS_GAS = zeros(T_inv,EDGES_GAS)
    unitsretired_TRANS_GAS = zeros(T_inv,EDGES_GAS)
end
if TRANSMISSION_EXPANSION_GAS == 1
    @variable(m, unitsbuilt_TRANS_GAS[I = 1:T_inv, e = 1:EDGES_GAS], Bin)
    @variable(m, unitsretired_TRANS_GAS[I = 1:T_inv, e = 1:EDGES_GAS], Bin)
    @constraint(m, [e = 1:EDGES_GAS], sum(unitsbuilt_TRANS_GAS[I,e] for I = 1:T_inv) <=  MaxNewUnits_GasTrans[e])
    @constraint(m, [I = 1:T_inv, e = 1:EDGES_GAS], sum(unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) <=  ExistingUnits_GasTrans[e] + sum(unitsbuilt_TRANS_GAS[i0,e] for i0 = 1:I))
end

if GASFLOW_DIRECTIONS == 0
    GasFlowDirection = ones(T_inv,T_ops,EDGES_GAS)
end
if GASFLOW_DIRECTIONS == 1
    @variable(m, GasFlowDirection[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Bin)
end

@constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], -1*Flows_Gas[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
@constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e] <= GasFlowDirection[I,T,e]*MAXFLOW_GAS[e])
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e] >= (GasFlowDirection[I,T,e]-1)*(MAXFLOW_GAS[e]))

if (STEADYSTATE_GAS == 1) & (NODES_GAS > 1)
    @variable(m, (PRESSURE_MIN) <= NodalPressureSqrd[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] <= (PRESSURE_MAX))
    @variable(m, PRESSURE_MIN <= CompressionPressureSqrd[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd[I,T,e] >= sum(max(A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd[I,T,e] <= CompressionRatio_MAX_Branch[e]*sum(max(A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS))

    @constraint(m, [I = 1:T_inv, T = 1:T_ops], NodalPressureSqrd[I,T,SLACK_NODE] == (SLACK_NODE_Pressure))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS],  CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS) <= GasFlowDirection[I,T,e]*((PRESSURE_MAX)-(PRESSURE_MIN)))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS) >= (1-GasFlowDirection[I,T,e])*((PRESSURE_MIN)-(PRESSURE_MAX)))
    
    
    @variable(m, 0 <= lambda[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX-PRESSURE_MIN)
    a1 = -1
    a2 = PRESSURE_MIN-(PRESSURE_MAX)
    b1 = 1
    b2 = (PRESSURE_MAX)-(PRESSURE_MIN)

    # Lambda is a tight constraint such that when y = 1, lambda = P1-P2 and when y = 0, lambda = P2-P1
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] >= a2*(2*GasFlowDirection[I,T,e]-1) + a1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - a1*a2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] >= b2*(2*GasFlowDirection[I,T,e]-1) + b1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - b1*b2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] <= b2*(2*GasFlowDirection[I,T,e]-1) + a1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - a1*b2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda[I,T,e] <= a2*(2*GasFlowDirection[I,T,e]-1) + b1*(CompressionPressureSqrd[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS)) - b1*a2)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas[I,T,e]*Flows_Gas[I,T,e] <= K2[e]/Length_Pipes[e]*(lambda[I,T,e]))
end


################################################################################
### Gaseous fuel slack supply
# See Eq. 2.41 in Von Wald thesis
###############################################################################
@variable(m, SUPPLY_GAS_slack[I = 1:T_inv, T = 1:T_ops,  n = 1:NODES_GAS] >= 0)
@constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], SUPPLY_GAS_slack[I,T,n] <= MAXSLACK[n])

## Introduce nominal gas component tracking formulation
# See Eq. 2.42 in Von Wald thesis
###############################################################################
@variable(m, NominalGasOfftakes[I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS, g = 1:GAS_COMPONENTS] >= 0)     # kmol/sec
@variable(m, NominalGasFlows[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS])                          # kmol/sec
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] <= GasFlowDirection[I,T,e]*(MAXFLOW_GAS[e]))
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] >= (GasFlowDirection[I,T,e]-1)*(MAXFLOW_GAS[e]))

## Constrain nominal values to match physical values
# See Eq. 2.43-2.44 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], sum(MolarMass[g]*NominalGasFlows[I,T,e,g] for g = 1:GAS_COMPONENTS) == M_CH4*V_m*Flows_Gas[I,T,e])
@constraint(m, [I = 1:T_inv,  T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], sum(LHV[g]*MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) == Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN))

## Molar balance
# See Eq. 2.45 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS, g = 1:GAS_COMPONENTS], MoleFracs_SLACK[n,g]/(LHV_SLACK[n]*MolarMass_SLACK[n])*SUPPLY_GAS_slack[I,T,n] + sum(-1*A_GAS[n,e]*NominalGasFlows[I,T,e,g] for e = 1:EDGES_GAS) - sum(NominalGasOfftakes[I,T,t,n,g] for t = 1:t_ops)/t_ops - sum(sum(MoleFracs_STORAGE[s,g]/(LHV_STORAGE[s]*MolarMass_STORAGE[s])*STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) for t = 1:t_ops)/t_ops + sum(sum(MoleFracs_P2G[d,g]/(LHV_P2G[d]*MolarMass_P2G[d])*P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)

## Energy balance
# See Eq. 2.47 in Von Wald thesis
###############################################################################
@constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], sum(sum(-1*A_GAS[n,e]*LHV[g]*MolarMass[g]*NominalGasFlows[I,T,e,g] for e = 1:EDGES_GAS) for g = 1:GAS_COMPONENTS) + SUPPLY_GAS_slack[I,T,n]  - sum(sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) for t = 1:t_ops)/t_ops - sum(Demand_GAS[I,T,t,n] + sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) - sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops)/t_ops == 0)


################################################################################
### Gas quality constraints
# See Eq. 2.48-2.50 in Von Wald thesis
################################################################################
if GasQuality == "Annual"
    ## Mole fraction limit imposed on annual, system-wide basis
    # See Eq. 2.48 in Von Wald thesis
    @constraint(m, [I = 1:T_inv, g = 1:GAS_COMPONENTS], sum(weights[i,T]*8760/t_ops*sum(sum(NominalGasOfftakes[I,T,t,n,g] for t = 1:t_ops) for n = 1:NODES_GAS) for T = 1:T_ops) <= MoleFrac_MAX[g]*sum(sum(weights[i,T]*8760/t_ops*sum(sum(NominalGasOfftakes[I,T,t,n,h] for t = 1:t_ops) for n = 1:NODES_GAS) for T = 1:T_ops) for h = 1:GAS_COMPONENTS))
end
if GasQuality == "Nodal"
    ## Mole fraction limit imposed on daily, nodal basis
    # See Eq. 2.49 in Von Wald thesis
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] <= MoleFrac_MAX[g]*sum(NominalGasFlows[I,T,e,h] for h = 1:GAS_COMPONENTS) + (1-GasFlowDirection[I,T,e])*(MAXFLOW_GAS[e]))
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS, g = 1:GAS_COMPONENTS], NominalGasFlows[I,T,e,g] >= MoleFrac_MAX[g]*sum(NominalGasFlows[I,T,e,h] for h = 1:GAS_COMPONENTS) - GasFlowDirection[I,T,e]*(MAXFLOW_GAS[e]))
    
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, n = 1:NODES_GAS, g = 1:GAS_COMPONENTS], sum(NominalGasOfftakes[I,T,t,n,g] for t = 1:t_ops) <= MoleFrac_MAX[g]*sum(sum(NominalGasOfftakes[I,T,t,n,h] for t = 1:t_ops) for h = 1:GAS_COMPONENTS))

    
    ## Heating value limit imposed on daily, nodal basis
    # See Eq. 2.50 in Von Wald thesis
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, n = 1:NODES_GAS], sum(sum(LHV[g]*MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops) <= HV_MAX*sum(sum(MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops))
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, n = 1:NODES_GAS], sum(sum(LHV[g]*MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops) >= HV_MIN*sum(sum(MolarMass[g]*NominalGasOfftakes[I,T,t,n,g] for g = 1:GAS_COMPONENTS) for t = 1:t_ops))
end


###############################################################################
### Bounding steady-states
# Not included in Von Wald thesis
###############################################################################
if bounding_steady_states == 1
    @variable(m, Flows_Gas_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS])                  # [m3/sec]
    @variable(m, Flows_Gas_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS])                  # [m3/sec]
    @variable(m, SUPPLY_GAS_slack_Max[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] >= 0)
    @variable(m, SUPPLY_GAS_slack_Min[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] >= 0)
    @variable(m, NetGasDemands_Max[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS])
    @variable(m, NetGasDemands_Min[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS])
    
    
    if GASFLOW_DIRECTIONS == 0
        GasFlowDirection_Min = ones(T_inv,T_ops,EDGES_GAS)
        GasFlowDirection_Max = ones(T_inv,T_ops,EDGES_GAS)
    end
    if GASFLOW_DIRECTIONS == 1
        @variable(m, GasFlowDirection_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Bin)
        @variable(m, GasFlowDirection_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Bin)
    end
    
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], -1*Flows_Gas_Min[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e] <= GasFlowDirection_Min[I,T,e]*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e] >= (GasFlowDirection_Min[I,T,e]-1)*(MAXFLOW_GAS[e]))
    
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], -1*Flows_Gas_Max[I,T,e] <= (sum(unitsbuilt_TRANS_GAS[i0,e] - unitsretired_TRANS_GAS[i0,e] for i0 = 1:I) + ExistingUnits_GasTrans[e])*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e] <= GasFlowDirection_Max[I,T,e]*MAXFLOW_GAS[e])
    @constraint(m, [I = 1:T_inv,  T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e] >= (GasFlowDirection_Max[I,T,e]-1)*(MAXFLOW_GAS[e]))
    
    if (STEADYSTATE_GAS == 1) & (NODES_GAS > 1)
        @variable(m, (PRESSURE_MIN) <= NodalPressureSqrd_Min[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] <= (PRESSURE_MAX))
        @variable(m, PRESSURE_MIN <= CompressionPressureSqrd_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX)
        @variable(m, (PRESSURE_MIN) <= NodalPressureSqrd_Max[I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS] <= (PRESSURE_MAX))
        @variable(m, PRESSURE_MIN <= CompressionPressureSqrd_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Min[I,T,e] >= sum(max(A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Min[I,T,e] <= CompressionRatio_MAX_Branch[e]*sum(max(A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Max[I,T,e] >= sum(max(A_GAS[n,e],0)*NodalPressureSqrd[I,T,n] for n = 1:NODES_GAS))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Max[I,T,e] <= CompressionRatio_MAX_Branch[e]*sum(max(A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS))
    
        @constraint(m, [I = 1:T_inv, T = 1:T_ops], NodalPressureSqrd_Min[I,T,SLACK_NODE] == (SLACK_NODE_Pressure))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS],  CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS) <= GasFlowDirection_Min[I,T,e]*((PRESSURE_MAX)-(PRESSURE_MIN)))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS) >= (1-GasFlowDirection_Min[I,T,e])*((PRESSURE_MIN)-(PRESSURE_MAX)))
    
        @constraint(m, [I = 1:T_inv, T = 1:T_ops], NodalPressureSqrd_Max[I,T,SLACK_NODE] == (SLACK_NODE_Pressure))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS],  CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS) <= GasFlowDirection_Max[I,T,e]*((PRESSURE_MAX)-(PRESSURE_MIN)))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS) >= (1-GasFlowDirection_Max[I,T,e])*((PRESSURE_MIN)-(PRESSURE_MAX)))
        @variable(m, 0 <= lambda_Min[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX-PRESSURE_MIN)
        @variable(m, 0 <= lambda_Max[I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS] <= PRESSURE_MAX-PRESSURE_MIN)
        a1 = -1
        a2 = PRESSURE_MIN-(PRESSURE_MAX)
        b1 = 1
        b2 = (PRESSURE_MAX)-(PRESSURE_MIN)
    
        # Lambda is a tight constraint such that when y = 1, lambda = P1-P2 and when y = 0, lambda = P2-P1
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] >= a2*(2*GasFlowDirection_Min[I,T,e]-1) + a1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - a1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] >= b2*(2*GasFlowDirection_Min[I,T,e]-1) + b1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - b1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] <= b2*(2*GasFlowDirection_Min[I,T,e]-1) + a1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - a1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Min[I,T,e] <= a2*(2*GasFlowDirection_Min[I,T,e]-1) + b1*(CompressionPressureSqrd_Min[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Min[I,T,n] for n = 1:NODES_GAS)) - b1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] >= a2*(2*GasFlowDirection_Max[I,T,e]-1) + a1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - a1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] >= b2*(2*GasFlowDirection_Max[I,T,e]-1) + b1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - b1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] <= b2*(2*GasFlowDirection_Max[I,T,e]-1) + a1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - a1*b2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], lambda_Max[I,T,e] <= a2*(2*GasFlowDirection_Max[I,T,e]-1) + b1*(CompressionPressureSqrd_Max[I,T,e]-sum(max(-1*A_GAS[n,e],0)*NodalPressureSqrd_Max[I,T,n] for n = 1:NODES_GAS)) - b1*a2)
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Min[I,T,e]*Flows_Gas_Min[I,T,e] <= K2[e]/Length_Pipes[e]*(lambda_Min[I,T,e]))
        @constraint(m, [I = 1:T_inv, T = 1:T_ops, e = 1:EDGES_GAS], Flows_Gas_Max[I,T,e]*Flows_Gas_Max[I,T,e] <= K2[e]/Length_Pipes[e]*(lambda_Max[I,T,e]))
    end
    
    
    ###############################################################################
    ### Gaseous fuel Supply = Demand
    ###############################################################################
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], SUPPLY_GAS_slack_Min[I,T,n] <= MAXSLACK[n])
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], SUPPLY_GAS_slack_Max[I,T,n] <= MAXSLACK[n])
    
#    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Max[I,T,n] >= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
#    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Min[I,T,n] <= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,s]-discharging_GAS[I,T,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Max[I,T,n] >= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_GAS], NetGasDemands_Min[I,T,n] <= -1*sum(STORAGE_GAS_NodalLoc_GAS[n,s]*(charging_GAS[I,T,t,s]-discharging_GAS[I,T,t,s]) for s = 1:STORAGE_GAS) - Demand_GAS[I,T,t,n] - sum(GEN_NodalLoc_GAS[n,g]*(generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) + sum(P2G_NodalLoc_GAS[n,d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G))
    
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas_Min[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack_Min[I,T,n] + NetGasDemands_Min[I,T,n] == 0)
    @constraint(m, [I = 1:T_inv, T = 1:T_ops, n = 1:NODES_GAS], LHV_CH4*(sum(-1*A_GAS[n,e]*Flows_Gas_Max[I,T,e] for e = 1:EDGES_GAS)) + SUPPLY_GAS_slack_Max[I,T,n] + NetGasDemands_Max[I,T,n] == 0)

end


################################################################################
### Policy-driven constraints
################################################################################

### Nominal allocation of net-zero emissions gas consumption to each sector, 
# not to exceed the amount of gaseous energy consumed by that sector
# See Eq. 2.62 in Von Wald thesis
################################################################################
@variable(m, CleanGas_gassector[I = 1:T_inv] >= 0)
@variable(m, CleanGas_powersector[I = 1:T_inv] >= 0)
@constraint(m, [I = 1:T_inv], CleanGas_gassector[I] + CleanGas_powersector[I] == sum(weights[I,T]*8760/t_ops*sum(sum(P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops))
@constraint(m, [I = 1:T_inv], CleanGas_powersector[I] <= sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + startup_GEN[I,T,t,g]*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops))
@constraint(m, [I = 1:T_inv], CleanGas_gassector[I] <= sum(weights[I,T]*8760/t_ops*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops))


### Slack variables for emissions constraint violations that could represent:
#   > use of negative emissions offsets at a fixed cost 
#   > excess carbon emissions evaluated in objective function at a "social cost of carbon"
# The power and gas sector emissions that exceed their respective emissions intensity constraint are constrained by the maxOffsets share of total emissions liabilities
# See Eq. 2.61 in Von Wald thesis
################################################################################
@variable(m, excess_powerEmissions[I = 1:T_inv] >= 0)
@variable(m, excess_gasEmissions[I = 1:T_inv] >= 0)
@constraint(m, [I = 1:T_inv], excess_powerEmissions[I] <= maxOffsets_elec[I]*(sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g])*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)))
@constraint(m, [I = 1:T_inv], excess_gasEmissions[I] <= maxOffsets_gas[I]*(sum(weights[I,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)))

### Constraints on the emissions intensity of the electric and gas sectors.
# Here, some emissions from gaseous fuel consumption are offset by the nominal allocation of net-zero emission gas to each sector.
# See Eq. 2.60 in Von Wald thesis
################################################################################
# @constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g])*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*(CleanGas_powersector[I])  <= EI_ElecSector[I]/1000*sum(weights[I,T]*8760/t_ops*sum(sum(generation[I,T,t,g0] for g0 = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) + excess_powerEmissions[I])
@constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g])*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)  <= EI_ElecSector[I]/1000*sum(weights[I,T]*8760/t_ops*sum(sum(generation[I,T,t,g0] for g0 = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) + excess_powerEmissions[I])
@constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[I] <= EI_GasSector[I]/1000*sum(weights[I,T]*8760/t_ops*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) + excess_gasEmissions[I])


### Maximum biomethane production and use of sustainable bio-energy. 
# Total bio-energy constraint is included in the case where net-zero emissions fuel production units can be used to generate methane or LPG fuel
# but must compete for sustainable biomass feedstocks.
# See Eq. 2.63 in Von Wald thesis
################################################################################
@constraint(m, [I = 1:T_inv], maxBiomethane[I] >= sum(weights[I,T]*8760/t_ops*sum(sum(ISBIOMETHANE[d]*P2G_dispatch[I,T,t,d]*eta_P2G[d] for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops))
@constraint(m, [I = 1:T_inv], maxSustainableBiomass[I] >= sum(weights[I,T]*8760/t_ops*sum(sum(ISBIOMASS[d]*P2G_dispatch[I,T,t,d]*(eta_P2G[d]+eta_P2L[d]) for d = 1:P2G) for t = 1:t_ops) for T = 1:T_ops))


###############################################################################
### Gas distribution retirement constraint set
# See Eq. 2.69 in Von Wald thesis
###############################################################################
# User can select whether to allow the model to decide when the retire the gas system
# by setting gasdistretirement_allowed = 1.
# In this case, binary variables are introduced to indicate for each distribution system
# when that system is shut down.
if gasdistretirement_allowed == 1
    @variable(m, distSysRetirement_GAS[I = 1:T_inv, d = 1:DIST_GAS], Bin)
    @constraint(m, [I = 1:T_inv, d = 1:DIST_GAS], (1-sum(distSysRetirement_GAS[j,d] for j = 1:I))*sum(sum(InitialAppliancePopulation[a]/1000*ApplianceProfilesGAS[t,a] for t = 1:8760) for a = 1:APPLIANCES) >= sum(sum(sum(APP_DistSystemLoc_GAS[d,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES) for t= 1:t_ops) for T = 1:T_ops))
    @constraint(m, [d = 1:DIST_GAS], sum(distSysRetirement_GAS[j,d] for j = 1:T_inv) <= 1)
end
# If no gas distribution retirement is contemplated by the model, the distSysRetirement_GAS indicators are fixed parameters
if gasdistretirement_allowed == 0
    distSysRetirement_GAS = zeros(T_inv,DIST_GAS)
end
# Use may also select whether to force the model to shut down gas distribution systems
# by using gasdistretirement_forced with a 1 in the year where the system must be shut down.
# Currently all distribution systems must be shut down in the same year.
if sum(gasdistretirement_forced) >= 1
    distSysRetirement_GAS = zeros(T_inv,DIST_GAS)
    for ii = 1:T_inv
        distSysRetirement_GAS[ii,:] .= gasdistretirement_forced[ii]
    end
    # And delivered gas volumes are constrained to be 0 during and after this designated shut-down year:
    @constraint(m, [I = 1:T_inv, d = 1:DIST_GAS], (1-sum(distSysRetirement_GAS[j,d] for j = 1:I))*sum(sum(InitialAppliancePopulation[a]/1000*ApplianceProfilesGAS[t,a] for t = 1:8760) for a = 1:APPLIANCES) >= sum(sum(sum(APP_DistSystemLoc_GAS[d,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES) for t= 1:t_ops) for T = 1:T_ops))
end

## The gas distribution system fixed costs are then computed based on these shut-down decisions:
# Not explicitly included in Von Wald thesis
###############################################################################
@variable(m, gasdistsyst_Cost[I = 1:T_inv] >= 0)
# Gas distribution system cost includes several terms that will be either active or zero depending on whether the retirement decision is made:
@constraint(m, [I = 1:T_inv], gasdistsyst_Cost[I] == sum(sum(AccDepGasSyst_FixedCosts[j,I]/1000*distSysRetirement_GAS[j,d] for d = 1:DIST_GAS) for j = 1:T_inv) + sum(BAUGasSyst_FixedCosts[I]/1000*(1-sum(distSysRetirement_GAS[j,d] for j = 1:T_inv)) for d = 1:DIST_GAS))


###############################################################################
### Ancillary customer electrification costs 
# See Eq. 2.66 in Von Wald thesis
###############################################################################
@variable(m, applianceInfrastructureCosts[I = 1:T_inv, a = 1:APPLIANCES] >= 0)
@constraint(m,[I = 1, a = 1:APPLIANCES], applianceInfrastructureCosts[I,a] >= unitsbuilt_APPS[I,a]*1000*upgrade_cost[a])
if T_inv > 1
    @constraint(m,[I = 2:T_inv, a = 1:APPLIANCES], applianceInfrastructureCosts[I,a] >= 1000*(unitsbuilt_APPS[I,a] - sum(round(cumulativefailurefrac[a,v,I]-cumulativefailurefrac[a,v,I-1],digits = 4)*unitsbuilt_APPS[v,a] for v = 1:I-1))*upgrade_cost[a])
end

###############################################################################
# Generalized distribution capital costs associated with peak electrical demand
# See Eq. 2.67/2.68 in Von Wald thesis
# Currently implement a few different approaches to compute:
# (a) the total peak power demand at each node
# (b) the peak distribution-level demand at each node
# (c) the peak incremental distribution-level demand due to appliance electrification (i.e., above baseline demand)
# Current version uses PeakDistDemandInc in objective function, but an argument could be made
# that the peak costs should be evaluated with respect to system-wide coincident peak, as opposed to
# the sum of individual nodal peaks.
###############################################################################
@variable(m, PeakDemand[I = 1:T_inv, n = 1:NODES_ELEC] >= 0)
@variable(m, PeakDistDemand[I = 1:T_inv, n = 1:NODES_ELEC] >= 0)
@variable(m, PeakDistDemandInc[I = 1:T_inv, n = 1:NODES_ELEC] >= 0)
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], PeakDemand[I,n] >= Demand_ELEC[I,T,t,n] + sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(charging_ELEC[I,T,t,s]-discharging_ELEC[I,T,t,s]) for s = 1:STORAGE_ELEC) + sum(P2G_NodalLoc_ELEC[n,d]*P2G_dispatch[I,T,t,d]*(1-ISBIOMETHANE[d]) for d = 1:P2G))
@constraint(m, [I = 1:T_inv, t = 1:8760, n = 1:NODES_ELEC], PeakDistDemand[I,n] >= D_Elec[t,n] + 1000*sum(APPLIANCES_NodalLoc_ELEC[n,a]*(unitsremaining_APPS[I,a])*ApplianceProfilesELEC[t,a] for a = 1:APPLIANCES))
@constraint(m, [I = 1:T_inv, T = 1:T_ops, t = 1:t_ops, n = 1:NODES_ELEC], PeakDistDemandInc[I,n] >= Demand_ELEC[I,T,t,n] - BaselineDemand_ELEC[I,T,t,n])


###############################################################################
### Objective function = total societal costs [$/yr]
# See Eq. 2.1 in Von Wald thesis
###############################################################################
# @objective(m, Min, sum(discountfactor[i]*(sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*CRF_GEN[g]*CAPEX_GEN[i0,g] for i0 = 1:i) + UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*FOM_GEN[i,g] for g = 1:GEN) + sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])/1000*generation[i,T,t,g] for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])/1000*sum(startup_GEN[i,T,t,g] for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s]-unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC) + sum(sum(CRF_APPLIANCES[a]*max(min((Years[i0]+ApplianceLifetime[a])-Years[i],1),0)*(CAPEX_APPLIANCES[i0,a]*unitsbuilt_APPS[i0,a]*1000 + applianceInfrastructureCosts[i0,a])  for i0 =1:i)/1000 for a = 1:APPLIANCES) + sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*max(min((Years[i0]+EconomicLifetime_P2G[d])-Years[i],1),0)*CRF_P2G[d]*CAPEX_P2G[i0,d] for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d]-unitsretired_P2G[i0,d] for i0 = 1:i))*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]/1000*P2G_dispatch[i,T,t,d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops) + sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_GAS[s])-Years[i],1),0)*CRF_STORAGE_GAS[s]*CAPEX_STORAGE_GAS[i0,s] for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s]-unitsretired_STORAGE_GAS[i0,s] for i0 = 1:i))*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS) + Cost_DistributionInfrastructure*sum(PeakDistDemandInc[i,n] for n = 1:NODES_ELEC) + sum(weights[i,T]*8760/t_ops*CommodityCost_NG[i]/1000*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - CommodityCost_NG[i]*(CleanGas_gassector[i] + CleanGas_powersector[i])/1000 + gasdistsyst_Cost[i] + offsets_Cost[i]/1000*(excess_powerEmissions[i] + excess_gasEmissions[i]) + sum(sum(max(min((Years[i0]+EconomicLifetime_ELECTrans)-Years[i],1),0)*CRF_ELECTrans*CAPEX_ELECTrans[e]*addflow_TRANS_ELEC[i0,e] for i0 = 1:i) for e = 1:EDGES_ELEC)) for i = 1:T_inv))

@objective(m, Min, sum(discountfactor[i]*(sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*CRF_GEN[g]*CAPEX_GEN[i0,g] for i0 = 1:i) + UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*FOM_GEN[i,g] for g = 1:GEN) +  sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])/1000*generation[i,T,t,g] for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])/1000*sum(startup_GEN[i,T,t,g] for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s]-unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC) + sum(sum(CRF_APPLIANCES[a]*max(min((Years[i0]+ApplianceLifetime[a])-Years[i],1),0)*(CAPEX_APPLIANCES[i0,a]*unitsbuilt_APPS[i0,a]*1000 + applianceInfrastructureCosts[i0,a])  for i0 =1:i)/1000 for a = 1:APPLIANCES) + sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*max(min((Years[i0]+EconomicLifetime_P2G[d])-Years[i],1),0)*CRF_P2G[d]*CAPEX_P2G[i0,d] for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d]-unitsretired_P2G[i0,d] for i0 = 1:i))*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]/1000*P2G_dispatch[i,T,t,d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops) + sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_GAS[s])-Years[i],1),0)*CRF_STORAGE_GAS[s]*CAPEX_STORAGE_GAS[i0,s] for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s]-unitsretired_STORAGE_GAS[i0,s] for i0 = 1:i))*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS) +  Cost_DistributionInfrastructure*sum(PeakDistDemandInc[i,n] for n = 1:NODES_ELEC) 
            + sum(weights[i,T]*8760/t_ops*sum(costOfGasStorage/1000*sum(charging_GAS[i,T,t,s] for s = 1:STORAGE_GAS) for t = 1:t_ops) for T = 1:T_ops) 
            + sum(weights[i,T]*8760*CommodityCost_NG[i,T]/1000*sum(SUPPLY_GAS_slack[i,T,n] for n = 1:NODES_GAS) for T = 1:T_ops) 
#           - CommodityCost_NG[i,1]*(CleanGas_gassector[i] + CleanGas_powersector[i])/1000
            + gasdistsyst_Cost[i] + offsets_Cost[i]/1000*(excess_powerEmissions[i] + excess_gasEmissions[i])
            + sum(sum(max(min((Years[i0]+EconomicLifetime_ELECTrans)-Years[i],1),0)*CRF_ELECTrans*CAPEX_ELECTrans[e]*addflow_TRANS_ELEC[i0,e] for i0 = 1:i) for e = 1:EDGES_ELEC)) for i = 1:T_inv))


status = optimize!(m)