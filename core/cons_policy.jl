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
# @constraint(m, [I = 1:T_inv], excess_powerEmissions[I] <= maxOffsets_elec[I]*(sum(weights[I,T]*8760/t_ops*sum(sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g])*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)))
# @constraint(m, [I = 1:T_inv], excess_gasEmissions[I] <= maxOffsets_gas[I]*(sum(weights[I,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)))
@constraint(m, [I = 1:T_inv], excess_powerEmissions[I] + excess_gasEmissions[I] <= maxOffsets[I]*initialEmissions)

### Constraints on the emissions intensity of the electric and gas sectors.
# Here, some emissions from gaseous fuel consumption are offset by the nominal allocation of net-zero emission gas to each sector.
# See Eq. 2.60 in Von Wald thesis
################################################################################
@variable(m, CleanGas_GEN[I = 1:T_inv, g = 1:GEN] >= 0)     # MWh NG
@constraint(m, [I = 1:T_inv], sum(CleanGas_GEN[I, g] for g = 1:GEN) == CleanGas_powersector[I])
@constraint(m, [I = 1:T_inv, g = 1:GEN], CleanGas_GEN[I, g]/MWh_PER_MMBTU <= sum(weights[I,T]*8760/t_ops*sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g])*NG_fueled[g] for t = 1:t_ops) for T = 1:T_ops))
@constraint(m, [I = 1:T_inv], sum((sum(weights[I,T]*8760/t_ops*sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g]) for t = 1:t_ops) for T = 1:T_ops) - CleanGas_GEN[I, g]/MWh_PER_MMBTU) *emissions_factors[g] for g = 1:GEN) <= EI_ElecSector[I]/1000*sum(weights[I,T]*8760/t_ops*sum(sum(generation[I,T,t,g0] for g0 = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) + excess_powerEmissions[I])
@constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[I] <= EI_GasSector[I]/1000*sum(weights[I,T]*8760/t_ops*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) + excess_gasEmissions[I])
# @variable(m, powerEmissions[I = 1:T_inv] >= 0)          # MMTCO2
# @variable(m, gasEmissions[I = 1:T_inv] >= 0)            # MMTCO2
# @constraint(m, [I = 1:T_inv], powerEmissions[I] + gasEmissions[I] <= TotalEmissions_Allowed[I])     # MMTCO2
# @constraint(m, [I = 1:T_inv], sum((sum(weights[I,T]*8760/t_ops*sum((generation[I,T,t,g]*HeatRate[g] + StartupFuel[g]*startup_GEN[I,T,t,g]) for t = 1:t_ops) for T = 1:T_ops) - CleanGas_GEN[I, g]/MWh_PER_MMBTU) *emissions_factors[g] for g = 1:GEN) <= powerEmissions[I]*1e6 + excess_powerEmissions[I])
# @constraint(m, [I = 1:T_inv], sum(weights[I,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[I,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[I] <= gasEmissions[I]*1e6 + excess_gasEmissions[I])

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
    # @variable(m, distSysRetirement_GAS[I = 1:T_inv, d = 1:DIST_GAS], Bin)
    # @constraint(m, [I = 1:T_inv, d = 1:DIST_GAS], (1-sum(distSysRetirement_GAS[j,d] for j = 1:I))*sum(sum(InitialAppliancePopulation[a]/1000*ApplianceProfilesGAS[t,a] for t = 1:8760) for a = 1:APPLIANCES) >= sum(sum(sum(APP_DistSystemLoc_GAS[d,a]*(unitsremaining_APPS[I,a])*ApplianceProfiles_GAS[T,t,a] for a = 1:APPLIANCES) for t= 1:t_ops) for T = 1:T_ops))
    @variable(m, distSysRetirement_GAS[I = 1:T_inv, d = 1:DIST_GAS] >= 0)
    @constraint(m, [I = 1:T_inv, d = 1:DIST_GAS], (1-sum(distSysRetirement_GAS[j,d] for j = 1:I))*sum(APP_DistSystemLoc_GAS[d,a]*InitialAppliancePopulation[a]/1000 for a = gasApps) >= sum(APP_DistSystemLoc_GAS[d,a]*(unitsremaining_APPS[I,a]) for a = gasApps))
    @constraint(m, [d = 1:DIST_GAS], sum(distSysRetirement_GAS[j,d] for j = 1:T_inv) <= 1)

# Use may also select whether to force the model to shut down gas distribution systems
# by using gasdistretirement_forced with a 1 in the year where the system must be shut down.
# Currently all distribution systems must be shut down in the same year.
elseif sum(gasdistretirement_forced) >= 1
    region_index = 1:DIST_GAS
    distSysRetirement_GAS = zeros(T_inv,DIST_GAS)
    if region_retire == "North"
        region_index = [1,2,3,4,11,12,16]
    elseif region_retire == "South"
        region_index = [5,6,7,8,9,10,13,14,15]
    elseif region_retire == "Cities"
        region_index = [3,4,6,7,8,9,12]
    end
    for ii = 1:T_inv
        distSysRetirement_GAS[ii,region_index] .= gasdistretirement_forced[ii]
    end
    # And delivered gas volumes are constrained to be 0 during and after this designated shut-down year:
    @constraint(m, [I = 1:T_inv, d = 1:DIST_GAS], (1-sum(distSysRetirement_GAS[j,d] for j = 1:I))*sum(APP_DistSystemLoc_GAS[d,a]*InitialAppliancePopulation[a]/1000 for a = gasApps) >= sum(APP_DistSystemLoc_GAS[d,a]*(unitsremaining_APPS[I,a]) for a = gasApps))

# If no gas distribution retirement is contemplated by the model, the distSysRetirement_GAS indicators are fixed parameters
else
    distSysRetirement_GAS = zeros(T_inv,DIST_GAS)
end

## The gas distribution system fixed costs are then computed based on these shut-down decisions:
# Not explicitly included in Von Wald thesis
###############################################################################
@variable(m, gasdistsyst_Cost[I = 1:T_inv] >= 0)
# Gas distribution system cost includes several terms that will be either active or zero depending on whether the retirement decision is made:
@constraint(m, [I = 1:T_inv], gasdistsyst_Cost[I] == sum(sum(AccDepGasSyst_FixedCosts[d,j,I]/1000*distSysRetirement_GAS[j,d] for d = 1:DIST_GAS) for j = 1:T_inv) + sum(BAUGasSyst_FixedCosts[d,I]/1000*(1-sum(distSysRetirement_GAS[j,d] for j = 1:T_inv)) for d = 1:DIST_GAS))


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