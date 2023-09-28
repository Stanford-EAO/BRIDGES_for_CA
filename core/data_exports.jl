################################################################################
################################################################################
## Export results for visualization
################################################################################
################################################################################
unitsbuilt_GEN = JuMP.value.(unitsbuilt_GEN)
unitsretired_GEN = JuMP.value.(unitsretired_GEN)
unitsbuilt_STORAGE_ELEC = JuMP.value.(unitsbuilt_STORAGE_ELEC)
unitsretired_STORAGE_ELEC = JuMP.value.(unitsretired_STORAGE_ELEC)
unitsbuilt_P2G = JuMP.value.(unitsbuilt_P2G)
unitsretired_P2G = JuMP.value.(unitsretired_P2G)
unitsbuilt_STORAGE_GAS = JuMP.value.(unitsbuilt_STORAGE_GAS)
unitsretired_STORAGE_GAS = JuMP.value.(unitsretired_STORAGE_GAS)
unitsbuilt_TRANS_GAS = JuMP.value.(unitsbuilt_TRANS_GAS)
unitsretired_TRANS_GAS = JuMP.value.(unitsretired_TRANS_GAS)
unitsbuilt_TRANS_ELEC = JuMP.value.(unitsbuilt_TRANS_ELEC)
unitsretired_TRANS_ELEC = JuMP.value.(unitsretired_TRANS_ELEC)
addflow_TRANS_ELEC = JuMP.value.(addflow_TRANS_ELEC)


Demand_GAS = JuMP.value.(Demand_GAS)
Demand_ELEC = JuMP.value.(Demand_ELEC)

if gasdistretirement_allowed == 1
    distSysRetirement_GAS = JuMP.value.(distSysRetirement_GAS)
end

CleanGas_gassector = JuMP.value.(CleanGas_gassector)
CleanGas_powersector = JuMP.value.(CleanGas_powersector)
unitsbuilt_APPS= JuMP.value.(unitsbuilt_APPS)

EmissionsAndCosts = zeros(20,T_inv)
for i = 1:T_inv
    # Terms for computing average electricity and gas rates
    xA = sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*CRF_GEN[g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*1000*CAPEX_GEN[i0,g] for i0 = 1:i) + UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*1000*FOM_GEN[i,g] for g = 1:GEN) +  sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])*JuMP.value.(generation[i,T,t,g]) for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*1000*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*1000*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC) + Cost_DistributionInfrastructure*1000*sum(JuMP.value.(PeakDistDemand[i,n]) for n = 1:NODES_ELEC) + offsets_Cost[i]*JuMP.value.(excess_powerEmissions[i]) - CommodityCost_NG[i]*CleanGas_powersector[i]
    xB = CleanGas_powersector[i]
    xC = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
    xD = sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*CRF_P2G[d]*1000*CAPEX_P2G[i0,d]  for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d] - unitsretired_P2G[i0,d] for i0 = 1:i))*1000*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]*JuMP.value.(P2G_dispatch[i,T,t,d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
    xE = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(P2G_dispatch[i,T,t,d])*(1-ISBIOMETHANE[d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
    xF = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(P2G_dispatch[i,T,t,d])*eta_P2G[d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
    xG = CommodityCost_NG[i]*(sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(Demand_GAS[i,T,t,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)) + JuMP.value.(gasdistsyst_Cost[i])*1000 + sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*CRF_STORAGE_GAS[s]*1000*CAPEX_STORAGE_GAS[i0,s]  for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:i))*1000*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS) + offsets_Cost[i]*JuMP.value.(excess_gasEmissions[i])  - CommodityCost_NG[i]*CleanGas_gassector[i]
    xH = CleanGas_gassector[i]
    xI = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(Demand_GAS[i,T,t,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
    xJ = sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])*JuMP.value.(generation[i,T,t,g]) for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + offsets_Cost[i]*JuMP.value.(excess_powerEmissions[i]) - CommodityCost_NG[i]*CleanGas_powersector[i]
    
    # Average electricity rate
    EmissionsAndCosts[1,i] = (xA*xF+xB*xD)/(xF*xC-xB*xE)
    # Average cost of zero-emission gas
    EmissionsAndCosts[3,i] = (xD+xE*EmissionsAndCosts[1,i])/xF
    # Average gas rate
    EmissionsAndCosts[2,i] = (xG+xH*EmissionsAndCosts[3,i])/xI

    # Average cost of zero-emission gas (exposed to marginal cost of electricity)
    EmissionsAndCosts[19,i] =  (xC*xD+xE*xJ)/(xF*xC-xB*xE)
    # Average electricity rate (assessed w.r.t remaining electricity (not used for P2G) and after the revenues provided by P2G)
    EmissionsAndCosts[17,i] = (xA + xB*EmissionsAndCosts[19,i] - xE*(xJ+xB*EmissionsAndCosts[19,i])/xC)/(xC-xE)
    # Average gas rate
    EmissionsAndCosts[18,i] = (xG+xH*EmissionsAndCosts[19,i])/xI
        
    # Gen Capex
    EmissionsAndCosts[4,i] = sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*CRF_GEN[g]*1000*CAPEX_GEN[i0,g] for i0 = 1:i) for g = 1:GEN)
    # Gen FOM
    EmissionsAndCosts[5,i] = sum(UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*1000*FOM_GEN[i,g] for g = 1:GEN)
    # Gen VOM and fuel
    EmissionsAndCosts[6,i] = sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])*JuMP.value.(generation[i,T,t,g]) for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) - (CommodityCost_NG[i])*CleanGas_powersector[i]
    # Storage ELEC
    EmissionsAndCosts[7,i] = sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*1000*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)  + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s] - unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*1000*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC)
    # T&D
    EmissionsAndCosts[8,i] = Cost_DistributionInfrastructure*1000*sum(JuMP.value.(PeakDistDemand[i,n]) for n = 1:NODES_ELEC)
    EmissionsAndCosts[20,i] = sum(sum(max(min((Years[i0]+EconomicLifetime_ELECTrans)-Years[i],1),0)*CRF_ELECTrans*CAPEX_ELECTrans[e]*addflow_TRANS_ELEC[i0,e] for i0 = 1:i) for e = 1:EDGES_ELEC)

    # Gas sector costs
    # Commodity
    EmissionsAndCosts[9,i] = CommodityCost_NG[i]*(sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(Demand_GAS[i,T,t,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)) - CommodityCost_NG[i]*CleanGas_gassector[i]
    # Distribution systems
    EmissionsAndCosts[10,i] = JuMP.value.(gasdistsyst_Cost[i])*1000
    # Storage GAS
    EmissionsAndCosts[11,i] = sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_GAS[s])-Years[i],1),0)*CRF_STORAGE_GAS[s]*1000*CAPEX_STORAGE_GAS[i0,s]  for i0 = 1:i) + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s] - unitsretired_STORAGE_GAS[i0,s]  for i0 = 1:i))*1000*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS)

    # Appliances
    EmissionsAndCosts[12,i] = sum(sum(CRF_APPLIANCES[a]*max(min((Years[i0]+ApplianceLifetime[a])-Years[i],1),0)*(CAPEX_APPLIANCES[i0,a]*unitsbuilt_APPS[i0,a]*1000 + JuMP.value.(applianceInfrastructureCosts[i0,a])) for i0 =1:i) for a = 1:APPLIANCES)

    # P2G costs (without electricity costs, these are included in the electricity generation sectoral costs)
    EmissionsAndCosts[13,i] = sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*max(min((Years[i0]+EconomicLifetime_P2G[d])-Years[i],1),0)*CRF_P2G[d]*1000*CAPEX_P2G[i0,d] for i0 = 1:i) + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d] - unitsretired_P2G[i0,d] for i0 = 1:i))*1000*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum((VOM_P2G[i,d])*JuMP.value.(P2G_dispatch[i,T,t,d]) for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)

    # Negative emissions offsets
    EmissionsAndCosts[14,i] = offsets_Cost[i]*(JuMP.value.(excess_powerEmissions[i]) + JuMP.value.(excess_gasEmissions[i]))
    
    # Emissions intensity of electricity generated and gas delivered (check to make sure constraint is satisfied)
    EmissionsAndCosts[15,i] = (sum(weights[i,T]*8760/t_ops*sum(StartupFuel[g]*emissions_factors[g]*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g])*HeatRate[g]*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_powersector[i])/sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
    EmissionsAndCosts[16,i] = (sum(weights[i,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[i])/sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
end

CapacityBuilt = zeros(T_inv,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1)
for i = 1:T_inv
    CapacityBuilt[i,1:GEN] = unitsbuilt_GEN[i,:].*UnitSize_GEN
    CapacityBuilt[i,GEN+1:GEN+STORAGE_ELEC] = unitsbuilt_STORAGE_ELEC[i,:].*UnitSize_STORAGE_ELEC
    CapacityBuilt[i,GEN+STORAGE_ELEC+1:GEN+STORAGE_ELEC+P2G] = unitsbuilt_P2G[i,:].*UnitSize_P2G
    CapacityBuilt[i,GEN+STORAGE_ELEC+P2G+1:GEN+STORAGE_ELEC+P2G+STORAGE_GAS] = unitsbuilt_STORAGE_GAS[i,:].*UnitSize_STORAGE_GAS
    CapacityBuilt[i,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1] = sum(distSysRetirement_GAS[i,d] for d = 1:DIST_GAS)
end

CapacityRetired = zeros(T_inv,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1)
for i = 1:T_inv
    CapacityRetired[i,1:GEN] = unitsretired_GEN[i,:].*UnitSize_GEN
    CapacityRetired[i,GEN+1:GEN+STORAGE_ELEC] = unitsretired_STORAGE_ELEC[i,:].*UnitSize_STORAGE_ELEC
    CapacityRetired[i,GEN+STORAGE_ELEC+1:GEN+STORAGE_ELEC+P2G] = unitsretired_P2G[i,:].*UnitSize_P2G
    CapacityRetired[i,GEN+STORAGE_ELEC+P2G+1:GEN+STORAGE_ELEC+P2G+STORAGE_GAS] = unitsretired_STORAGE_GAS[i,:].*UnitSize_STORAGE_GAS
    CapacityRetired[i,GEN+STORAGE_ELEC+P2G+STORAGE_GAS+1] = sum(distSysRetirement_GAS[i,d] for d = 1:DIST_GAS)
end

GenerationSave = zeros(T_inv,GEN+P2G+9)
for i = 1:T_inv
     GenerationSave[i,1:GEN] = sum(weights[i,T]*8760/t_ops*sum(JuMP.value.(generation[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+1:GEN+P2G] = sum(weights[i,T]*8760/t_ops*sum(JuMP.value.(P2G_dispatch[i,T,t,:]).*eta_P2G for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+1] = sum(weights[i,T]*8760/t_ops*sum(sum(InitialAppliancePopulation[:].*ApplianceProfiles_GAS[T,t,:]) + sum(BaselineDemand_GAS[i,T,t,:])  for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+2] = sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+3] = CleanGas_powersector[i]
     GenerationSave[i,GEN+P2G+4] = CleanGas_gassector[i]
     GenerationSave[i,GEN+P2G+5] = sum(weights[i,T]*8760/t_ops*sum(sum(InitialAppliancePopulation[:].*ApplianceProfiles_ELEC[T,t,:]) + sum(BaselineDemand_ELEC[i,T,t,:])  for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+6] = sum(weights[i,T]*8760/t_ops*sum(sum(Demand_ELEC[i,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+7] = sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(SUPPLY_GAS_slack[i,T,n]) for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+8] = sum(weights[i,T]*8760/t_ops*sum(sum(sum(GEN_NodalLoc_ELEC[n,g]*JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) + sum(-1*A_ELEC[n,e]*JuMP.value.(Flows_Elec[i,T,t,e]) for e = 1:EDGES_ELEC) - sum(STORAGE_ELEC_NodalLoc_ELEC[n,s]*(JuMP.value.(charging_ELEC[i,T,t,s])-JuMP.value.(discharging_ELEC[i,T,t,s])) for s = 1:STORAGE_ELEC) - Demand_ELEC[i,T,t,n] - sum(P2G_NodalLoc_ELEC[n,d]*JuMP.value.(P2G_dispatch[i,T,t,d])*(1-ISBIOMETHANE[d]) for d = 1:P2G) for n = 1:NODES_ELEC)  for t = 1:t_ops) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(curtailmentRE[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
     GenerationSave[i,GEN+P2G+9] = sum(weights[i,T]*8760/t_ops*sum(sum((JuMP.value.(generation[i,T,t,g])*HeatRate[g] + JuMP.value.(startup_GEN[i,T,t,g])*StartupFuel[g])*MWh_PER_MMBTU*NG_fueled[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
end

HourlyGenFullSave = zeros(T_inv*8760,GEN+STORAGE_ELEC+P2G+1)
HourlyLoadFullSave = zeros(T_inv*8760,4)
HourlyStoredElecFullSave = zeros(T_inv*8760,STORAGE_ELEC)
HourlyTransmissionFullSave = zeros(T_inv*8760,EDGES_ELEC)
DailyGasSOCFullSave = zeros(T_inv*365,STORAGE_GAS+2)
HourlyStoredGasEnergyFullSave = zeros(T_inv*8760,STORAGE_GAS)
DailyGasTransmissionFullSave = zeros(T_inv*Periods_Per_Year,EDGES_GAS)
for i = 1:T_inv
    for c = 1:Periods_Per_Year
        j = Int(RepDays[i,c])
        count = Int((i-1)*8760+(c-1)*t_ops)+1
        HourlyGenFullSave[count:count+t_ops-1,1:GEN] = JuMP.value.(generation[i,j,:,:])
        HourlyGenFullSave[count:count+t_ops-1, GEN+1:GEN+STORAGE_ELEC] = (JuMP.value.(charging_ELEC[i,j,:,:])-JuMP.value.(discharging_ELEC[i,j,:,:]))
        HourlyGenFullSave[count:count+t_ops-1, GEN+STORAGE_ELEC+1:GEN+STORAGE_ELEC+P2G] = JuMP.value.(P2G_dispatch[i,j,:,:].*transpose(ones(P2G)-ISBIOMETHANE))
        HourlyGenFullSave[count:count+t_ops-1, GEN+STORAGE_ELEC+P2G+1] = sum(JuMP.value.(curtailmentRE[i,j,:,:]), dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 1] = sum(Demand_ELEC[i,j,:,:], dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 2] = sum(Demand_GAS[i,j,:,:], dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 3] = sum(BaselineDemand_ELEC[i,j,:,:], dims = 2)
        HourlyLoadFullSave[count:count+t_ops-1, 4] = sum(BaselineDemand_GAS[i,j,:,:], dims = 2)
        HourlyStoredElecFullSave[count:count+t_ops-1,:] = JuMP.value.(storedEnergy_ELEC[i,j,1:t_ops,:])
        HourlyStoredGasEnergyFullSave[count:count+t_ops-1,:] = JuMP.value.(storedEnergy_GAS[i,j,1:t_ops,:])
        if LINKED_PERIODS_STORAGE == 1
            DailyGasSOCFullSave[Int((i-1)*Periods_Per_Year + c),1:STORAGE_GAS] = JuMP.value.(SOCTracked_GAS[i,c,:])
        end
        DailyGasSOCFullSave[Int((i-1)*Periods_Per_Year + c), STORAGE_GAS+1] = JuMP.value.(SUPPLY_GAS_slack[i,j,18])
        DailyGasSOCFullSave[Int((i-1)*Periods_Per_Year + c), STORAGE_GAS+2] = JuMP.value.(SUPPLY_GAS_slack[i,j,20])
        HourlyTransmissionFullSave[count:count+t_ops-1,:] = JuMP.value.(Flows_Elec[i,j,:,:])
        DailyGasTransmissionFullSave[Int((i-1)*Periods_Per_Year + c),:] = JuMP.value.(Flows_Gas[i,j,:])
    end
end

ApplianceDecisions = zeros(3*T_inv,APPLIANCES)
for i = 1:T_inv
    count = (3*i-2)
    ApplianceDecisions[count,:] = unitsbuilt_APPS[i,:]
    ApplianceDecisions[count+1,:] = JuMP.value.(unitsretired_APPS[i,:])
    ApplianceDecisions[count+2,:] = JuMP.value.(unitsremaining_APPS[i,:])
end

# Creates output folder autmatically

function mk_output_dir()
    timestamp = Dates.format(now(), "YYYYmmdd-HHMMSS")
    dir_name = joinpath(@__DIR__, "Output", "$timestamp")
    @assert !ispath(dir_name) "File name already taken"
    mkpath(dir_name)
    return dir_name
end
top_dir = mk_output_dir()
println("Saving to: ",last(top_dir, 15))

CSV.write("$(top_dir)/APPLIANCE_DECISIONS.csv",Tables.table(ApplianceDecisions'), writeheader = true)
CSV.write("$(top_dir)/EMISSIONS_COSTS.csv",Tables.table(EmissionsAndCosts'), writeheader = true)
CSV.write("$(top_dir)/CAPACITY_BUILT.csv",Tables.table(CapacityBuilt'), writeheader = true)
CSV.write("$(top_dir)/CAPACITY_RETIRED.csv",Tables.table(CapacityRetired'), writeheader = true)
CSV.write("$(top_dir)/GENERATION.csv",Tables.table(GenerationSave'), writeheader = true)
CSV.write("$(top_dir)/HOURLY_GENERATION.csv",Tables.table(HourlyGenFullSave'), writeheader = true)
CSV.write("$(top_dir)/HOURLY_LOAD.csv",Tables.table(HourlyLoadFullSave'), writeheader = true)
CSV.write("$(top_dir)/HOURLY_ELECTRIC_TRANSMISSION.csv",Tables.table(HourlyTransmissionFullSave'), writeheader = true)

CSV.write("$(top_dir)/HOURLY_GAS_STOREDENERGY.csv",Tables.table(HourlyStoredGasEnergyFullSave'), writeheader = true)
CSV.write("$(top_dir)/DAILY_GAS_SOC.csv",Tables.table(DailyGasSOCFullSave'), writeheader = true)
CSV.write("$(top_dir)/HOURLY_ELEC_STOREDENERGY.csv",Tables.table(HourlyStoredElecFullSave'), writeheader = true)
CSV.write("$(top_dir)/DAILY_GAS_TRANSMISSION.csv",Tables.table(DailyGasTransmissionFullSave'), writeheader = true)

CSV.write("$(top_dir)/REPDAYS.csv",Tables.table(RepDays[1,:]), writeheader = true)
CSV.write("$(top_dir)/ADDELECFLOW.csv",Tables.table(addflow_TRANS_ELEC), writeheader = true)

