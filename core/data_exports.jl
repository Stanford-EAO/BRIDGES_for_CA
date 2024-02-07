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

EmissionsAndCosts = zeros(25,T_inv)
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
    # T&D (Peak Distribution, New Transmission CAPEX + FOM, Existing CAPEX + FOM)
    EmissionsAndCosts[8,i] = Cost_DistributionInfrastructure*1000*sum(JuMP.value.(PeakDistDemand[i,n]) for n = 1:NODES_ELEC)
    EmissionsAndCosts[20,i] = sum(sum(max(min((Years[i0]+EconomicLifetime_ELECTrans)-Years[i],1),0)*(CRF_ELECTrans+ElecTransmissionOperatingCosts)*CAPEX_ELECTrans[e]*addflow_TRANS_ELEC[i0,e] for i0 = 1:i) for e = 1:EDGES_ELEC)
    EmissionsAndCosts[25,i] = sum(AMMORTIZED_ELECTrans[e] for e = 1:EDGES_ELEC)

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
    EmissionsAndCosts[15,i] = (sum(weights[i,T]*8760/t_ops*sum(StartupFuel[g]*emissions_factors[g]*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g])*HeatRate[g]*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - sum(JuMP.value.(CleanGas_GEN[i, g])/MWh_PER_MMBTU*emissions_factors[g] for g = 1:GEN))/sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g]) for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops)
    EmissionsAndCosts[16,i] = (sum(weights[i,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[i])/sum(weights[i,T]*8760/t_ops*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops)

    EmissionsAndCosts[21,i] = sum(weights[i,T]*8760/t_ops*sum(StartupFuel[g]*emissions_factors[g]*sum(JuMP.value.(startup_GEN[i,T,t,g]) for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops) + sum(weights[i,T]*8760/t_ops*sum(sum(JuMP.value.(generation[i,T,t,g])*HeatRate[g]*emissions_factors[g] for g = 1:GEN) for t = 1:t_ops) for T = 1:T_ops) - sum(JuMP.value.(CleanGas_GEN[i, g])/MWh_PER_MMBTU*emissions_factors[g] for g = 1:GEN)
    EmissionsAndCosts[22,i] = sum(weights[i,T]*8760/t_ops*EF_NG*sum(sum(Demand_GAS[i,T,t,n] for n = 1:NODES_GAS) for t = 1:t_ops) for T = 1:T_ops) - EF_NG*CleanGas_gassector[i]
    EmissionsAndCosts[23,i] = JuMP.value.(excess_powerEmissions[i])
    EmissionsAndCosts[24,i] = JuMP.value.(excess_gasEmissions[i])
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
CSV.write("$(top_dir)/PEAKDEMAND.csv",Tables.table(JuMP.value.(PeakDistDemand)), writeheader = true)
CSV.write("$(top_dir)/PEAKDEMANDINC.csv",Tables.table(JuMP.value.(PeakDistDemandInc)), writeheader = true)

if gasdistretirement_allowed == 1
    CSV.write("$(top_dir)/DISTRETIRE.csv",Tables.table(JuMP.value.(distSysRetirement_GAS)), writeheader = true)
end


## then let's find the indeces of the different P2G
idx_PowerToCH4  = PrimeMover_P2G .!= fill("Electrolysis", length(PrimeMover_P2G))
idx_PowerToH2   = PrimeMover_P2G .== fill("Electrolysis", length(PrimeMover_P2G))

# concatenate
capPowerToCH4_array  = vcat(transpose(preCapacity_P2G[idx_PowerToCH4]), invCapacity_P2G[:,idx_PowerToCH4])
capPowerToH2_array   = vcat(transpose(preCapacity_P2G[idx_PowerToH2]), invCapacity_P2G[:,idx_PowerToH2])

# find sum
capPowerToCH4 = sum(capPowerToCH4_array, dims=2)
capPowerToH2  = sum(capPowerToH2_array, dims=2)

# Specify column names as strings
column_names_P2G = ["Power-to-CH4", "Power-to-H2"]

# Convert data to 1-dimensional arrays
capPowerToCH4  = vec(capPowerToCH4)
capPowerToH2   = vec(capPowerToH2)







#### Storage Capacity in MWh
preCapacity_ELEC = (NumUnits_STORAGE_ELEC .* UnitSize_STORAGE_ELEC .* duration_ELEC)
# define array based on inv periods
invCapacity_ELEC = zeros(T_inv, length(preCapacity_ELEC))
# loop
for i=1:T_inv
    #
    netForInvP = UnitSize_STORAGE_ELEC .* duration_ELEC .* sum( JuMP.value.(unitsbuilt_STORAGE_ELEC[invP,:]) - JuMP.value.(unitsretired_STORAGE_ELEC[invP,:]) for invP = 1:i)
    invCapacity_ELEC[i,:] = preCapacity_ELEC[:] + netForInvP[:]
end
# gas
preCapacity_GAS = UnitSize_STORAGE_GAS .* duration_GAS

## then let's find the indeces of the different storage mechanisms
idx_PHS       = PrimeMover_STORAGE_ELEC .== fill("Pumped hydro storage", length(PrimeMover_STORAGE_ELEC))
idx_LiBattery = PrimeMover_STORAGE_ELEC .== fill("Li-ion battery", length(PrimeMover_STORAGE_ELEC))
idx_H2Storage = PrimeMover_STORAGE_ELEC .== fill("Long-duration storage", length(PrimeMover_STORAGE_ELEC))
idx_FeBattery = PrimeMover_STORAGE_ELEC .== fill("Multi-day storage", length(PrimeMover_STORAGE_ELEC))

# concatenate
capPHS_array       = vcat(transpose(preCapacity_ELEC[idx_PHS]), invCapacity_ELEC[:,idx_PHS])
capLiBattery_array = vcat(transpose(preCapacity_ELEC[idx_LiBattery]), invCapacity_ELEC[:,idx_LiBattery])
capH2Storage_array = vcat(transpose(preCapacity_ELEC[idx_H2Storage]), invCapacity_ELEC[:,idx_H2Storage])
capFeBattery_array = vcat(transpose(preCapacity_ELEC[idx_FeBattery]), invCapacity_ELEC[:,idx_FeBattery])

# find sum
capPHS       = sum(capPHS_array, dims=2)
capLiBattery = sum(capLiBattery_array, dims=2)
capH2Storage = sum(capH2Storage_array, dims=2)
capFeBattery = sum(capFeBattery_array, dims=2)

# Specify column names as strings
column_names = ["Li-ion battery", "Pumped hydro storage", "Hydrogen Storage", "Fe-air Battery"]

# Convert data to 1-dimensional arrays
capPHS       = vec(capPHS)
capLiBattery = vec(capLiBattery)
capH2Storage = vec(capH2Storage)
capFeBattery = vec(capFeBattery)

# Create a DataFrame with column names and data
df = DataFrame(column_names[1] => capLiBattery, column_names[2] => capPHS, column_names[3] => capH2Storage, column_names[4] => capFeBattery)

#
resultName   = "$(top_dir)/STORAGE_ENERGY_CAPACITIES"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
outputName = string(resultName ,temporalName, scenarioName,".csv")
#
CSV.write(outputName, df, writeheader = true)




#
#### Storage Capacity in MW
preCapacity_ELEC = (NumUnits_STORAGE_ELEC .* UnitSize_STORAGE_ELEC)
# define array based on inv periods
invCapacity_ELEC = zeros(T_inv, length(preCapacity_ELEC))
# loop
for i=1:T_inv
    #
    netForInvP = UnitSize_STORAGE_ELEC .* sum( JuMP.value.(unitsbuilt_STORAGE_ELEC[invP,:]) - JuMP.value.(unitsretired_STORAGE_ELEC[invP,:]) for invP = 1:i)
    invCapacity_ELEC[i,:] = preCapacity_ELEC[:] + netForInvP[:]
end
# gas
preCapacity_GAS = UnitSize_STORAGE_GAS

## then let's find the indeces of the different storage mechanisms
idx_PHS       = PrimeMover_STORAGE_ELEC .== fill("Pumped hydro storage", length(PrimeMover_STORAGE_ELEC))
idx_LiBattery = PrimeMover_STORAGE_ELEC .== fill("Li-ion battery", length(PrimeMover_STORAGE_ELEC))
idx_H2Storage = PrimeMover_STORAGE_ELEC .== fill("Long-duration storage", length(PrimeMover_STORAGE_ELEC))
idx_FeBattery = PrimeMover_STORAGE_ELEC .== fill("Multi-day storage", length(PrimeMover_STORAGE_ELEC))

# concatenate
capPHS_array       = vcat(transpose(preCapacity_ELEC[idx_PHS]), invCapacity_ELEC[:,idx_PHS])
capLiBattery_array = vcat(transpose(preCapacity_ELEC[idx_LiBattery]), invCapacity_ELEC[:,idx_LiBattery])
capH2Storage_array = vcat(transpose(preCapacity_ELEC[idx_H2Storage]), invCapacity_ELEC[:,idx_H2Storage])
capFeBattery_array = vcat(transpose(preCapacity_ELEC[idx_FeBattery]), invCapacity_ELEC[:,idx_FeBattery])

# find sum
capPHS       = sum(capPHS_array, dims=2)
capLiBattery = sum(capLiBattery_array, dims=2)
capH2Storage = sum(capH2Storage_array, dims=2)
capFeBattery = sum(capFeBattery_array, dims=2)

# Specify column names as strings
column_names = ["Li-ion battery", "Pumped hydro storage", "Hydrogen Storage", "Fe-air Battery"]

# Convert data to 1-dimensional arrays
capPHS       = vec(capPHS)
capLiBattery = vec(capLiBattery)
capH2Storage = vec(capH2Storage)
capFeBattery = vec(capFeBattery)


# Create a DataFrame with column names and data
df = DataFrame(column_names[1] => capLiBattery, column_names[2] => capPHS, column_names[3] => capH2Storage, column_names[4] => capFeBattery,
               column_names_P2G[1] => capPowerToCH4, column_names_P2G[2] => capPowerToH2)

#
resultName   = "$(top_dir)/STORAGE_POWER_CAPACITIES"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
outputName = string(resultName, temporalName, scenarioName,".csv")
#
CSV.write(outputName, df, writeheader = true)



#### Objective Function Value
outputName = string("$(top_dir)/ObjectiveFunctionValue", temporalName, scenarioName,".csv")
CSV.write(outputName, DataFrame(Value = [objective_value(m)]), writeheader = true)



##########################################################################
########################### GAS FLOWS EXPORTS ############################
##########################################################################
# function for hourly output
function process_gasOutput(gas_output, outputName)
    # Determine the column names for each slice in the third dimension
    column_names = String[]  # Initialize as an empty string array
    for i = 1:size(gas_output, 1)
        prefix = "InvPeriod_$(i)"
        for j = 1:size(gas_output, 2)
            suffix = "RepDay_$(j)"
            push!(column_names, string(prefix, "+", suffix))
        end
    end

    # Reshape the 3D matrix into a 2D matrix
    gas_output2D = gas_output[1, :,:]'
    if size(gas_output, 1) > 1
        for i = 2:size(gas_output, 1)
            gas_output2D = hcat(gas_output2D, gas_output[i, :,:]')
        end
    end

    # Create a DataFrame with column names
    df = DataFrame(gas_output2D, column_names)

    # Rename the columns to the specified names
    rename!(df, Symbol.(column_names))

    # Save the DataFrame to a CSV file
    CSV.write(outputName, df)
end

### P2G
P2G_output = JuMP.value.( sum(P2G_dispatch[:,:,:,d]*eta_P2G[d] for d=1:P2G) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/P2G_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(P2G_output, outputName)

### Imports
Imports_output = JuMP.value.( sum(SUPPLY_GAS_slack[:,:,n] for n=1:NODES_GAS) )
Imports_output = cat([Imports_output for _ in 1:24]..., dims=3)
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/Imports_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(Imports_output, outputName)

### Charging
Charging_output = JuMP.value.( sum(charging_GAS[:,:,:,s] for s = 1:STORAGE_GAS) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/Charging_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(Charging_output, outputName)

### Discharging
Discharging_output = JuMP.value.( sum(discharging_GAS[:,:,:,s] for s = 1:STORAGE_GAS) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/Discharging_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(Discharging_output, outputName)

### Demand
Demand_output = JuMP.value.( sum(sum(NominalGasOfftakes[:,:,:,n,g]*LHV[g]*MolarMass[g] for g=1:GAS_COMPONENTS) for n=1:NODES_GAS) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/Demand_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(Demand_output, outputName)





##########################################################################
########################## ELEC FLOWS EXPORTS ############################
##########################################################################
#
### using the same code as gasOutput
#
### P2G
ELEC_P2G_output = JuMP.value.( sum(P2G_dispatch[:,:,:,d]*(1-ISBIOMETHANE[d]) for d = 1:P2G) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/ELEC_P2G_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_P2G_output, outputName)

### Generation
ELEC_Generation_output = JuMP.value.( sum(generation[:,:,:,g] for g = 1:GEN) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/ELEC_Gen_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Generation_output, outputName)

### Curtailment
ELEC_Curtailment_output = JuMP.value.( sum(curtailmentRE[:,:,:,g] for g = 1:GEN) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/ELEC_Curtailment_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Curtailment_output, outputName)

### Charging
ELEC_Charging_output = JuMP.value.( sum(charging_ELEC[:,:,:,s] for s = 1:STORAGE_ELEC) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/ELEC_Charging_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Charging_output, outputName)

### Discharging
ELEC_Discharging_output = JuMP.value.( sum(discharging_ELEC[:,:,:,s] for s = 1:STORAGE_ELEC) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/ELEC_Discharging_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Discharging_output, outputName)

### Demand
ELEC_Demand_output = JuMP.value.( sum(Demand_ELEC[:,:,:,n] for n = 1:NODES_ELEC) )
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/ELEC_Demand_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Demand_output, outputName)





##########################################################################
############### STORAGE CHARGE/DISCHARGE ELEC FLOWS EXPORTS ##############
##########################################################################
# the data should sum up to the charging/discharging from the previous block

# data form before:
### Charging
ELEC_Charging_vec = JuMP.value.( charging_ELEC[:,:,:,:] )
### Discharging
ELEC_Discharging_vec = JuMP.value.( discharging_ELEC[:,:,:,:] )

# for each storage technology:
idx_LiBattery = PrimeMover_STORAGE_ELEC .== fill("Li-ion battery", length(PrimeMover_STORAGE_ELEC))
idx_PHS       = PrimeMover_STORAGE_ELEC .== fill("Pumped hydro storage", length(PrimeMover_STORAGE_ELEC))
idx_FeBattery = PrimeMover_STORAGE_ELEC .== fill("Multi-day storage", length(PrimeMover_STORAGE_ELEC))
idx_H2Storage = PrimeMover_STORAGE_ELEC .== fill("Long-duration storage", length(PrimeMover_STORAGE_ELEC))


### Li-ion
#
ELEC_Charging_LiIon    = sum(ELEC_Charging_vec[:,:,:,idx_LiBattery], dims=4)
ELEC_Discharging_LiIon = sum(ELEC_Discharging_vec[:,:,:,idx_LiBattery], dims=4)
## charging
resultName   = "$(top_dir)/ELEC_Charging_LiIon_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Charging_LiIon, outputName)

## discharging
resultName   = "$(top_dir)/ELEC_Discharging_LiIon_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Discharging_LiIon, outputName)



### PHS
#
ELEC_Charging_PHS    = sum(ELEC_Charging_vec[:,:,:,idx_PHS], dims=4)
ELEC_Discharging_PHS = sum(ELEC_Discharging_vec[:,:,:,idx_PHS], dims=4)
## charging
resultName   = "$(top_dir)/ELEC_Charging_PHS_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Charging_PHS, outputName)

## discharging
resultName   = "$(top_dir)/ELEC_Discharging_PHS_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Discharging_PHS, outputName)



### FeBattery
#
ELEC_Charging_FeBattery    = sum(ELEC_Charging_vec[:,:,:,idx_FeBattery], dims=4)
ELEC_Discharging_FeBattery = sum(ELEC_Discharging_vec[:,:,:,idx_FeBattery], dims=4)
## charging
resultName   = "$(top_dir)/ELEC_Charging_FeBattery_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Charging_FeBattery, outputName)

## discharging
resultName   = "$(top_dir)/ELEC_Discharging_FeBattery_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Discharging_FeBattery, outputName)



### H2Storage
#
ELEC_Charging_H2Storage    = sum(ELEC_Charging_vec[:,:,:,idx_H2Storage], dims=4)
ELEC_Discharging_H2Storage = sum(ELEC_Discharging_vec[:,:,:,idx_H2Storage], dims=4)
## charging
resultName   = "$(top_dir)/ELEC_Charging_H2Storage_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Charging_H2Storage, outputName)

## discharging
resultName   = "$(top_dir)/ELEC_Discharging_H2Storage_HOURLY"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_gasOutput(ELEC_Discharging_H2Storage, outputName)





##########################################################################
############################# GEN EXPORTS ################################
##########################################################################

# determine pre capacity
preCapacity_GEN = (NumUnits_GEN .* UnitSize_GEN)

# define array based on inv periods
invCapacity_GEN = zeros(T_inv, length(preCapacity_GEN))
# loop
for i=1:T_inv
    #
    netForInvP = UnitSize_GEN .* (sum( JuMP.value.(unitsbuilt_GEN[invP,:]) - JuMP.value.(unitsretired_GEN[invP,:]) for invP = 1:i))
    invCapacity_GEN[i,:] = preCapacity_GEN[:] + netForInvP[:]
end

# Create an empty dictionary
idx_GEN = Dict{String, Vector{Bool}}()
#
GenTechnology = unique(PrimeMover_GEN)

# loop to get indeces
for i in 1:length(GenTechnology)
    key = GenTechnology[i]
    value = PrimeMover_GEN .== fill(key, length(PrimeMover_GEN))

    # Add the key-value pair to the dictionary
    idx_GEN[key] = BitVector(value)
end


# # loop to concatenate
cap_GEN_matrix = Dict{String, Matrix}()
energy_GEN = Dict{String, Vector}()
for i in 1:length(GenTechnology)
    key = GenTechnology[i]
    idx = idx_GEN[key]
    cap_GEN_matrix[key] = vcat(transpose(preCapacity_GEN[idx]), invCapacity_GEN[:,idx])
end

# add the technologies' power cap
cap_GEN = Dict{String, Vector}()
for i in 1:length(GenTechnology)
    key = GenTechnology[i]
    cap_GEN[key] = vec(sum(cap_GEN_matrix[key], dims=2))
end


# create dataframe
df = DataFrame(cap_GEN)
#
resultName   = "$(top_dir)/GEN_POWER_CAPACITIES"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
outputName = string(resultName, temporalName, scenarioName,".csv")
#
CSV.write(outputName, df, writeheader = true)


######### now for energy generation by generator prime mover:
EnergyGen = sum(weights[:,T].*8760/t_ops.*sum(JuMP.value.(generation[:,T,t,:]) for t = 1:t_ops) for T = 1:T_ops)
#
energy_GEN = Dict{String, Vector}()
# loop to combine different sources to 1
for i in 1:length(GenTechnology)
    key = GenTechnology[i]
    idx = idx_GEN[key]
    energy_GEN[key] = vec(sum(EnergyGen[:,idx], dims=2))
end

# create dataframe
df = DataFrame(energy_GEN)
#
resultName   = "$(top_dir)/GEN_ENERGY_OUTPUT"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
outputName = string(resultName, temporalName, scenarioName,".csv")
#
CSV.write(outputName, df, writeheader = true)



##########################################################################
######################### YEARLY STORAGE EXPORTS #########################
##########################################################################
#
## then let's find the indeces of the different storage mechanisms
idx_LiBattery = PrimeMover_STORAGE_ELEC .== fill("Li-ion battery", length(PrimeMover_STORAGE_ELEC))
idx_PHS       = PrimeMover_STORAGE_ELEC .== fill("Pumped hydro storage", length(PrimeMover_STORAGE_ELEC))
idx_FeBattery = PrimeMover_STORAGE_ELEC .== fill("Multi-day storage", length(PrimeMover_STORAGE_ELEC))
idx_H2Storage = PrimeMover_STORAGE_ELEC .== fill("Long-duration storage", length(PrimeMover_STORAGE_ELEC))

totalStorageSOC_LiBattery = sum(JuMP.value.(SOCTracked_ELEC[:,:,idx_LiBattery]), dims=3)
totalStorageSOC_PHS       = sum(JuMP.value.(SOCTracked_ELEC[:,:,idx_PHS]), dims=3)
totalStorageSOC_FeBattery = sum(JuMP.value.(SOCTracked_ELEC[:,:,idx_FeBattery]), dims=3)
totalStorageSOC_H2Storage = sum(JuMP.value.(SOCTracked_ELEC[:,:,idx_H2Storage]), dims=3)
totalStorageSOC_GAS = sum(JuMP.value.(SOCTracked_GAS), dims=3)

findmax(totalStorageSOC_LiBattery)

storageTypes = ["Li-ion battery", "Pumped hydro storage", "Multi-day storage", "Long-duration storage", "Natural gas storage"]

# Determine the column names for each slice in the third dimension
column_names = String[]  # Initialize as an empty string array
for i = 1:size(storageTypes, 1)
    prefix = storageTypes[i]
    for j = 1:size(totalStorageSOC_LiBattery, 1)
        suffix = "InvPeriod_$(j)"
        push!(column_names, string(prefix, "+", suffix))
    end
end

# create yearly flows matrix, assign it to LiBattery 
# yearlyStorageFlows = totalStorageSOC_LiBattery[1, :,:]
# Li-ion
# if size(totalStorageSOC_LiBattery, 1) > 1
#     for i = 2:size(totalStorageSOC_LiBattery, 1)
#         yearlyStorageFlows = hcat(yearlyStorageFlows, totalStorageSOC_LiBattery[i, :,:])
#     end
# end
# # PHS
# for i = 1:size(totalStorageSOC_PHS, 1)
#     yearlyStorageFlows = hcat(yearlyStorageFlows, totalStorageSOC_PHS[i, :,:])
# end
# # MDS
# for i = 1:size(totalStorageSOC_FeBattery, 1)
#     yearlyStorageFlows = hcat(yearlyStorageFlows, totalStorageSOC_FeBattery[i, :,:])
# end
# # H2
# for i = 1:size(totalStorageSOC_H2Storage, 1)
#     yearlyStorageFlows = hcat(yearlyStorageFlows, totalStorageSOC_H2Storage[i, :,:])
# end
# # GAS
# for i = 1:size(totalStorageSOC_GAS, 1)
#     yearlyStorageFlows = hcat(yearlyStorageFlows, totalStorageSOC_GAS[i, :,:])
# end

function concatenateYearlyOutput(input, output)
    # add
    for i = 1:size(input, 1)
        output = hcat(output, input[i, :,:])
    end
    return output
end

yearlyStorageFlows = totalStorageSOC_LiBattery[1, :,:]
#
yearlyStorageFlows = concatenateYearlyOutput(totalStorageSOC_LiBattery, yearlyStorageFlows)
yearlyStorageFlows = concatenateYearlyOutput(totalStorageSOC_PHS, yearlyStorageFlows)
yearlyStorageFlows = concatenateYearlyOutput(totalStorageSOC_FeBattery, yearlyStorageFlows)
yearlyStorageFlows = concatenateYearlyOutput(totalStorageSOC_H2Storage, yearlyStorageFlows)
yearlyStorageFlows = concatenateYearlyOutput(totalStorageSOC_GAS, yearlyStorageFlows)
yearlyStorageFlows = yearlyStorageFlows[:,2:end]

# Create a DataFrame with column names
df = DataFrame(yearlyStorageFlows, column_names)

# Rename the columns to the specified names
rename!(df, Symbol.(column_names))

#
resultName   = "$(top_dir)/YEARLY_STORAGE_FLOWS"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
outputName = string(resultName, temporalName, scenarioName,".csv")
# Save the DataFrame to a CSV file
CSV.write(outputName, df)



##########################################################################
######################### STORAGE UNITS BY NODE ##########################
##########################################################################
#
### Storage units by node and type
#
#### Storage Units
preCapacity_ELEC = (NumUnits_STORAGE_ELEC)
# define array based on inv periods
invCapacity_Built_ELEC = zeros(T_inv, length(preCapacity_ELEC))
invCapacity_Ret_ELEC = zeros(T_inv, length(preCapacity_ELEC))
# loop
for i = 1:T_inv
    #
    invCapacity_Built_ELEC[i,:] = JuMP.value.(unitsbuilt_STORAGE_ELEC[i,:])
    invCapacity_Ret_ELEC[i,:] = JuMP.value.(unitsretired_STORAGE_ELEC[i,:])
end

# find types 
nodeNumber = unique(ElectricalStorage[:,1])
storageType = unique(PrimeMover_STORAGE_ELEC)
# find number of storage types and nodes
numNodes = length(nodeNumber)
numStorageTypes = length(storageType)
#
invCapacity_Built_CZ_ELEC = zeros(T_inv, numStorageTypes, numNodes)
invCapacity_Ret_CZ_ELEC   = zeros(T_inv, numStorageTypes, numNodes)
#
for i = 1:T_inv
    # generate local version of built/ret for inv period I
    built_InvPeriod = invCapacity_Built_ELEC[i,:]
    ret_InvPeriod   = invCapacity_Ret_ELEC[i,:]
    #
    for j = 1:numStorageTypes
        # find indexing based on storage type
        idx_Storage = PrimeMover_STORAGE_ELEC .== fill(storageType[j], length(PrimeMover_STORAGE_ELEC))
        for k = 1:numNodes
            # find indexing based on node
            idx_Node = ElectricalStorage[:,1] .== fill(nodeNumber[k], length(ElectricalStorage[:,1]))
            # intersection
            idx_intersection = idx_Storage .& idx_Node
            # find intersection and then add to 
            resultsBuilt_intersection = built_InvPeriod[idx_intersection]
            resultsRet_intersection   = ret_InvPeriod[idx_intersection]
            #
            resultsBuilt_array = zeros(length(idx_intersection))
            resultsRet_array   = zeros(length(idx_intersection))
            resultsBuilt_array[findall(idx_intersection .== 1)] = resultsBuilt_intersection
            resultsRet_array[findall(idx_intersection .== 1)]   = resultsRet_intersection
            #
            # add to matrix
            invCapacity_Built_CZ_ELEC[i,j,k] = sum(resultsBuilt_array)
            invCapacity_Ret_CZ_ELEC[i,j,k]   = sum(resultsRet_array)
        end
    end
end



function process_storageUnitsOutput(gas_output, outputName)
    # Determine the column names for each slice in the third dimension
    column_names = String[]  # Initialize as an empty string array
    for i = 1:size(gas_output, 1)
        prefix = "InvPeriod_$(i)"
        for j = 1:size(gas_output, 2)
            suffix = storageType[j]
            push!(column_names, string(prefix, "+", suffix))
        end
    end

    # Reshape the 3D matrix into a 2D matrix
    gas_output2D = gas_output[1, :,:]'
    if size(gas_output, 1) > 1
        for i = 2:size(gas_output, 1)
            gas_output2D = hcat(gas_output2D, gas_output[i, :,:]')
        end
    end

    # Create a DataFrame with column names
    df = DataFrame(gas_output2D, column_names)

    # Rename the columns to the specified names
    rename!(df, Symbol.(column_names))

    # Save the DataFrame to a CSV file
    CSV.write(outputName, df)
end


### Built
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/NumUnitsBuilt_StorageELEC_PerNode"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_storageUnitsOutput(invCapacity_Built_CZ_ELEC, outputName)

### Built
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/NumUnitsRet_StorageELEC_PerNode"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_storageUnitsOutput(invCapacity_Ret_CZ_ELEC, outputName)




##########################################################################
# 
### Storage units by node and type in final year **ONLY**
##
#### Storage Capacity in MW
preCapacity_ELEC = (NumUnits_STORAGE_ELEC .* UnitSize_STORAGE_ELEC)
# define array based on inv periods
invCapacity_ELEC = zeros(1, length(preCapacity_ELEC))
# loop
for i=T_inv
    #
    netForInvP = UnitSize_STORAGE_ELEC .* sum( JuMP.value.(unitsbuilt_STORAGE_ELEC[invP,:]) - JuMP.value.(unitsretired_STORAGE_ELEC[invP,:]) for invP = 1:i)
    invCapacity_ELEC[1,:] = preCapacity_ELEC[:] + netForInvP[:]
end


# find types 
nodeNumber = unique(ElectricalStorage[:,1])
storageType = unique(PrimeMover_STORAGE_ELEC)
# find number of storage types and nodes
numNodes = length(nodeNumber)
numStorageTypes = length(storageType)
#
invCapacity_Total_CZ_ELEC = zeros(1, numStorageTypes, numNodes)
#
for i = T_inv
    # # generate local version of built/ret for inv period I
    net_InvPeriod = invCapacity_ELEC[1,:]
    #
    for j = 1:numStorageTypes
        # find indexing based on storage type
        idx_Storage = PrimeMover_STORAGE_ELEC .== fill(storageType[j], length(PrimeMover_STORAGE_ELEC))
        for k = 1:numNodes
            # find indexing based on node
            idx_Node = ElectricalStorage[:,1] .== fill(nodeNumber[k], length(ElectricalStorage[:,1]))
            # intersection
            idx_intersection = idx_Storage .& idx_Node
            # find intersection and then add to 
            resultsTotal_intersection = net_InvPeriod[idx_intersection]
            #
            resultsTotal_array = zeros(length(idx_intersection))
            resultsTotal_array[findall(idx_intersection .== 1)] = resultsTotal_intersection
            #
            # add to matrix
            invCapacity_Total_CZ_ELEC[1,j,k] = sum(resultsTotal_array)
        end
    end
end



### NetCapacity
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/NetCapacity_StorageELEC_PerNode"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_storageUnitsOutput(invCapacity_Total_CZ_ELEC, outputName)





##########################################################################
########################### GEN UNITS BY NODE ############################
##########################################################################
#
### Storage capacities by node and type
#
#### Storage Units
preCapacity_GEN = (NumUnits_GEN)
# define array based on inv periods
invCapacity_Built_GEN = zeros(T_inv, length(preCapacity_GEN))
invCapacity_Ret_GEN = zeros(T_inv, length(preCapacity_GEN))
# loop
for i = 1:T_inv
    #
    invCapacity_Built_GEN[i,:] = JuMP.value.(unitsbuilt_GEN[i,:])
    invCapacity_Ret_GEN[i,:] = JuMP.value.(unitsretired_GEN[i,:])
end

# find types 
nodeNumber = unique(Generators[:,1])
genType = unique(PrimeMover_GEN)
# find number of storage types and nodes
numNodes = length(nodeNumber)
numGenTypes = length(genType)
#
invCapacity_Built_CZ_GEN = zeros(T_inv, numGenTypes, numNodes)
invCapacity_Ret_CZ_GEN   = zeros(T_inv, numGenTypes, numNodes)
#
for i = 1:T_inv
    # generate local version of built/ret for inv period I
    built_InvPeriod = invCapacity_Built_GEN[i,:]
    ret_InvPeriod   = invCapacity_Ret_GEN[i,:]
    #
    for j = 1:numGenTypes
        # find indexing based on storage type
        idx_Gen = PrimeMover_GEN .== fill(genType[j], length(PrimeMover_GEN))
        for k = 1:numNodes
            # find indexing based on node
            idx_Node = Generators[:,1] .== fill(nodeNumber[k], length(Generators[:,1]))
            # intersection
            idx_intersection = idx_Gen .& idx_Node
            # find intersection and then add to 
            resultsBuilt_intersection = built_InvPeriod[idx_intersection]
            resultsRet_intersection   = ret_InvPeriod[idx_intersection]
            #
            resultsBuilt_array = zeros(length(idx_intersection))
            resultsRet_array   = zeros(length(idx_intersection))
            resultsBuilt_array[findall(idx_intersection .== 1)] = resultsBuilt_intersection
            resultsRet_array[findall(idx_intersection .== 1)]   = resultsRet_intersection
            #
            # add to matrix
            invCapacity_Built_CZ_GEN[i,j,k] = sum(resultsBuilt_array)
            invCapacity_Ret_CZ_GEN[i,j,k]   = sum(resultsRet_array)
        end
    end
end



function process_genUnitsOutput(output, outputName)
    # Determine the column names for each slice in the third dimension
    column_names = String[]  # Initialize as an empty string array
    for i = 1:size(output, 1)
        prefix = "InvPeriod_$(i)"
        for j = 1:size(output, 2)
            suffix = genType[j]
            push!(column_names, string(prefix, "+", suffix))
        end
    end

    # Reshape the 3D matrix into a 2D matrix
    output2D = output[1, :,:]'
    if size(output, 1) > 1
        for i = 2:size(output, 1)
            output2D = hcat(output2D, output[i, :,:]')
        end
    end

    # Create a DataFrame with column names
    df = DataFrame(output2D, column_names)

    # Rename the columns to the specified names
    rename!(df, Symbol.(column_names))

    # Save the DataFrame to a CSV file
    CSV.write(outputName, df)
end


### Built
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/NumUnitsBuilt_GEN_PerNode"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_genUnitsOutput(invCapacity_Built_CZ_GEN, outputName)

### Built
# Specify the path to the CSV file where you want to save the data
resultName   = "$(top_dir)/NumUnitsRet_GEN_PerNode"
temporalName = string("_$(T_inv)","Inv","+","$(N_Periods)","RepDays")
scenarioName = string("_MDS","$(FormEnergy_allowed)","+","PHS","$(PHS_allowed)")
#
outputName = string(resultName, temporalName, scenarioName,".csv")
##
process_genUnitsOutput(invCapacity_Ret_CZ_GEN, outputName)

