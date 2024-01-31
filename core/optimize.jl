###############################################################################
### Objective function = total societal costs [$/yr]
# See Eq. 2.1 in Von Wald thesis
###############################################################################

@objective(m, Min, sum(discountfactor[i]*(
    sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[i],1),0)*CRF_GEN[g]*CAPEX_GEN[i0,g] for i0 = 1:i)
    + UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:i))*FOM_GEN[i,g] for g = 1:GEN)
    + sum(weights[i,T]*8760/t_ops*sum(sum((VOM_GEN[i,g]+HeatRate[g]*FuelCosts[i,g])/1000*generation[i,T,t,g] for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops)
    + sum(weights[i,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[i,g])/1000*sum(startup_GEN[i,T,t,g] for t = 1:t_ops)  for g = 1:GEN) for T = 1:T_ops)
    + sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[i],1),0)*CRF_STORAGE_ELEC[s]*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:i)
    + UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s]-unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:i))*FOM_STORAGE_ELEC[i,s] for s = 1:STORAGE_ELEC)
    + sum(sum(CRF_APPLIANCES[a]*max(min((Years[i0]+ApplianceLifetime[a])-Years[i],1),0)*(CAPEX_APPLIANCES[i0,a]*unitsbuilt_APPS[i0,a]*1000 + applianceInfrastructureCosts[i0,a])  for i0 =1:i)/1000 for a = 1:APPLIANCES)
    + sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*max(min((Years[i0]+EconomicLifetime_P2G[d])-Years[i],1),0)*CRF_P2G[d]*CAPEX_P2G[i0,d] for i0 = 1:i)
    + UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d]-unitsretired_P2G[i0,d] for i0 = 1:i))*FOM_P2G[i,d] for d = 1:P2G) + sum(weights[i,T]*8760/t_ops*sum(sum(VOM_P2G[i,d]/1000*P2G_dispatch[i,T,t,d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops)
    + sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_GAS[s])-Years[i],1),0)*CRF_STORAGE_GAS[s]*CAPEX_STORAGE_GAS[i0,s] for i0 = 1:i)
    + UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s]-unitsretired_STORAGE_GAS[i0,s] for i0 = 1:i))*FOM_STORAGE_GAS[i,s] for s = 1:STORAGE_GAS)
    + Cost_DistributionInfrastructure*sum(PeakDistDemandInc[i,n] for n = 1:NODES_ELEC) 
    + sum(AMMORTIZED_ELECTrans[e] for e = 1:EDGES_ELEC)/1000
    + sum(weights[i,T]*8760/t_ops*sum(costOfGasStorage/1000*sum(charging_GAS[i,T,t,s] for s = 1:STORAGE_GAS) for t = 1:t_ops) for T = 1:T_ops) 
    + sum(weights[i,T]*8760*CommodityCost_NG[i,T]/1000*sum(SUPPLY_GAS_slack[i,T,n] for n = 1:NODES_GAS) for T = 1:T_ops)
    + gasdistsyst_Cost[i] + offsets_Cost[i]/1000*(excess_powerEmissions[i] + excess_gasEmissions[i])
    + sum(sum(max(min((Years[i0]+EconomicLifetime_ELECTrans)-Years[i],1),0)*(CRF_ELECTrans+ElecTransmissionOperatingCosts)*CAPEX_ELECTrans[e]*addflow_TRANS_ELEC[i0,e] for i0 = 1:i)/1000 for e = 1:EDGES_ELEC)) for i = 1:T_inv))

status = optimize!(m)