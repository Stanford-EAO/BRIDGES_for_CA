###############################################################################
### Objective function = total societal costs [$000s/yr]
# See Eq. 2.1 in Von Wald thesis
###############################################################################

### COSTS
# transmission costs, operation
@variable(m, costs_transmission[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_transmission[I] == sum(sum(max(min((Years[i0]+EconomicLifetime_ELECTrans)-Years[I],1),0)*(CRF_ELECTrans+ElecTransmissionOperatingCosts)*CAPEX_ELECTrans[e]*addflow_TRANS_ELEC[i0,e] for i0 = 1:I)/1000 for e = 1:EDGES_ELEC) + sum(AMMORTIZED_ELECTrans[e] for e = 1:EDGES_ELEC)/1000 )
# gas storage costs
@variable(m, costs_gasStorage[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_gasStorage[I] == sum(weights[I,T]*8760/t_ops*sum(costOfGasStorage/1000*sum(charging_GAS[I,T,t,s] for s = 1:STORAGE_GAS) for t = 1:t_ops) for T = 1:T_ops) )
# CO2 offsets costs
@variable(m, costs_CO2offsets[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_CO2offsets[I] == offsets_Cost[I]/1000*(excess_powerEmissions[I] + excess_gasEmissions[I] + excess_refEmissions[I]) )
# imported natgas costs
@variable(m, costs_NGimports[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_NGimports[I] == sum(weights[I,T]*8760*CommodityCost_NG[I,T]/1000*sum(SUPPLY_GAS_slack[I,T,n] for n = 1:NODES_GAS) for T = 1:T_ops) )
# generator capital costs, ammortized to each investment period
@variable(m, costs_GENcapital[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_GENcapital[I] == sum(UnitSize_GEN[g]*sum(unitsbuilt_GEN[i0,g]*max(min((Years[i0]+EconomicLifetime_GEN[g])-Years[I],1),0)*CRF_GEN[g]*CAPEX_GEN[i0,g] for i0 = 1:I) for g = 1:GEN) )
# generator operating costs, fuel + VOM + FOM
@variable(m, costs_GENoperating[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_GENoperating[I] == sum(UnitSize_GEN[g]*(NumUnits_GEN[g]+sum(unitsbuilt_GEN[i0,g]-unitsretired_GEN[i0,g] for i0 = 1:I))*FOM_GEN[I,g] for g = 1:GEN)
                                                     + sum(weights[I,T]*8760/t_ops*sum(sum((VOM_GEN[I,g]+HeatRate[g]*FuelCosts[I,g]*(1-NG_fueled[g]))/1000*generation[I,T,t,g] for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops)
                                                     + sum(weights[I,T]*8760/t_ops*sum((StartUpCosts[g]+StartupFuel[g]*FuelCosts[I,g]*(1-NG_fueled[g]))/1000*sum(startup_GEN[I,T,t,g] for t = 1:t_ops) for g = 1:GEN) for T = 1:T_ops))
# elec storage capital costs, ammortized to each investment period
@variable(m, costs_ELECSTORAGEcapital[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_ELECSTORAGEcapital[I] == sum(UnitSize_STORAGE_ELEC[s]*sum(unitsbuilt_STORAGE_ELEC[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_ELEC[s])-Years[I],1),0)*CRF_STORAGE_ELEC[s]*CAPEX_STORAGE_ELEC[i0,s] for i0 = 1:I) for s = 1:STORAGE_ELEC) )
# elec storage operating costs, FOM + VOM
@variable(m, costs_ELECSTORAGEoperating[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_ELECSTORAGEoperating[I] == sum(UnitSize_STORAGE_ELEC[s]*(NumUnits_STORAGE_ELEC[s]+sum(unitsbuilt_STORAGE_ELEC[i0,s]-unitsretired_STORAGE_ELEC[i0,s] for i0 = 1:I))*FOM_STORAGE_ELEC[I,s] for s = 1:STORAGE_ELEC) )
#
### RONDO EDIT
# heat storage capital costs, ammortized to each investment period
@variable(m, costs_HEATSTORAGEcapital[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_HEATSTORAGEcapital[I] == sum(UnitSize_STORAGE_HEAT[s]*sum(unitsbuilt_STORAGE_HEAT[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_HEAT[s])-Years[I],1),0)*CRF_STORAGE_HEAT[s]*CAPEX_STORAGE_HEAT[i0,s] for i0 = 1:I) for s = 1:STORAGE_HEAT) )
# heat storage operating costs, FOM + VOM
@variable(m, costs_HEATSTORAGEoperating[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_HEATSTORAGEoperating[I] == sum(UnitSize_STORAGE_HEAT[s]*(NumUnits_STORAGE_HEAT[s]+sum(unitsbuilt_STORAGE_HEAT[i0,s]-unitsretired_STORAGE_HEAT[i0,s] for i0 = 1:I))*FOM_STORAGE_HEAT[I,s] for s = 1:STORAGE_HEAT) )
#
# POWER TO HEAT
# P2H operating costs
@variable(m, costs_P2Hoperating[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_P2Hoperating[I] == sum(UnitSize_P2H[d]*(NumUnits_P2H[d] + sum(unitsbuilt_P2H[i0,d]-unitsretired_P2H[i0,d] for i0 = 1:I))*FOM_P2H[I,d] for d = 1:P2H) 
                                                     + sum(weights[I,T]*8760/t_ops*sum(sum(VOM_P2H[I,d]/1000*P2H_dispatch[I,T,t,d] for t = 1:t_ops) for d = 1:P2H) for T = 1:T_ops) )
# P2H capital costs, ammortized to each investment period
@variable(m, costs_P2Hcapital[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_P2Hcapital[I] == sum(UnitSize_P2H[d]*sum(unitsbuilt_P2H[i0,d]*max(min((Years[i0]+EconomicLifetime_P2H[d])-Years[I],1),0)*CRF_P2H[d]*CAPEX_P2H[i0,d] for i0 = 1:I) for d = 1:P2H) )
#
### RONDO EDIT
#
# appliances costs
@variable(m, costs_appliances[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_appliances[I] ==  sum(sum(CRF_APPLIANCES[a]*max(min((Years[i0]+ApplianceLifetime[a])-Years[I],1),0)*(CAPEX_APPLIANCES[i0,a]*unitsbuilt_APPS[i0,a]*1000 + applianceInfrastructureCosts[i0,a]) for i0 = 1:I)/1000 for a = 1:APPLIANCES) )
# transport costs
@variable(m, costs_transport[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_transport[I]  ==  sum(sum(CRF_TRANSPORT[tr]*max(min((Years[i0]+TransportLifetime[tr])-Years[I],1),0)*(CAPEX_TRANSPORT[i0,tr]*unitsbuilt_TRANSPORT[i0,tr]*1000 + transportInfrastructureCosts[i0,tr] + nonElectricFuelCost[i0,tr])  for i0 =1:I) for tr = 1:TRANSPORTS))

# P2G operating costs
@variable(m, costs_P2Goperating[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_P2Goperating[I] == sum(UnitSize_P2G[d]*(NumUnits_P2G[d] + sum(unitsbuilt_P2G[i0,d]-unitsretired_P2G[i0,d] for i0 = 1:I))*FOM_P2G[I,d] for d = 1:P2G) + sum(weights[I,T]*8760/t_ops*sum(sum(VOM_P2G[I,d]/1000*P2G_dispatch[I,T,t,d] for t = 1:t_ops) for d = 1:P2G) for T = 1:T_ops) )
# P2G capital costs, ammortized to each investment period
@variable(m, costs_P2Gcapital[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_P2Gcapital[I] == sum(UnitSize_P2G[d]*sum(unitsbuilt_P2G[i0,d]*max(min((Years[i0]+EconomicLifetime_P2G[d])-Years[I],1),0)*CRF_P2G[d]*CAPEX_P2G[i0,d] for i0 = 1:I) for d = 1:P2G) )
# gas storage capital costs, ammortized to each investment period
@variable(m, costs_GASSTORAGEcapital[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_GASSTORAGEcapital[I] == sum(UnitSize_STORAGE_GAS[s]*sum(unitsbuilt_STORAGE_GAS[i0,s]*max(min((Years[i0]+EconomicLifetime_STORAGE_GAS[s])-Years[I],1),0)*CRF_STORAGE_GAS[s]*CAPEX_STORAGE_GAS[i0,s] for i0 = 1:I) for s = 1:STORAGE_GAS) )
# gas storage operating costs
@variable(m, costs_GASSTORAGEoperating[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_GASSTORAGEoperating[I] == sum(UnitSize_STORAGE_GAS[s]*(NumUnits_STORAGE_GAS[s]+sum(unitsbuilt_STORAGE_GAS[i0,s]-unitsretired_STORAGE_GAS[i0,s] for i0 = 1:I))*FOM_STORAGE_GAS[I,s] for s = 1:STORAGE_GAS) )
# distribution costs
@variable(m, costs_distribution[I = 1:T_inv])
@constraint(m, [I = 1:T_inv], costs_distribution[I] == Cost_DistributionInfrastructure*sum(PeakDistDemandInc[I,n] for n = 1:NODES_ELEC) )




### OBJECTIVE FUNCTION
@objective(m, Min, sum(discountfactor[i]*(
    # GEN
      costs_GENcapital[i]
    + costs_GENoperating[i]
    # ELEC Storage
    + costs_ELECSTORAGEcapital[i]
    + costs_ELECSTORAGEoperating[i]
    ### RONDO EDIT
    # Heat Storage
    + costs_HEATSTORAGEcapital[i]
    + costs_HEATSTORAGEoperating[i]
    # P2H
    + costs_P2Hcapital[i]
    + costs_P2Hoperating[i]
    ### RONDO EDIT
    # Appliances
    + costs_appliances[i]
    # Transport
    + costs_transport[i]
    # P2G
    + costs_P2Gcapital[i]
    + costs_P2Goperating[i]
    # GAS Storage
    + costs_GASSTORAGEcapital[i]
    + costs_GASSTORAGEoperating[i]
    # Distribution costs
    + costs_distribution[i]
    #
    + costs_NGimports[i]
    + gasdistsyst_Cost[i] + costs_CO2offsets[i]
    + costs_transmission[i] + costs_gasStorage[i]
    ) for i = 1:T_inv)
    )



status = optimize!(m)