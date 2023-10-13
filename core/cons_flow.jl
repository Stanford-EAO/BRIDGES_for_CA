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
elseif TRANSMISSION_EXPANSION_ELEC == 3
    unitsbuilt_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
    unitsretired_TRANS_ELEC = zeros(T_inv,EDGES_ELEC)
    @variable(m, addflow_TRANS_ELEC[I = 1:T_inv,e = 1:EDGES_ELEC] >= 0)        # MW
    @constraint(m, [I = 1:T_inv, e = 1:EDGES_ELEC-2], addflow_TRANS_ELEC[I,e] ==  0)
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
## to track energy balances:
# P2G
@variable(m, GAS_energy_P2G[I = 1:T_inv, T = 1:T_ops])
@constraint(m, [I = 1:T_inv, T = 1:T_ops], GAS_energy_P2G[I,T] == sum(sum(P2G_dispatch[I,T,t,d]*eta_P2G[d] for d=1:P2G) for t=1:t_ops) )
# Imports
@variable(m, GAS_energy_Imports[I = 1:T_inv, T = 1:T_ops])
@constraint(m, [I = 1:T_inv, T = 1:T_ops], GAS_energy_Imports[I,T] == t_ops * sum(SUPPLY_GAS_slack[I,T,n] for n=1:NODES_GAS) )
# Charging
@variable(m, GAS_energy_Charging[I = 1:T_inv, T = 1:T_ops])
@constraint(m, [I = 1:T_inv, T = 1:T_ops], GAS_energy_Charging[I,T] == sum(sum(charging_GAS[I,T,t,s] for s = 1:STORAGE_GAS) for t = 1:t_ops) )
# Discharging
@variable(m, GAS_energy_Discharging[I = 1:T_inv, T = 1:T_ops])
@constraint(m, [I = 1:T_inv, T = 1:T_ops], GAS_energy_Discharging[I,T] == sum(sum(discharging_GAS[I,T,t,s] for s = 1:STORAGE_GAS) for t = 1:t_ops) )
# Demand
@variable(m, GAS_energy_Demand[I = 1:T_inv, T = 1:T_ops])
@constraint(m, [I = 1:T_inv, T = 1:T_ops], GAS_energy_Demand[I,T] == sum(sum(sum(NominalGasOfftakes[I,T,t,n,g]*LHV[g]*MolarMass[g] for g=1:GAS_COMPONENTS) for n=1:NODES_GAS) for t = 1:t_ops) )






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