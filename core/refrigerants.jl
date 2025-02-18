# Refrigerant leakages #
GWP_R410A = 2088
GWP_R32 = 675
GWP_propane = 3
GWP_CO2 = 1
lb_to_ton = 0.0004535924
annual_leak_pu = APP_avg_charge.*APP_annual_leak.*lb_to_ton       # per unit annual leak [ton refrigerant/year]
eol_loss_pu = APP_eol_charge.*APP_eol_loss.*lb_to_ton             # per unit eol loss [ton refrigerant]

# heatpump wh, res unitary AC, res heat pump space heaters
@variable(m, appliance_leak[I = 1:T_inv,a = 1:APPLIANCES] >= 0)     # MMTCO2

@variable(m, APP_historical_leak[I = 1:T_inv,a = 1:APPLIANCES] >= 0)
@variable(m, APP_historical_retired[I = 1:T_inv,a = 1:APPLIANCES] >= 0)
@constraint(m, [I = 1,a = 1:APPLIANCES], APP_historical_retired[I,a] == InitialAppliancePopulation[a]/1000 - unitsremaining_APPS_historical[I,a])
if T_inv > 1
    @constraint(m, [I = 2:T_inv,a = 1:APPLIANCES], APP_historical_retired[I,a] == unitsremaining_APPS_historical[I-1,a] - unitsremaining_APPS_historical[I,a])
end
@constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], APP_historical_leak[I,a] == GWP_R410A*(annual_leak_pu[a]*unitsremaining_APPS_historical[I,a] + eol_loss_pu[a]*APP_historical_retired[I,a])/1000) # MMT CO2e/year

@variable(m, APP_new_leak[I = 1:T_inv,a = 1:APPLIANCES] >= 0)
@variable(m, APP_new_remaining[I = 1:T_inv,a = 1:APPLIANCES] >= 0)
@constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], APP_new_remaining[I,a] == sum(unitsbuilt_APPS[i0,a] - unitsretired_APPS[i0,a] for i0 = 1:I))

if consider_refrigerants == 1       # CARB aligned
    @constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], APP_new_leak[I,a] == GWP_R32*(annual_leak_pu[a]*APP_new_remaining[I,a] + eol_loss_pu[a]*unitsretired_APPS[I,a])/1000) # MMT CO2e/year
    @constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], appliance_leak[I,a] == APP_historical_leak[I,a] + APP_new_leak[I,a])  # MMTCO2/unit
elseif consider_refrigerants == 2   # beyond CARB
    @constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], APP_new_leak[I,a] == GWP_propane*(annual_leak_pu[a]*APP_new_remaining[I,a] + eol_loss_pu[a]*unitsretired_APPS[I,a])/1000) # MMT CO2e/year
    @constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], appliance_leak[I,a] == APP_historical_leak[I,a] + APP_new_leak[I,a])  # MMTCO2/unit
elseif consider_refrigerants == 3   # stagnant case
    @constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], APP_new_leak[I,a] == GWP_R410A*(annual_leak_pu[a]*APP_new_remaining[I,a] + eol_loss_pu[a]*unitsretired_APPS[I,a])/1000) # MMT CO2e/year
    @constraint(m, [I = 1:T_inv,a = 1:APPLIANCES], appliance_leak[I,a] == APP_historical_leak[I,a] + APP_new_leak[I,a])  # MMTCO2/unit
else
    error("consider_refrigerants must be between 0 and 3!")
end