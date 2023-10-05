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
    # ELECTRIC
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
    
    # GAS
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
    # end of year == beginning of year
    @constraint(m, [I = 1:T_inv, s = 1:STORAGE_GAS], SOCTracked_GAS[I,Int(Periods_Per_Year),s] == storedEnergy_GAS[I,Int(RepDays[I,1]),1,s])

    # define initial conditions for storage
    # gas: for day 1 in investment period 1, start with that
    @constraint(m, [I = 1, d = 1, s = 1:STORAGE_GAS], SOCTracked_GAS[I,d,s] == initialStorage_GAS[s])

    # electric: start at zero or at max capacity
    @constraint(m, [I = 1, d = 1, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,d,s] == 0)
    # @constraint(m, [I = 1, d = 1, s = 1:STORAGE_ELEC], SOCTracked_ELEC[I,d,s] == (UnitSize_STORAGE_ELEC .* duration_ELEC)[s])

end