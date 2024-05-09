DemandClustering = copy(D_Elec)
for n = 1:NODES_ELEC
    DemandClustering[:,n] = DemandClustering[:,n] +  sum(APPLIANCES_NodalLoc_ELEC[n,a]*ApplianceProfilesELEC[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)
end
DemandClustering = sum(DemandClustering, dims = 2)
LOAD_ELEC = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))
vals = zeros(365)
for i = 1:365
    vals[i],day = findmax(LOAD_ELEC[:,i])
end
# p1 = plot(vals)
x,max_day_elec = findmax(vals[:])
# vline!(p1, [max_day_elec],line=(:dash, 2),color=(:black))
# vline!(p1, medoids,color=(:orange))


DemandClustering = copy(D_Gas)
for n = 1:NODES_GAS
    DemandClustering[:,n] = DemandClustering[:,n] +  sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)
end
DemandClustering = sum(DemandClustering, dims = 2)
LOAD_GAS = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))
valsGAS = zeros(365)
for i = 1:365
    valsGAS[i],day = findmax(LOAD_GAS[:,i])
end
# p2 = plot(valsGAS)
x,max_day_gas = findmax(valsGAS[:])
# vline!(p2, [max_day_gas],line=(:dash, 2),color=(:black))
# vline!(p2, medoids,color=(:orange))

# plot(p1, p2, layout = (2, 1),legend=false)

println("Peak electric demand day = ", max_day_elec)
println("Peak gas demand day = ", max_day_gas)

println("Check weights with extreme days:")
unnorm_weights = weights.*365
valweights,day_ind = findmax(unnorm_weights[1,:])
unnorm_weights[:,day_ind] = ones(5)*(valweights-2)
weights = [unnorm_weights ones(5) ones(5)]./365
println(weights[1,:])
println("Check medoids with extreme days:")
medoids = [medoids ones(5)*max_day_gas ones(5)*max_day_elec]
println(medoids[1,:])

N_Periods = N_Periods + 2
T_ops = N_Periods   

# update RepDays
mm = findall(x->x==day_ind, RepDays[1,:])
ind_summer = mm[Int64(round(length(mm)/2))]      # summer
ind_winter = mm[Int64(round(length(mm)))]        # winter
RepDays[:,ind_summer] = ones(5)*N_Periods
RepDays[:,ind_winter] = ones(5)*(N_Periods-1)