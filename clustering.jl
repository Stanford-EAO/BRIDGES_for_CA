################################################################################
### Clustering time series data for operational simulations
################################################################################
# Here we hope to allow for computationally tractable simulations of grid operations
# by providing the option to cluster the full 8760 hourly profile down into
# a set of representative days. K-mediods is employed.
# Reshape your load vector for your zone of interest into lengths of HOURS_PER_PERIOD

DemandClustering = copy(D_Elec)
for n = 1:NODES_ELEC
    DemandClustering[:,n] = DemandClustering[:,n] +  sum(APPLIANCES_NodalLoc_ELEC[n,a]*ApplianceProfilesELEC[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)
end
DemandClustering = sum(DemandClustering, dims = 2)
DemandClustering = (DemandClustering.- minimum(DemandClustering))./(maximum(DemandClustering) - minimum(DemandClustering))
LOAD_ELEC = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))
DemandClustering = copy(D_Gas)
for n = 1:NODES_GAS
    DemandClustering[:,n] = DemandClustering[:,n] +  sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES)
end
DemandClustering = sum(DemandClustering, dims = 2)
DemandClustering = (DemandClustering.- minimum(DemandClustering))./(maximum(DemandClustering) - minimum(DemandClustering))
LOAD_GAS = reshape(DemandClustering, (HOURS_PER_PERIOD, Periods_Per_Year))
ProfilesClustering = unique(HourlyVRE, dims = 2)
HourlyVREProfilesClustering = reshape(ProfilesClustering, (HOURS_PER_PERIOD, Periods_Per_Year,length(ProfilesClustering[1,:])))
global ClusteringData = vcat(LOAD_ELEC,  LOAD_GAS)
for x = 1:length(ProfilesClustering[1,:])
    global ClusteringData = vcat(ClusteringData, HourlyVREProfilesClustering[:,:,x])
end


ClustNatgas_mat = zeros(1,Periods_Per_Year)
# remove the magic numbers
price_Natgas = CSV.read("$(foldername)/NatgasPrice_perSector_2019.csv",DataFrame)
NGprices_timeSeries = zeros(HOURS_PER_YEAR,1)
if toggle_variableNatGasPrice ==  true
    for i = 3:length(price_Natgas[1,:])
        localClustNatgas_vector = transpose( price_Natgas[:,i] )
        localClustNatgas_mat = localClustNatgas_vector
        for j = 1:23
            localClustNatgas_mat = vcat(localClustNatgas_mat, localClustNatgas_vector)
        end
        # print(size(reshape(localClustNatgas_mat, (HOURS_PER_YEAR,1))))
        NGprices_timeSeries[:,i-2] = reshape(localClustNatgas_mat, (HOURS_PER_YEAR,1))
        global ClustNatgas_mat = vcat(ClustNatgas_mat, localClustNatgas_mat)
    end    
elseif toggle_variableNatGasPrice == false
    for i = 3:length(price_Natgas[1,:])
        i_fixed = 3
        priceNatGas_fixed = 7
        localClustNatgas_vector = transpose(priceNatGas_fixed * ones(Periods_Per_Year,1))
        localClustNatgas_mat = localClustNatgas_vector
        for j = 1:23
            localClustNatgas_mat = vcat(localClustNatgas_mat, localClustNatgas_vector)
        end
        global ClustNatgas_mat = vcat(ClustNatgas_mat, localClustNatgas_mat)
    end    
end
ClustNatgas_mat = (ClustNatgas_mat.- minimum(ClustNatgas_mat))./(maximum(ClustNatgas_mat) - minimum(ClustNatgas_mat))
ClustNatgas_mat = ClustNatgas_mat[1:size(ClustNatgas_mat,1) .!= 2,: ]
ClusteringData = vcat(ClusteringData, ClustNatgas_mat)


medoids = zeros(T_inv, T_ops)
RepDays = zeros(T_inv, Periods_Per_Year)
weights = zeros(T_inv, T_ops)
m = zeros(length(ClusteringData[:,1]),1)

D = pairwise(Euclidean(), ClusteringData, dims = 2)
H = hclust(D, linkage = :average)
if clustering_case =="ward"
    global H = hclust(D, linkage = :ward)
end
a = cutree(H,k = T_ops)

if clustering_case =="kmeans"
    global H = kmeans(ClusteringData, T_ops)
    global a = assignments(H)
end

for i = 1:T_inv
    for c = 1:(T_ops)
        sub = ClusteringData[:,a .== c]
        FindMedoids = pairwise(Euclidean(),sub,dims = 2)
        Klustering = kmedoids(FindMedoids, 1)
        m[:,1] = sub[:,Klustering.medoids[1]]
        C2 = pairwise(Euclidean(),ClusteringData,m)
        medoids[i,c] = argmin(C2,dims =1)[1][1]
        weights[i,c] = length(sub[1,:])/Periods_Per_Year
    end
    RepDays[i,:] = copy(a)
end
println("Check weights :")
println(weights[1,:])
println("Check medoids :")
println(medoids[1,:])

## After identifying the representative days, the raw data are
# re-shaped into the required indexing format for optimization
#################################################################
BaselineDemand_ELEC = zeros(T_inv, T_ops, t_ops, NODES_ELEC)
BaselineDemand_GAS = zeros(T_inv, T_ops, t_ops, NODES_GAS)
HourlyVRE_full = copy(HourlyVRE)
HourlyVRE = zeros(T_ops, t_ops, GEN)
ApplianceProfiles_GAS = zeros(T_ops, t_ops, APPLIANCES)
ApplianceProfiles_ELEC = zeros(T_ops, t_ops, APPLIANCES)

fuelPrice_res_GAS = zeros(T_inv, T_ops, t_ops)
fuelPrice_com_GAS = zeros(T_inv, T_ops, t_ops)

for t = 1:T_inv
    for i = 1:T_ops
        start_hour = Int((medoids[1,i]-1)*HOURS_PER_PERIOD+1)
        end_hour = start_hour + HOURS_PER_PERIOD-1
        for n = 1:NODES_ELEC
            BaselineDemand_ELEC[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Elec[start_hour:end_hour,n]
        end
        for n = 1:NODES_GAS
            BaselineDemand_GAS[t,i,:,n] = (1+LoadGrowthRate)^(Years[t]-BaseYear).*D_Gas[start_hour:end_hour,n]
        end
        for g = 1:GEN
            HourlyVRE[i,:,g] = HourlyVRE_full[start_hour:end_hour,g]
        end
        for a =1:APPLIANCES
            ApplianceProfiles_ELEC[i,:,a] = round.(ApplianceProfilesELEC[start_hour:end_hour,a], digits = 8)
            ApplianceProfiles_GAS[i,:,a] = round.(ApplianceProfilesGAS[start_hour:end_hour,a], digits = 8)
        end
        fuelPrice_res_GAS[t,i,:] = NGprices_timeSeries[start_hour:end_hour,1]
    end
end



################################################################################
### Produce residence matrices to identify which nodes contain which resources
################################################################################
# For each set of resources (GEN, DEM, STORAGE_ELEC, STORAGE_GAS, RNG, FS, CCS)
# produce a matrix that is of dimensions [NODES_ELEC x X] and [NODES_GAS x X]
# containing a 1 where each resource resides.
GEN_NodalLoc_ELEC = zeros(NODES_ELEC, GEN)
GEN_NodalLoc_GAS = zeros(NODES_GAS, GEN)
P2G_NodalLoc_ELEC = zeros(NODES_ELEC, P2G)
P2G_NodalLoc_GAS = zeros(NODES_GAS, P2G)
STORAGE_ELEC_NodalLoc_ELEC = zeros(NODES_ELEC, STORAGE_ELEC)
STORAGE_ELEC_NodalLoc_GAS = zeros(NODES_GAS, STORAGE_ELEC)
STORAGE_GAS_NodalLoc_ELEC = zeros(NODES_ELEC, STORAGE_GAS)
STORAGE_GAS_NodalLoc_GAS = zeros(NODES_GAS, STORAGE_GAS)
APPLIANCES_NodalLoc_ELEC = zeros(NODES_ELEC, APPLIANCES)
APPLIANCES_NodalLoc_GAS = zeros(NODES_GAS, APPLIANCES)

NodalLoc_ELEC = Generators[:,1]
NodalLoc_GAS = Generators[:,2]
for g = 1:GEN
    GEN_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[g])],REGIONS_ELEC)),g] = 1
    GEN_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[g])],REGIONS_GAS)),g] = 1
end
NodalLoc_ELEC = PowerToGas[:,1]
NodalLoc_GAS = PowerToGas[:,2]
for d = 1:P2G
 P2G_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[d])],REGIONS_ELEC)),d] = 1
 P2G_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[d])],REGIONS_GAS)),d] = 1
end
NodalLoc_ELEC = ElectricalStorage[:,1]
NodalLoc_GAS = ElectricalStorage[:,2]
for s = 1:STORAGE_ELEC
 STORAGE_ELEC_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[s])],REGIONS_ELEC)),s] = 1
 STORAGE_ELEC_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[s])],REGIONS_GAS)),s] = 1
end
NodalLoc_ELEC = GasStorage[:,1]
NodalLoc_GAS = GasStorage[:,2]
for s = 1:STORAGE_GAS
 STORAGE_GAS_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[s])],REGIONS_ELEC)),s] = 1
 STORAGE_GAS_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[s])],REGIONS_GAS)),s] = 1
end
NodalLoc_ELEC = EndUseAppliances[:,1]
NodalLoc_GAS = EndUseAppliances[:,2]
for a = 1:APPLIANCES
 APPLIANCES_NodalLoc_ELEC[findfirst(occursin.([string(NodalLoc_ELEC[a])],REGIONS_ELEC)),a] = 1
 APPLIANCES_NodalLoc_GAS[findfirst(occursin.([string(NodalLoc_GAS[a])],REGIONS_GAS)),a] = 1
end

################################################################################
### Produce topology matrices to identify which nodes are connected by the edges
################################################################################
# For each Flow (EDGES_ELEC and EDGES_GAS) produce an edge-nodal-incidence matrix A
# If simulating a single-node system, set A to all zeros, there is no transfer outside of the modeled region
A_ELEC = zeros(NODES_ELEC, EDGES_ELEC)
for e = 1:EDGES_ELEC
    A_ELEC[findfirst(occursin.([string(TransmissionLinks_ELEC[e,1])],REGIONS_ELEC)),e] = 1
    A_ELEC[findfirst(occursin.([string(TransmissionLinks_ELEC[e,2])],REGIONS_ELEC)),e] = -1
end
if NODES_ELEC == 1
    A_ELEC = A_ELEC*0
end
A_GAS = zeros(NODES_GAS, EDGES_GAS)
for e = 1:EDGES_GAS
    A_GAS[findfirst(occursin.([string(TransmissionLinks_GAS[e,1])],REGIONS_GAS)),e] = 1
    A_GAS[findfirst(occursin.([string(TransmissionLinks_GAS[e,2])],REGIONS_GAS)),e] = -1
end
if NODES_GAS == 1
    A_GAS = A_GAS*0
end

# Separately specify the cost of commodity natural gas for core demands
# using the same assumption applied to gas-fired electricity generators from FuelCosts.csv
################################################################################
CommodityCost_NG = 0.0*ones(T_inv,T_ops)      # $/MWh NG
for i = 1:T_inv
    for T = 1:T_ops
        scaleForRepDay_NGprice = fuelPrice_res_GAS[i,T,1] / 3.5
#         CommodityCost_NG[i,T] = sum(FuelCosts[i,:].*NG_fueled)/sum(NG_fueled)/MWh_PER_MMBTU * scaleForRepDay_NGprice
#         CommodityCost_NG[i,T] = sum(FuelCosts[i,:].*NG_fueled)/sum(NG_fueled)/MWh_PER_MMBTU
        CommodityCost_NG[i,T] = fuelPrice_res_GAS[i,T,1] / MWh_PER_MMBTU
    end
end

################################################################################
### Biomethane limitations
################################################################################
# Compute the maximum biomethane and bio-energy use [MWh/year] as a share of initial core gas demands
maxBiomethane = max_biomethane_share*sum(sum(sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES) + D_Gas[:,n] for n = 1:NODES_GAS))*ones(T_inv)               # MWh/yr
# The more generic use of sustainable biomass allows for eventual incorporation of bio-LPG which competes with biomethane for access to limited bioenergy feedstocks
maxSustainableBiomass = max_biomethane_share*sum(sum(sum(APPLIANCES_NodalLoc_GAS[n,a]*ApplianceProfilesGAS[:,a]*InitialAppliancePopulation[a] for a = 1:APPLIANCES) + D_Gas[:,n] for n = 1:NODES_GAS))*ones(T_inv)               # MWh/yr
