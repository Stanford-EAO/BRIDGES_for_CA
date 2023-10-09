# import Pkg
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Clustering")
# Pkg.add("Distances")
# Pkg.add("Gurobi")
# Pkg.add("Tables")
# Pkg.add("DelimitedFiles")
# Pkg.add("Dates")
# Pkg.add("Grisu")
# Pkg.add("Random")
# Pkg.add("JuMP")
# Pkg.add(Pkg.PackageSpec(;name="JuMP", version="1.7.0"))

# Pkg.add("NBInclude")
# Pkg.add("Plots")
using DataFrames, CSV, Tables, Clustering, Distances, Gurobi, Dates, Random
# using NBInclude, Plots
using JuMP

# Read parameter file
include("core/parameters.jl")

# Read data import file
include("core/data_imports.jl")

# Read clustering file
include("core/clustering.jl")

# Define optimization program
m = Model(optimizer_with_attributes(Gurobi.Optimizer,"Threads" => 46,"BarHomogeneous" => 1,"ScaleFlag"=>2, "FeasibilityTol"=> 0.005, 
    "OptimalityTol" => 0.001, "BarConvTol"=> 0.0001, "Method"=> 2, "Crossover"=> 0, "NumericFocus"=>2, "Presolve"=>2))

# Read constraint and optimize file
include("core/cons_capacity.jl")
include("core/cons_dispatch.jl")
include("core/cons_flow.jl")
include("core/cons_policy.jl")

include("core/optimize.jl")

# Read export file
include("core/data_exports.jl")

println("Success!")

