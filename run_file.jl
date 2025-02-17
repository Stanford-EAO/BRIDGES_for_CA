import Pkg
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Clustering")
Pkg.add("Distances")
Pkg.add("Gurobi")
Pkg.add(Pkg.PackageSpec(;name="Gurobi", version="1.1.0"))
Pkg.add("Tables")
Pkg.add("DelimitedFiles")
Pkg.add("Dates")
Pkg.add("Grisu")
Pkg.add("Random")
Pkg.add("JuMP")
Pkg.add(Pkg.PackageSpec(;name="JuMP", version="1.7.0"))
Pkg.add("YAML")
# Pkg.add("NBInclude")
# Pkg.add("Plots")

using DataFrames, CSV, Tables, Clustering, Distances, Dates, Random, JuMP, YAML
using Gurobi
# using NBInclude, Plots

# Choose desired config file (CONFIG_FILE_NAME is used in Snakefile and parameters.jl)
if length(ARGS) == 0
    CONFIG_FILE_NAME = "config_default"
else
    CONFIG_FILE_NAME = ARGS[1]
end
lines = readlines("DataPreprocessing/Snakefile")
lineindex = findfirst(startswith("CONFIG_FILE_NAME ="), lines)
lines[lineindex] = "CONFIG_FILE_NAME = \"" * CONFIG_FILE_NAME * "\""
open("DataPreprocessing/Snakefile", "w") do file
    foreach(line -> println(file, line), lines)
end

# Read parameters from config file
include("core/Parameters/parameters.jl")

# Read data import file
include("core/data_imports.jl")

# Read clustering file
include("core/clustering.jl")

# Define optimization program
m = Model(optimizer_with_attributes(Gurobi.Optimizer,"Threads" => 30,"BarHomogeneous" => 1,"ScaleFlag"=>2, "FeasibilityTol"=> 0.005, 
<<<<<<< HEAD
"LogToConsole" => 1, "ScaleFlag" => 1, "OptimalityTol" => 0.001, "BarConvTol"=> 0.0001, "Method"=> 2, "Crossover"=> 0)) #"Presolve"=>2)) #, "NumericFocus"=>2, "Presolve"=>2))
=======
    "LogToConsole" => 1, "ScaleFlag" => 1,
    "OptimalityTol" => 0.001, "BarConvTol"=> 0.0001, "Method"=> 2, "Crossover"=> 0)) #"Presolve"=>2)) #, "NumericFocus"=>2, "Presolve"=>2))
>>>>>>> 33743b9 (Initial commit)

# Read constraint and optimize file
include("core/cons_capacity.jl")
include("core/cons_dispatch.jl")
include("core/cons_flow.jl")
include("core/cons_policy.jl")

include("core/optimize.jl")

# Read export file
include("core/data_exports.jl")

<<<<<<< HEAD
println("Success!")
=======
println("Success!")

>>>>>>> 33743b9 (Initial commit)
