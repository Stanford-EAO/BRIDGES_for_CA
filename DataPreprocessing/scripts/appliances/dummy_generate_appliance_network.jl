using CSV
using DataFrames

# I tried to follow this julia style guide: https://github.com/invenia/BlueStyle

"""
    copy_profile(input_csv_path::AbstractString, output_csv_path::AbstractString, base_year::Int) -> None

Reads the csv file from input_csv_path into a dataframe and writes this dataframe to a csv file output_csv_path.

# Arguments:
- `input_csv_path::AbstractString`: the csv file to import
- `output_csv_path::AbstractString`: the path for saving the dataframe
- `base_year`: simulation base year (not needed for anything, just illustrates parameter import)

# Returns:

"""
function copy_profile(; input_csv_path::AbstractString, output_csv_path::AbstractString, base_year::Int)
    
    # Read CSV file into a DataFrame
    df = CSV.File(input_csv_path) |> DataFrame

    # Write DataFrame to CSV file in the output directory
    CSV.write(output_csv_path, df)

    # Display a message
    println("Some output from the dummy script for appliances:")
    println("Base year is: ", base_year)
    println("DataFrame written to: ", output_csv_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    copy_profile(

        input_csv_path=snakemake.input[1], # this would be [0] in python scripts, damn that took me 30 minutes :P.
        output_csv_path=snakemake.output[1],
        base_year=snakemake.params["base_year"]
    )
end
