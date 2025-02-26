import pandas as pd

# I tried to follow this python style guide: https://peps.python.org/pep-0008/ 


def copy_profile(input_csv_path: str,
                 output_csv_path: str,
                 base_year: int) -> None:

    """
    Reads the csv file from input_csv_path into a dataframe and writes this dataframe to a csv file output_csv_path.

    Parameters:
        input_csv_path (str):
            the csv file to import
        output_csv_path (str):
            the path for saving the dataframe
        base_year (int):
            simulation base year (not needed for anything, just illustrates parameter import)
    """

    # Read CSV file into a DataFrame
    df = pd.read_csv(input_csv_path)

    # Write DataFrame to CSV file
    df.to_csv(output_csv_path, index=False)

    # Display a message
    print(f"Some output from the dummy script for storage:")
    print(f"Base year is: {base_year}")
    print(f"DataFrame written to: {output_csv_path}")

if __name__ == "__main__":
    copy_profile(
        input_csv_path=snakemake.input[0],
        output_csv_path=snakemake.output[0],
        base_year=snakemake.params["base_year"]
    )
