# Running BRIDGES

The goal of this section is to provide **instructions on running BRIDGES**. BRIDGES can be run on your **local machine** or on a **high-performance computing cluster**. For running the full optimization model at a high temporal resolution (> 2 investment periods, > 2 representative days) we strongly recommend using a computing cluster. Moreover, we distinguish between running the **data preprocessing pipeline**, i.e., the Snakemake workflow that generates the input data for the optimization model, and running the **julia optimization model** (see section "Project structure").

## Running BRIDGES on your computer

Start with copying the `BRIDGES_for_CA` repository to your computer.

### Running the data preprocessing pipeline

To run the data preprocessing pipeline and create the input data to the optimization model in `BRIDGES_for_CA/Data/`, follow the following steps.

1. Make sure you have a recent version of **Python** (Snakemake requires Python >= 3.5) and **Julia** installed. Install otherwise.
2. Activate (or create) the virtual/conda environment you are using for BRIDGES (if you use one).
3. Open a terminal and install Snakemake using one of the following commands, depending on whether you prefer Anaconda or pip as Python package manager. For additional tips and optional further installations, read the section "Introduction to Snakemake > Installation".
```
conda install -c bioconda -c conda-forge snakemake-minimal
```
or
```
pip install snakemake
```
4. Have a look at the folder `BRIDGES_for_CA/DataPreprocessing`. This folder contains python and julia scripts. Snakemake takes care of creating isolated environments in which the python scripts are executed. However, for julia scripts such a feature does not exist. Therefore, have a look at the **julia** scripts in `BRIDGES_for_CA/DataPreprocessing` and check at the top of the script files for required packages. Install the packages using the following commands in your terminal. 
```
julia                   # Starts julia.
]                       # Type ] to enter package mode.
status                  # Prints current package list.
add SomePackageName     # Type this to install a package. Replace "SomePackageName" with the package name.
status                  # To double-check that the package has been installed.
<Backspace>             # Hit backspace to exit package mode.
exit()                  # Exits julia.
```
(If you wish, julia also has environments similar to virtual/conda environments to install packages in. See [this video](https://www.youtube.com/watch?v=S91mAyow2tw) (or any other) for a short intro.)

5. Next, have a look at the file named `Snakefile`. It is located in the folder `BRIDGES_for_CA/DataPreprocessing`. At the top of the Snakefile, adjust the `BRIDGES_path` which is the location of the folder `BRIDGES_for_CA` on your local computer. In my case the path is `BRIDGES_path = "C:/Users/mheyer/Code_Stanford/BRIDGES_for_CA/"`.
6. In the terminal, navigate to the folder `BRIDGES_for_CA/DataPreprocessing` using "cd" commands. This folder should contain the file named `Snakefile`.
7. To execute the entire data preprocessing pipeline, run
```
snakemake run_all --cores all --use-conda --conda-frontend conda
```
Expect this to take a few moments or minutes. This should create the input files to the optimization model in `BRIDGES_for_CA/Data`. If you are prompted "Nothing to be done.", the input files already exist.

For a deeper understanding of Snakemake, the Snakemake command above, and other commands to run, please read "Introduction to Snakemake".

!!! note

    At the time of writing this documentation, an incompatibility between snakemake and the newest version of the package "pulp" exists. For now, proceed without worrying but if you run into issues, come back and simply downgrade pulp using the command `pip install --force-reinstall -v "pulp==2.7.0"`.

### Running the julia optimization model 

This assumes that all required input files exist in the folder `BRIDGES_for_CA/Data`, i.e., the data preprocessing pipeline has been executed (see above).

1. Make sure you have a recent version of **Julia** installed. Install otherwise. (If you want to execute the actual optimization on your personal computer, also install Gurobi. See "Tip" below.)
2. In a terminal window, navigate to the root folder `BRIDGES_for_CA`.
3. Start julia by typing `julia`.
4. Have a look at the file `BRIDGES_for_CA/run_file.jl`. Make sure the `Pkg.add("SomePackageName")` commands at the top of the file are active and not commented. (If you are certain that the listed packages are already installed in your julia environment, they don't need to be active - this saves computation time. You can check for installed julia packages by running `]` to enter the package manager and `status` to get a list of the installed packages.)
5. To run the optimization model, type the following in the julia command line:
```
include("run_file.jl")
```
This will run the optimization and generate the output files.

!!! Tip

    In practice, we almost never run the actual optimization on our personal computer but rather on the cluster. On our personal computer we test the `run_file.jl` up to `include("core/clustering.jl")` and deactivate the following commands using comments. In this case the installation of Gurobi (and the line `Pkg.add("Gurobi")` in the `run_file.jl`) is not required.

## Running BRIDGES on Stanford's Sherlock cluster

There exist two ways to connect to Sherlock. Either through **Sherlock's web interface** or by **connecting VS Code and Sherlock** so that the folder structure on Sherlock appears as "normal" workspace folder in VS code. Learn how to connect to Sherlock by reading "Introduction to Sherlock". Then, follow the steps below to run BRIDGES.

### Running the data preprocessing pipeline and the optimization model

1. Connect to Sherlock as described in "Introduction to Sherlock".
2. Copy the `BRIDGES_for_CA` repository to your Sherlock Home Directory.
3. Have a look at the file `BRIDGES_for_CA/run_file.jl`. Make sure all lines are active. After running the model once, the `Pkg.add("")` lines can be deactivated using comments - this saves computation time.
4. Now, have a look at the file `BRIDGES_for_CA/my_job.script`. In the line `#SBATCH --mail-user=XYZ@stanford.edu` include your email address. Return to this file and adjust the `#SBATCH` options, if your run fails because of insufficient computing power (you will know from the "slurm files" as described below). The file should look something like this (**note that there might have been updates to the file since this documentation was written**):
``` title="BRIDGES_for_CA/my_job.script"
#!/bin/bash

#SBATCH --job-name=BRIDGES
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your_username@stanford.edu

#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=40GB
#SBATCH --cpus-per-task=24
#SBATCH --partition=serc


# Load modules here
ml load julia               # required for both the data preprocessing and the optimization
ml load gurobi              # required for the optimization
ml load python/3.9.0        # Snakemake needs Python >3.5. Activate this line, if your Sherlock uses an older version of Python (check with "python --version").

# Run the data preprocessing pipeline
snakemake -s "./DataPreprocessing/Snakefile" --cores all --use-conda --conda-frontend conda     # executes the data preprocessing pipeline that creates the input files for the optimization

# Run the optimization model
srun julia run_file.jl
```
5. In the same file: If you want to run only the data preprocessing pipeline or only the optimization model, use comments to deactivate the corresponding commands: 
```
# Run the data preprocessing pipeline
snakemake -s "./DataPreprocessing/Snakefile" --cores all --use-conda --conda-frontend conda     # deactivate this line to only run the optimization model (the input files must be present)

# Run the optimization model
srun julia run_file.jl                                                                          # deactivate this line to only run the data preprocessing pipeline
```
Be advised that the files generated by the data preprocessing pipeline are required by the optimization model.
6. Open a terminal on Sherlock as described in "Introduction to Sherlock".
7. Install Snakemake to your Sherlock home directory by running:
```
module load python/3.9.0                    # or another python version >3.7. Check which ones are available on Sherlock by typing "module av".
pip3 install --user snakemake==7.32.4
``` 
8. In the terminal, navigate to the `BRIDGES_for_CA` folder using "cd" commands. (Only if you copied the BRIDGES repository from a Windows computer to Sherlock: Run the command `dos2unix my_job.script` to convert the line break formatting of this file to unix style.)
7. Execute BRIDGES by running the command `sbatch my_job.script`. This will send the job to the computing cluster. You will receive an email when your job begins and terminates. Upon termination you will find a "slurm file" in `BRIDGES_for_CA` containing potential error messages and all print-outs of the scripts. Moreover, in the folder `BRIDGES_for_CA/Output` the output files of the model are created. Use these for further analyses and plotting.

!!! note

    * At the time of writing this documentation, an incompatibility between snakemake and the newest version of the package "pulp" exists. For now, proceed without worrying but if you run into issues, come back and simply downgrade pulp using the command `pip install --force-reinstall -v "pulp==2.7.0"`.
    * Running Snakemake on the cluster with the flag `--use-conda` (see "Introduction to Snakemake" why that is necessary) requires the installation of `conda`, `miniconda`, or `mamba` on Sherlock. Generally, the high-performance computing team [discourages the installation of conda/miniconda/mamba on the cluster](https://www.sherlock.stanford.edu/docs/software/using/anaconda/), however, after reaching out, the team confirmed that for our use case the installation of `miniconda` (since it is smaller than `conda`) makes sense. They recommended installation in the `$GROUP_HOME` directory as there is more storage available and the software is available for all users of that group. **For Adam Brandt's group's Sherlock folder, miniconda is now installed - so no need for action for his students.**
    * Other energy system models are structured so that they are run by a single snakemake command instead of separate commands for data preprocessing and optimization. For this case, Snakemake provides a command to include all information from the "my_job.script" file into the snakemake command instead of having a script file with the snakemake command in it. See [here](https://hackmd.io/@bluegenes/BJPrrj7WB) for more information. 
