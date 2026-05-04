# radiative_transfer_classIIdisk_SFR
 Disk Modeling Pipeline
Input File Generation & Radiative Transfer Workflow
README — Complete User Guide
1. Overview
This pipeline generates disk model input files for a set of protoplanetary disk parameter configurations, and subsequently runs a multi-physics radiative transfer code on each disk. The workflow consists of two sequential stages:

    • Stage 1 — Input File Generation: reads disk parameters from a TSV file and produces structured input files for each disk using the disk_model executable.
    • Stage 2 — Radiative Transfer: runs the radiative_transfer executable inside each disk's folder, reading those input files and producing temperature/flux outputs.

2. Directory Structure
The pipeline expects and produces the following directory layout:

Path
Description
Disk_modelling_code/
Root working directory
Disk_modelling_code/disk_parameters.tsv
Parameter file — one row per disk
Disk_modelling_code/results/
Output root — one subfolder per disk
Disk_modelling_code/results/disk_N/
All input & output files for disk N
Disk_modelling_code/logs/
Logs from Stage 1 (input generation)
Disk_modelling_code/logs_rt/
Logs from Stage 2 (radiative transfer)
Disk_modelling_code/radiative_transfer/
Fortran source code & compiled binary
Disk_modelling_code/radiative_transfer/radiative_transfer
Compiled executable

3. Prerequisites
3.1 Required Software
    • gfortran (GNU Fortran compiler, version 9+)
    • make
    • bash (version 4+)

3.2 Required Files Before Running
    • disk_parameters.tsv — parameter file in the root directory
    • disk_model executable — compiled disk model binary in root directory
    • Fortran source files — in the radiative_transfer/ subdirectory
    • Opacity files (.Kappa, dustopac*.inp) — in the root directory

4. Stage 1 — Input File Generation
The script run_all_disks_simple.sh loops over a range of disk IDs, pipes each ID into the disk_model executable, and collects the generated input files into individual disk directories under results/.

For each disk ID, the following files are generated and moved to results/disk_N/:
    • dust_dens.inp — dust density structure
    • starinfo.inp — stellar parameters
    • frequency.inp — frequency grid


Compilation
The disk_model binary must be compiled before running Stage 1. From the root directory:

make

Running Stage 1
The script is located at the root of Disk_modelling_code/. Run it as follows:

chmod +x run_all_disks_simple.sh          # first time only
./run_all_disks_simple.sh                 # runs all disk IDs (default)
./run_all_disks_simple.sh 1 5            # run IDs 1 to 5 only
./run_all_disks_simple.sh 5 10          # run IDs 5 to 10 only

Note: Always test with a single disk first (e.g. ./run_all_disks_simple.sh 1 1) before running the full batch.

4.4 Script Parameters

Parameter
Description
START_ID
First disk ID to process (default: 1)
END_ID
Last disk ID to process (default: 10)
PARAM_FILE
Path to disk_parameters.tsv (default: disk_parameters.tsv)
EXECUTABLE
Path to disk_model binary (default: ./disk_model)
OUTPUT_DIR
Root output directory (default: results/)
LOG_DIR
Log directory (default: logs/)

4.5 Output After Stage 1
After Stage 1 completes, the results/ folder will contain one subfolder per disk:

results/disk_1/    radius.inp  dust_dens.inp  starinfo.inp  ...
results/disk_2/    radius.inp  dust_dens.inp  starinfo.inp  ...
results/disk_10/  radius.inp  dust_dens.inp  starinfo.inp  ...

Tip: If a disk fails, check its log at logs/disk_N.log for the error message. For some compact disks you have to run these disks manually one by one.

5. Stage 2 — Radiative Transfer
The script radiative_transfer_bash.sh enters each disk's results directory and runs the radiative_transfer executable in-place. Since the code reads input files from its current working directory, running it inside each disk folder allows it to find radius.inp, dust_dens.inp, dustopac*.inp, and all other required files automatically.

The Fortran code consists of several physics modules that run in sequence:
    • configure — reads model configuration and disk parameters
    • montecarlo — Monte Carlo radiative transfer for temperature computation
    • diffusion — flux diffusion in optically thick regions
    • vertical_structure — computes vertical hydrostatic equilibrium

 Compilation
Compile the radiative transfer code before running Stage 2:

cd /media/admin1/DATA/Input_Params_bashscript/Disk_modelling_code/radiative_transfer (change the path according to your directory)
make

This uses the following compiler settings from the Makefile:

Setting
Value
Compiler
gfortran
Flags
-O2 -g -Wall -fcheck=all -cpp
Output binary
radiative_transfer

Running Stage 2
Save the radiative_transfer_bash.sh script inside Disk_modelling_code/ and run:

chmod +x radiative_transfer_bash.sh          # first time only
./radiative_transfer_bash.sh 1 1             # test with disk 1 first
./radiative_transfer_bash.sh                 # run all disks if test passes

Key Paths in the Script
Variable
Value
RESULTS_DIR
/media/admin1/DATA/Input_Params_bashscript/Disk_modelling_code/results
RT_EXECUTABLE
/media/admin1/DATA/Input_Params_bashscript/Disk_modelling_code/radiative_transfer/radiative_transfer
LOG_DIR
/media/admin1/DATA/Input_Params_bashscript/Disk_modelling_code/logs_rt


Tip: You have to change the paths of the above three variables before compiling and run the code 
5.5 Expected Output Files
After a successful radiative transfer run, new output files will appear inside each disk folder alongside the input files. 

Typical output files produced by the radiative transfer code include:
    • dusttemp_final.dat — dust temperature as a function of grid position
    • spectrum_all.dat —  spectra information


6. Troubleshooting
6.1 Common Errors

Error
Cause & Fix
Directory disk_N not found
results/ folder missing or wrong path. Check RESULTS_DIR in the script.
Exit code 13 / Dust species count mismatch
Opacity files (dustopac_1.inp, dustopac_2.inp, *.Kappa) missing from disk folder.
Executable not found
Code not compiled yet. Run make inside radiative_transfer/.
Parameter file not found
disk_parameters.tsv missing. Place it in the root Disk_modelling_code/ directory.
Permission denied
Script not executable. Run: chmod +x <scriptname>.sh


6.2 Checking Logs
    • Stage 1 logs: logs/disk_N.log
    • Stage 2 logs: logs_rt/rt_disk_N.log






7. Quick Reference — Full Workflow

Step
Command
1. Compile disk model
cd Disk_modelling_code && make
2. Generate input files
./run_all_disks_simple.sh 1 10
3. Compile RT code
cd radiative_transfer && make
4. Run radiative transfer
./radiative_transfer_bash.sh 1 10
5. Check results
results/disk_1/
6. View logs on error
 logs_rt/rt_disk_N.log

Important: Always compile both executables before running either script. Always test on a single disk (ID 1) before launching the full 10-disk batch. You can specify any number of disks in the disk_parameters.tsv file and run the entire code together.
