#!/bin/bash
#SBATCH -J dntau                         # Job name
#SBATCH -n 1                            # Number of tasks (CPUs)
#SBATCH --mem=20G                       # Total memory
#SBATCH -o outfile_%j.log               # STDOUT (job ID in filename)
#SBATCH -e errorfile_%j.log             # STDERR
#SBATCH --partition=long              # Partition/queue (confirm it's valid)




# Compile
make clean
make


# Run the program
./radiative_transfer

