#!/bin/bash
#SBATCH --job-name=el_ep_r1
#SBATCH --partition shared
#SBATCH --nodes=1 --ntasks=1 -c 12
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2G

ml bio/Beast

beast -threads 12 -prefix redo_elep/epoch/ redo_elep/epoch/elep_epoch_r1.xml