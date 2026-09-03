#!/bin/bash

# email on start, end, and abortion
#SBATCH --job-name=reweight
##SBATCH --output=out.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=srustidonapati@rcc.uchicago.edu
#SBATCH --partition=andrewferguson
#SBATCH --account=pi-andrewferguson
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
##SBATCH --mem=0 # use max memory on node
#SBATCH --time=2:00:00
##SBATCH --exclude=midway2-0661

# Usage: sbatch reweight.sbatch <probe> <analyte>
# e.g.   sbatch reweight.sbatch 00528 pfoa
# Falls back to parsing the system directory name (<probe>-<analyte>).
PROBE="${1:-}"
ANALYTE="${2:-}"

if [[ -z "$PROBE" || -z "$ANALYTE" ]]; then
    sysdir=$(basename "$PWD")
    PROBE="${sysdir%-*}"
    ANALYTE="${sysdir##*-}"
    echo "No probe/analyte given; inferred from directory name: probe=$PROBE analyte=$ANALYTE"
fi

echo "Reweighting probe=$PROBE analyte=$ANALYTE in $PWD"

bash reweight.sh
# dG-plot.py takes: <probe> <analyte>
python ../dG-plot.py "$PROBE" "$ANALYTE"
