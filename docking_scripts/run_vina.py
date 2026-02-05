#! /Users/srusti/FergLab/md-tools/bin/python3

from vina import Vina # for autodock vina
import MDAnalysis as mda # for defining ligand box

def run_vina(receptor_file, ligand_file, output_file, center, size, exhaustiveness=8):
    # Create a Vina object
    v = Vina(sf_name='vina')

    # Load receptor and set up the docking grid
    v.set_receptor(receptor_file) # receptor = cyclodrextrin

    # Load ligand
    v.set_ligand_from_file(ligand_file)
    print("ligand worked")

    # Set the center and size of the docking box
    v.compute_vina_maps(center=center, box_size=size)
    
    # Perform the docking
    v.dock(exhaustiveness=exhaustiveness, n_poses=10)  # Perform docking with 10 poses
    
    # Save the docked poses to output
    v.write_poses(output_file, n_poses=5, overwrite=True)  # Save the top 5 poses to output
    
    print(f"Docking completed.")

# Test Vina under the assumption of "blind docking" (not knowing binding site)

receptor = "/Users/srusti/FergLab/data/clean_bcd_start.pdbqt"  # Path to receptor file (CD) in PDBQT format
ligand = "/Users/srusti/FergLab/data/sds.pdbqt"  # Path to ligand file (PFAS) in PDBQT format
output = "/Users/srusti/FergLab/test_autodock/results/box/sds/min_box.pdbqt"      # Path to save output poses in PDBQT format
center = [30.029000282287598, 29.998000144958496, 29.662999153137207]  # Center of the grid box (x, y, z coordinates)
                             # PyMol "center" coordinates for start CD based on box finding algorithm
size = [18, 18, 18]    # Size of the grid box (x, y, z dimensions)
                             # Start with 60, 60, 60, 60 as the maximum box size for system

run_vina(receptor, ligand, output, center, size)
# %%

## Results from PyMol for minimum box size of CD
# Center: 29.960000038146973, 30.285000801086426, 30.0600004196167
# Minimum Size: 14.719999313354492, 14.510000228881836, 6.059999465942383

## Results from PyMol for minimum box size of SDS
# Center: 30.029000282287598, 29.998000144958496, 29.662999153137207
# Size: 2.460000991821289, 2.469999313354492, 18.469999313354492
