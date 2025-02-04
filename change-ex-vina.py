#! /Users/srusti/FergLab/md-tools/bin/python3

## Remaining code from changing exhaustiveness parameters for run_vina function

from vina import Vina # for autodock vina
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda # for defining ligand box
import os

def run_vina(receptor_file, ligand_file, output_path, center, size, exhaustiveness=exhaustiveness): #8 is default ex
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
    output_file = os.path.join(output_path, f"sds_bcd_mmff94_out.pdbqt")
    v.write_poses(output_file, n_poses=5, overwrite=True)  # Save the top 5 poses to output
    
    print(f"Docking completed for box size {size}.")

    # Retrieve the binding affinity of the top pose
    affinities = v.energies()  # Get binding affinities for all poses
    
    if affinities is not None and len(affinities) > 0:
        top_affinity = affinities[0][0]  # Extract the top pose affinity
    else:
        top_affinity = None

    print(f"Top binding affinity: {top_affinity} kcal/mol")
    return top_affinity

# "True" dG from molecular simulations
dG_bcd = np.array([-27.96716669556723])

# Initialize lists to store binding affinities and box sizes
exs = []
binding_affinities = []

# Run Vina with increasing exhuastiveness
for ex in range(4, 16):  # Increase exhaustivess by 1 for each run
    if ex > 15:  
        break
    top_affinity = run_vina(receptor, ligand, output_path, center, size, ex)
    if top_affinity is not None:
        binding_affinities.append(top_affinity)
        exs.append(ex)

# Plotting the results
plt.figure(figsize=(10, 6))

# Scatter plot: Box size (x, y, z dimension increase) vs Binding Affinity
plt.scatter([e for e in exs], binding_affinities, color='b', label='Predicted Binding Affinity (Vina)', s=100)

# Labels and Title
plt.xlabel('Exhaustiveness', fontsize=14)
plt.ylabel('Binding Affinity (kcal/mol)', fontsize=14)
plt.title('SDS-bCD Binding Affinity vs Exhaustiveness', fontsize=16)

# Show plot
plt.grid(True)
plt.show()