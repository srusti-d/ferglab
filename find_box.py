from pymol import cmd

def calc_bounding_box(selection='all'):
    min_x, min_y, min_z = [9999] * 3 # Initial values, later replaced by script
    max_x, max_y, max_z = [-9999] * 3
    stored.coords = [] # PyMol-specific dictionary used to store data across iterations

    # Get atom coordinates
    # Iterating over all atoms in the selection and running a small Python code snippet
    cmd.iterate_state(1, selection, "stored.coords.append([x, y, z])")
    
    for coord in stored.coords:
        x, y, z = coord
        min_x, min_y, min_z = min(min_x, x), min(min_y, y), min(min_z, z) # Find smallest coord along all axes
        max_x, max_y, max_z = max(max_x, x), max(max_y, y), max(max_z, z) # Find largest coord along all axes

    # Find center: midpoint of the bounding box along the x, y, and z axes
    center_x = (min_x + max_x) / 2
    center_y = (min_y + max_y) / 2
    center_z = (min_z + max_z) / 2

    # Find minimum size of atom: subtracting the minimum from the maximum coordinate values
    size_x = max_x - min_x
    size_y = max_y - min_y
    size_z = max_z - min_z
    
    print(f"Center: {center_x}, {center_y}, {center_z}")
    print(f"Size: {size_x}, {size_y}, {size_z}")

# Usage on PyMol terminal
# load /Users/srusti/FergLab/data/bcd.pdb  # Load CD PDB file
# run /Users/srusti/FergLab/test_autodock/find_box.py
# calc_bounding_box('all')  # Run calculation with all atoms as input values

## Results from PyMol for minimum box size of CD
# Center: 29.960000038146973, 30.285000801086426, 30.0600004196167
# Size: 14.719999313354492, 14.510000228881836, 6.059999465942383

## Results from PyMOl for minimum box size of SDS
# Center: 30.029000282287598, 29.998000144958496, 29.662999153137207
# Size: 2.460000991821289, 2.469999313354492, 18.469999313354492

