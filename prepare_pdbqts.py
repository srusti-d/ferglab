import os
import subprocess

input_dir = "/Users/srusti/FergLab/data/md_bcd_mods/relaxed-structures"

for filename in os.listdir(input_dir):
    if filename.endswith(".pdb"):
        base = os.path.splitext(filename)[0]
        pdb_file = os.path.join(input_dir, filename)
        pdbqt_file = os.path.join(input_dir, f"{base}.pdbqt")
        atom_file = os.path.join(input_dir, f"{base}_ATOMs.txt")

        # Convert .pdb to .pdbqt
        subprocess.run(["obabel", pdb_file, "-O", pdbqt_file], check=True)

        # grep ATOM lines from .pdbqt and write to new file
        with open(pdbqt_file, "r") as infile, open(atom_file, "w") as outfile:
            for line in infile:
                if "ATOM" in line:
                    outfile.write(line)
            outfile.write("TER\n")

        print(f"Finished {filename}")

