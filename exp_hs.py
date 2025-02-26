#! /Users/srusti/FergLab/md-tools/bin/python3

from rdkit import Chem

# ## FIX ATOMS WITH IMPLICIT HYDROGENS (no OH avoidance)

# # Load the molecule, ensuring we don't remove explicit Hs
# mol = Chem.MolFromMolFile("carbon_only_h_sds.sdf", removeHs=False)
# # Force all implicit hydrogens to be converted to explicit
# mol = Chem.AddHs(mol, addCoords=True)
# # Save the corrected molecule
# Chem.MolToMolFile(mol, "exp_c_only_h_sds.sdf")
# print("Fixed SDF file.")


# # ADD EXPLICIT HYDROGENS ON SDS + AVOID HYDROGENS ON OXYGENS

# mol = Chem.MolFromMolFile("ligands/og_sds.sdf", removeHs=False)

# for atom in mol.GetAtoms():
#     if atom.GetSymbol() == "O":
#         atom.SetNoImplicit(True)  
#         atom.SetNumExplicitHs(0)  # Forces oxygen to have zero explicit Hs

# mol = Chem.AddHs(mol, addCoords=True)  # Ensures explicit Hs

# # Save SDS molecule with hydrogens only on carbon atoms
# Chem.MolToMolFile(mol, "ligands/final_sds.sdf")


## CHECK CHARGE ON MOLECULE

# Load molecule from SDF
mol = Chem.SDMolSupplier("ligands/final_sds.sdf", removeHs=False)[0]

# Get the total formal charge
total_charge = Chem.GetFormalCharge(mol)

print(f"Total charge on the molecule: {total_charge}")