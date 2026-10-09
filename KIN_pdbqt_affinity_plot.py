
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from pathlib import Path
import re

from pathlib import Path

# Main folder
project_dir = Path(
    "/Users/kanokpant.sriwong/Desktop/RIGI_docking"
)

# Input SDF file
input_sdf = project_dir / "RIGI_compounds_3D.sdf"

# Output folder for prepared ligands
output_dir = project_dir / "kin_pdbqt"
output_dir.mkdir(parents=True, exist_ok=True)

print("Input file exists:", input_sdf.exists())
print("Output directory:", output_dir)

# Read compounds
supplier = Chem.SDMolSupplier(str(input_sdf), removeHs=False)
molecules = [mol for mol in supplier if mol is not None]

print(f"Found {len(molecules)} compounds")

preparator = MoleculePreparation()
results = []

for i, mol in enumerate(molecules, start=1):
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"compound_{i}"
    safe_name = re.sub(r"[^A-Za-z0-9_.-]", "_", name)

    try:
        # Add hydrogens, retaining existing 3D coordinates
        mol = Chem.AddHs(mol, addCoords=True)

        # Optimize 3D structure if MMFF parameters are available
        if AllChem.MMFFHasAllMoleculeParams(mol):
            status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            if status == 1:
                print(f"Warning: Optimization not converged: {name}")
        else:
            print(f"Warning: MMFF parameters unavailable: {name}")

        # Prepare PDBQT
        setups = preparator.prepare(mol)
        if not setups:
            raise ValueError("No ligand setup generated")

        pdbqt, is_ok, error = PDBQTWriterLegacy.write_string(setups[0])
        if not is_ok:
            raise ValueError(str(error))

        # Save individual compound
        output_file = output_dir / f"{safe_name}.pdbqt"
        output_file.write_text(pdbqt)

        results.append((name, "Success"))
        print(f"Prepared: {name}")

    except Exception as e:
        results.append((name, f"Failed: {e}"))
        print(f"ERROR: {name}: {e}")

# Summary
print("\n--- Preparation Summary ---")
for name, status in results:
    print(f"{name}: {status}")

successful = sum(status == "Success" for _, status in results)
print(f"\nPrepared {successful}/{len(molecules)} compounds")
print(f"Output directory: {output_dir}")

import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10,
    "axes.linewidth": 1.0,
    "pdf.fonttype": 42,
})

df = pd.read_csv("/Users/kanokpant.sriwong/Downloads/affinity_summary_P0.csv")
df = df.sort_values("P0_Best_Affinity").reset_index(drop=True)  # strongest first

fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.bar(df["Compound"], df["P0_Best_Affinity"], width=0.65,
       color="#7f7f7f", edgecolor="black", linewidth=0.8)
bars = ax.bar(df["Compound"], df["P0_Best_Affinity"], width=0.65,
              color="#7f7f7f", edgecolor="black", linewidth=0.8)

for bar, val in zip(bars, df["P0_Best_Affinity"]):
    ax.text(bar.get_x() + bar.get_width() / 2, val - 0.15, f"{val:.2f}",
            ha="center", va="top", fontsize=8)

ax.set_ylim(-10, 0)
ax.set_yticks(range(-10, 1, 2))
ax.set_ylabel("Binding affinity (kcal/mol)")
ax.set_xlabel("Compound name")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.xaxis.set_ticks_position("bottom")
ax.tick_params(axis="x", pad=4)
# keep labels below the plot even though bars hang down from zero
ax.spines["bottom"].set_position(("data", -10))
ax.spines["left"].set_bounds(-10, 0)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

plt.show()
