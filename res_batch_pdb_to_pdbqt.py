"""Batch-convert DoGSiteScorer pocket PDB files to rigid-receptor PDBQT for AutoDock Vina,
and write a boxes.csv with the Vina box (center + size) of every pocket.
No external packages needed. Vina ignores partial charges (written as 0.000);
the AutoDock atom type is what matters."""
import csv
import re
from pathlib import Path

IN_DIR = Path("/Users/kanokpant.sriwong/Desktop/RIGI_docking/CTD_pocket/residues")
OUT_DIR = IN_DIR.parent / "pdbqt"
PADDING = 8.0          # box size = 2 * pocket radius + PADDING (Angstrom)
OUT_DIR.mkdir(parents=True, exist_ok=True)

AROMATIC_C = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "ND1", "CD2", "CE1", "NE2"},
}
ACCEPTOR_N = {("HIS", "ND1"), ("HIS", "NE2")}

def ad_type(resname, atom, element):
    if element == "C":
        return "A" if atom in AROMATIC_C.get(resname, ()) else "C"
    if element == "N":
        return "NA" if (resname, atom) in ACCEPTOR_N else "N"
    if element == "O":
        return "OA"
    if element == "S":
        return "SA"
    return element.capitalize()

def convert(pdb_path):
    center = radius = None
    out = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith("HEADER") and "Geometric pocket center" in line:
            m = re.search(r"center at\s+(\S+)\s+(\S+)\s+(\S+)\s+with max radius\s+(\S+)", line)
            if m:
                center, radius = tuple(map(float, m.groups()[:3])), float(m.group(4))
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom, resname = line[12:16].strip(), line[17:20].strip()
        element = (line[76:78].strip() or re.sub(r"[^A-Za-z]", "", atom)[0]).upper()
        if element == "H":
            continue
        out.append(f"{line[:66]:<66}    0.000 {ad_type(resname, atom, element):<2}")
    return out, center, radius

rows = []
for pdb in sorted(IN_DIR.glob("*.pdb")):
    lines, center, radius = convert(pdb)
    (OUT_DIR / f"{pdb.stem}.pdbqt").write_text("\n".join(lines) + "\n")
    size = round(2 * radius + PADDING) if radius else ""
    rows.append([pdb.stem, len(lines), *(center or ("", "", "")), radius, size])
    print(f"{pdb.stem:<10} {len(lines):>3} atoms -> {pdb.stem}.pdbqt")

with open(OUT_DIR / "boxes.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["pocket", "n_atoms", "center_x", "center_y", "center_z", "radius", "box_size"])
    w.writerows(rows)
print(f"\nBox table -> {OUT_DIR / 'boxes.csv'}")
