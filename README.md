# 3D Compound structure generator for molecular docking

This project automates the preparation of small-molecule structures for molecular docking by extracting selected compounds from a chemical library, generating optimized three-dimensional conformations, and creating an interactive web-based molecular viewer. The resulting structures can be used as input for docking studies while providing a convenient platform for visual quality assessment.

---

## Methodology

### 1. Compound Selection

The script reads an SDF chemical library and automatically identifies user-defined compounds based on their compound identifiers (KIN IDs). Multiple identifier formats (e.g., `KIN-0121148`, `KIN0121148`, and `KIN0001148`) are normalized to ensure accurate compound retrieval regardless of naming convention.

### 2. Three-Dimensional Structure Generation

Each selected molecule is processed using **RDKit**, an open-source cheminformatics toolkit.

The workflow consists of:

- Removal of existing hydrogen atoms
- Addition of explicit hydrogens
- Three-dimensional conformer generation using the **ETKDG v3 (Experimental-Torsion Knowledge Distance Geometry)** algorithm
- Geometry optimization using the **MMFF94 (Merck Molecular Force Field)** to obtain a low-energy molecular conformation suitable for docking studies

If ETKDG fails to generate a conformer, the script automatically falls back to alternative distance geometry embedding methods to maximize successful structure generation.

### 3. Output Generation

The workflow produces two outputs:

- **`RIGI_compounds_3D.sdf`** — optimized 3D molecular structures suitable for molecular docking software.
- **`RIGI_compounds_3D.html`** — an interactive browser-based viewer for inspecting each compound prior to docking.

---

## Interactive 3D Viewer

The generated HTML viewer enables rapid visual inspection of every optimized ligand before molecular docking.

### Features

- Browse compounds individually
- Interactive 3D rotation, zooming, and panning
- Multiple rendering styles
  - Stick
  - Sphere
  - Line
  - Cross
- Van der Waals surface visualization
- Light/Dark background switching
- Automatic 2D structure thumbnails
- Heavy atom count display
- Browser-based interface
- 
**Open the viewer here:**

`RIGI_compounds_3D.html`

---

## Workflow

```text
Library_export_11OCT.sdf
           │
           ▼
Compound ID Matching
           │
           ▼
Target Compound Extraction
           │
           ▼
RDKit Processing
(Remove H → Add H)
           │
           ▼
ETKDG v3 Conformer Generation
           │
           ▼
MMFF94 Energy Minimization
           │
     ┌─────┴──────────┐
     ▼                ▼
3D Docking SDF   Interactive HTML Viewer
```

---

# Results

The workflow successfully extracted the selected RIG-I compounds from the screening library and generated optimized three-dimensional conformations using RDKit. 
Initial conformers were generated with the ETKDG v3 algorithm, which combines distance geometry with experimentally derived torsional preferences to produce chemically realistic structures. 
Each conformer was subsequently refined using MMFF94 energy minimization to obtain low-energy geometries appropriate for downstream molecular docking. 
The final structures were exported as a docking-ready SDF file and visualized through an interactive HTML interface.

---

## Methods and Tools

- Python
- RDKit
- ETKDG v3 conformer generation
- MMFF94 force-field optimization
- HTML
- JavaScript
- 3Dmol.js
