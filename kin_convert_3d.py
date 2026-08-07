"""
Pipeline: Extract compounds → Generate 3D structures → Interactive HTML viewer

Steps:
  1. Read Library_export_11OCT.sdf
  2. Extract the 13 target compounds (handles KIN-XXXXXXX and KINXXXX formats)
  3. Generate / clean up 3D coordinates with RDKit (ETKDG + MMFF94 minimisation)
  4. Save RIGI_compounds_3D.sdf  (for docking input)
  5. Save RIGI_compounds_3D.html (interactive 3D viewer — open in any browser)

Usage:
    python generate_3D_and_html.py

Requirements:
    pip install rdkit
"""

import os, sys, io, base64, re

BASE_DIR  = "/Users/kanokpant.sriwong/Desktop/RIGI_docking"
INPUT_SDF = os.path.join(BASE_DIR, "Library_export_11OCT.sdf")
OUT_SDF   = os.path.join(BASE_DIR, "RIGI_compounds_3D.sdf")
OUT_HTML  = os.path.join(BASE_DIR, "RIGI_compounds_3D.html")

# Target compound IDs exactly as they appear in your list
RAW_TARGETS = [
    "KIN-0001148", "KIN-0121148", "KIN-0121408", "KIN-0121407",
    "KIN-0112000", "KIN-0110941", "KIN-0110865", "KIN-0000865",
    "KIN-0111501", "KIN-0001312", "KIN-0111312", "KIN-0121312",
    "KIN5059","KIN-0120865"
]

# ── ID normalisation ──────────────────────────────────────────────────────────
def normalise(s):
    """
    Strip formatting characters, the 'KIN' prefix (with or without dash),
    leading zeros, then uppercase.

    Examples
    --------
    "KIN-0120532"  →  "120532"
    "KIN5059"      →  "5059"
    "KIN0001148"   →  "1148"
    "0121148"      →  "121148"
    "5059"         →  "5059"
    """
    s = s.strip().upper().replace("-", "").replace("_", "").replace(" ", "")
    # Strip leading 'KIN' prefix if present
    if s.startswith("KIN"):
        s = s[3:]
    # Strip leading zeros so "0001148" == "1148"
    s = s.lstrip("0") or "0"
    return s

# Build normalised→original lookup for our targets
TARGETS = {normalise(t): t for t in RAW_TARGETS}

print("Target normalised keys:", TARGETS)  # debug

# ── SDF parsing helpers ───────────────────────────────────────────────────────
def iter_sdf_records(path):
    """Yield one complete SDF record (including $$$$) at a time."""
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        buf = []
        for line in fh:
            buf.append(line)
            if line.strip() == "$$$$":
                yield "".join(buf)
                buf = []
        if buf:
            yield "".join(buf)

def parse_props(record):
    """
    Parse all > <TAG> blocks from an SDF record.
    Handles multi-word tags like 'Molecule Name 2'.
    Returns dict of {tag: value}.
    """
    props = {}
    lines = record.splitlines()
    i = 0
    while i < len(lines):
        ln = lines[i]
        # Match lines like:  >  <Molecule Name 2>  or  > <CompoundID>
        m = re.match(r'>\s*<([^>]+)>', ln)
        if m:
            tag = m.group(1).strip()
            vals = []
            i += 1
            while i < len(lines) and lines[i].strip() not in ("", "$$$$"):
                vals.append(lines[i])
                i += 1
            props[tag] = "\n".join(vals).strip()
        else:
            i += 1
    return props

# Extended list of ID tags — including the ones seen in your SDF
ID_TAGS = [
    # Standard / common
    "Compound_ID", "CompoundID", "COMPOUND_ID",
    "Name", "NAME", "ID",
    "Catalog_ID", "Reg_No", "RegistryNumber", "molregno",
    "Catalog Number", "Sample_ID",
    # Your library-specific tags
    "Molecule Name 2",       # <── the tag that holds KIN-0121148 in your data
    "Molecule Name",
    "Molecule_Name",
    "External_Compound_ID",
    "External_ID",
    "ChemblID",
    "IDNUMBER",
    "Compound Name",
]

def get_ids(record):
    """
    Return all candidate ID strings from a record:
      - mol name (first line of the molfile block)
      - every property tag value
      - multi-value fields split by newline (some SDF writers put
        several synonyms in one tag)
    """
    ids = []

    # First line of record = molecule name in V2000/V3000 molfile
    first_line = record.splitlines()[0].strip()
    if first_line:
        ids.append(first_line)

    props = parse_props(record)

    # Debug: on first few records print all tags found
    # (comment out after confirming)
    # print("  tags:", list(props.keys()))

    for tag in ID_TAGS:
        v = props.get(tag, "").strip()
        if v:
            # Some fields have multiple IDs separated by newlines
            for part in v.splitlines():
                part = part.strip()
                if part:
                    ids.append(part)

    # Also scan ALL property values for anything that looks like a KIN ID
    for tag, v in props.items():
        for part in v.splitlines():
            part = part.strip()
            if re.match(r'^KIN[-_]?\d+', part, re.IGNORECASE):
                ids.append(part)

    return ids

def find_match(record):
    """Return the original target label if this record matches any target, else None."""
    for raw in get_ids(record):
        n = normalise(raw)
        if not n:
            continue
        # Exact match after normalisation
        if n in TARGETS:
            return TARGETS[n]
        # Substring match (handles extra prefixes / suffixes)
        for tn, tl in TARGETS.items():
            if tn and (tn in n or n in tn):
                return tl
    return None

# ── 3D generation ─────────────────────────────────────────────────────────────
def make_3d(mol, label):
    from rdkit.Chem import AllChem
    mol = AllChem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        print(f"    ⚠  ETKDG failed for {label}, trying ETDG fallback")
        result = AllChem.EmbedMolecule(mol, AllChem.ETDG())
    if result == -1:
        print(f"    ⚠  All embedding failed for {label}, trying distance geometry")
        AllChem.EmbedMolecule(mol)
    try:
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94")
    except Exception as e:
        print(f"    ⚠  MMFF minimisation skipped for {label}: {e}")
    return mol

# ── HTML builder ──────────────────────────────────────────────────────────────
def mol_to_sdf_string(mol):
    from rdkit import Chem
    buf = io.StringIO()
    w = Chem.SDWriter(buf)
    w.write(mol)
    w.close()
    return buf.getvalue()

def build_html(mols_3d):
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw

    print("Building interactive HTML viewer …")

    # 2D thumbnail images (base64 PNG)
    thumb_size = (300, 200)
    thumbnails = {}
    for label, mol in mols_3d:
        mol2d = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol2d)
        img = Draw.MolToImage(mol2d, size=thumb_size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        thumbnails[label] = base64.b64encode(buf.getvalue()).decode()

    # SDF blocks for 3Dmol.js
    sdf_blocks = {label: mol_to_sdf_string(mol) for label, mol in mols_3d}

    # Heavy-atom counts for info bar
    atom_counts = {label: Chem.RemoveHs(mol).GetNumAtoms() for label, mol in mols_3d}

    labels_js = str([l for l, _ in mols_3d])
    sdf_js    = "{\n" + ",\n".join(
        f'  {repr(l)}: {repr(sdf_blocks[l])}' for l, _ in mols_3d
    ) + "\n}"
    thumb_js  = "{\n" + ",\n".join(
        f'  {repr(l)}: "data:image/png;base64,{thumbnails[l]}"' for l, _ in mols_3d
    ) + "\n}"
    atoms_js  = "{" + ", ".join(
        f'{repr(l)}: {atom_counts[l]}' for l, _ in mols_3d
    ) + "}"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>RIGI Docking Compounds — 3D Preview</title>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.7.1/jquery.min.js"></script>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #0f1117; color: #e2e8f0; min-height: 100vh; }}

  header {{
    background: linear-gradient(135deg, #1e293b, #0f172a);
    padding: 18px 28px;
    border-bottom: 2px solid #334155;
    display: flex; align-items: center; gap: 16px;
  }}
  header h1 {{ font-size: 1.3rem; font-weight: 700; color: #38bdf8; letter-spacing: 0.4px; }}
  header span {{ font-size: 0.82rem; color: #94a3b8; }}
  .badge {{
    background: #0ea5e9; color: #fff; border-radius: 20px;
    padding: 3px 12px; font-size: 0.76rem; font-weight: 600; margin-left: auto;
  }}

  .layout {{ display: flex; height: calc(100vh - 66px); }}

  .sidebar {{
    width: 255px; min-width: 200px; background: #1e293b;
    border-right: 1px solid #334155; overflow-y: auto; padding: 10px 0;
  }}
  .sidebar-title {{
    font-size: 0.68rem; font-weight: 700; color: #64748b;
    text-transform: uppercase; letter-spacing: 1px; padding: 6px 14px 10px;
  }}
  .compound-card {{
    padding: 9px 12px; cursor: pointer; border-left: 3px solid transparent;
    transition: background 0.15s, border-color 0.15s;
    display: flex; align-items: center; gap: 9px;
  }}
  .compound-card:hover  {{ background: #263348; }}
  .compound-card.active {{ background: #1d3a52; border-left-color: #38bdf8; }}
  .compound-card img    {{ width: 58px; height: 38px; object-fit: contain;
                           background: #fff; border-radius: 4px; flex-shrink: 0; }}
  .card-label {{ font-size: 0.8rem; font-weight: 600; color: #cbd5e1; word-break: break-all; }}
  .card-idx   {{ font-size: 0.68rem; color: #64748b; }}

  .main {{ flex: 1; display: flex; flex-direction: column; min-width: 0; }}

  .toolbar {{
    background: #1e293b; border-bottom: 1px solid #334155;
    padding: 9px 18px; display: flex; align-items: center; gap: 8px; flex-wrap: wrap;
  }}
  .toolbar label {{ font-size: 0.78rem; color: #94a3b8; }}
  .toolbar select, .toolbar button {{
    background: #334155; color: #e2e8f0; border: 1px solid #475569;
    border-radius: 6px; padding: 4px 11px; font-size: 0.78rem; cursor: pointer;
  }}
  .toolbar button:hover {{ background: #0ea5e9; border-color: #0ea5e9; color: #fff; }}
  #compound-title {{ font-size: 0.95rem; font-weight: 700; color: #38bdf8; margin-right: auto; }}

  #viewer-container {{ flex: 1; position: relative; background: #000; }}
  #viewer {{ width: 100%; height: 100%; }}

  .info-bar {{
    background: #1e293b; border-top: 1px solid #334155;
    padding: 7px 18px; font-size: 0.76rem; color: #64748b;
    display: flex; gap: 22px; flex-wrap: wrap;
  }}
  .info-bar b {{ color: #94a3b8; }}
</style>
</head>
<body>

<header>
  <div>
    <h1>RIGI Docking — Compound 3D Preview</h1>
    <span>Pre-docking structure review · ETKDG v3 + MMFF94</span>
  </div>
  <div class="badge">{len(mols_3d)} compounds</div>
</header>

<div class="layout">
  <div class="sidebar">
    <div class="sidebar-title">Compound List</div>
    <div id="compound-list"></div>
  </div>

  <div class="main">
    <div class="toolbar">
      <span id="compound-title">Select a compound</span>
      <label>Style</label>
      <select id="style-select">
        <option value="stick">Stick</option>
        <option value="sphere">Sphere</option>
        <option value="line">Line</option>
        <option value="cross">Cross</option>
      </select>
      <label>Color</label>
      <select id="color-select">
        <option value="element">By Element</option>
        <option value="spectrum">Spectrum</option>
        <option value="chain">By Chain</option>
      </select>
      <button id="spin-btn">⟳ Spin ON</button>
      <button id="bg-btn">☀ Light BG</button>
      <button id="surface-btn">◈ Surface</button>
    </div>

    <div id="viewer-container">
      <div id="viewer"></div>
    </div>

    <div class="info-bar">
      <span><b>Compound:</b> <span id="info-id">—</span></span>
      <span><b>Heavy atoms:</b> <span id="info-atoms">—</span></span>
      <span><b>3D method:</b> ETKDG v3 + MMFF94 (RDKit)</span>
    </div>
  </div>
</div>

<script>
const LABELS     = {labels_js};
const SDF_DATA   = {sdf_js};
const THUMBS     = {thumb_js};
const ATOM_COUNTS = {atoms_js};

let viewer, spinning = false, surfaceShown = false, surfaceId = null, lightBg = false;
let currentLabel = null;

function initViewer() {{
  viewer = $3Dmol.createViewer($("#viewer"), {{
    backgroundColor: "black",
    antialias: true,
  }});
}}

function applyStyle() {{
  const style = document.getElementById("style-select").value;
  const color = document.getElementById("color-select").value;
  viewer.setStyle({{}}, {{}});
  viewer.setStyle({{}}, {{ [style]: {{ colorscheme: color === "element" ? "Jmol" : color }} }});
  viewer.render();
}}

function loadCompound(label) {{
  currentLabel = label;
  viewer.clear();
  surfaceShown = false; surfaceId = null;

  viewer.addModel(SDF_DATA[label], "sdf");
  applyStyle();
  viewer.zoomTo();
  viewer.render();

  document.getElementById("compound-title").textContent = label;
  document.getElementById("info-id").textContent = label;
  document.getElementById("info-atoms").textContent = ATOM_COUNTS[label] ?? "?";

  document.querySelectorAll(".compound-card").forEach(el =>
    el.classList.toggle("active", el.dataset.label === label)
  );
}}

function buildSidebar() {{
  const list = document.getElementById("compound-list");
  LABELS.forEach((label, i) => {{
    const card = document.createElement("div");
    card.className = "compound-card";
    card.dataset.label = label;
    card.innerHTML = `
      <img src="${{THUMBS[label]}}" alt="${{label}}">
      <div>
        <div class="card-label">${{label}}</div>
        <div class="card-idx">#${{i + 1}}</div>
      </div>`;
    card.addEventListener("click", () => loadCompound(label));
    list.appendChild(card);
  }});
}}

document.getElementById("style-select").addEventListener("change", applyStyle);
document.getElementById("color-select").addEventListener("change", applyStyle);

document.getElementById("spin-btn").addEventListener("click", function () {{
  spinning = !spinning;
  this.textContent = spinning ? "⟳ Spin OFF" : "⟳ Spin ON";
  viewer.spin(spinning);
}});

document.getElementById("bg-btn").addEventListener("click", function () {{
  lightBg = !lightBg;
  viewer.setBackgroundColor(lightBg ? "white" : "black");
  viewer.render();
  this.textContent = lightBg ? "🌙 Dark BG" : "☀ Light BG";
}});

document.getElementById("surface-btn").addEventListener("click", function () {{
  if (!currentLabel) return;
  if (surfaceShown) {{
    viewer.removeSurface(surfaceId);
    surfaceShown = false;
    this.textContent = "◈ Surface";
  }} else {{
    surfaceId = viewer.addSurface($3Dmol.SurfaceType.VDW, {{
      opacity: 0.5, colorscheme: "whiteCarbon"
    }});
    surfaceShown = true;
    this.textContent = "◈ Hide Surface";
  }}
  viewer.render();
}});

initViewer();
buildSidebar();
if (LABELS.length > 0) loadCompound(LABELS[0]);
</script>
</body>
</html>"""

    with open(OUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"HTML viewer saved → {OUT_HTML}")
    print("Open the HTML file in any browser to review compounds.")


# ── Diagnostic: scan first N records and print all tags found ─────────────────
def diagnose(path, n=5):
    """Print property tags and candidate IDs from first N records."""
    print(f"\n{'='*60}")
    print(f"DIAGNOSTIC — first {n} records of {path}")
    print(f"{'='*60}")
    for i, record in enumerate(iter_sdf_records(path)):
        if i >= n:
            break
        first_line = record.splitlines()[0].strip()
        props = parse_props(record)
        ids = get_ids(record)
        print(f"\n--- Record {i+1} ---")
        print(f"  Molfile name (line 1): {repr(first_line)}")
        print(f"  Property tags: {list(props.keys())}")
        print(f"  All candidate IDs: {ids}")
        # Check if any match
        match = find_match(record)
        if match:
            print(f"  ✓ MATCH → {match}")
    print(f"{'='*60}\n")


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    from rdkit import Chem

    if not os.path.isfile(INPUT_SDF):
        sys.exit(f"ERROR: Cannot find {INPUT_SDF}")

    # Run diagnostic on first 5 records so you can see the actual tag names
    diagnose(INPUT_SDF, n=5)

    matched, total = {}, 0
    print(f"Scanning: {INPUT_SDF}\n")
    for record in iter_sdf_records(INPUT_SDF):
        total += 1
        label = find_match(record)
        if label and label not in matched:
            matched[label] = record
            print(f"  ✓  [{label}]  record #{total}")

    missing = [t for t in RAW_TARGETS if t not in matched]
    print(f"\n{'─'*60}")
    print(f"Scanned: {total}  |  Found: {len(matched)}/13  |  Missing: {missing or 'none'}")
    print(f"{'─'*60}\n")

    if not matched:
        print("ERROR: No compounds found. Check the DIAGNOSTIC output above to see")
        print("what property tags are actually in your SDF, then add them to ID_TAGS.")
        sys.exit(1)

    # Generate 3D for each compound
    mols_3d = []
    writer = Chem.SDWriter(OUT_SDF)

    for label in RAW_TARGETS:
        if label not in matched:
            continue
        record = matched[label]
        mol = Chem.MolFromMolBlock(record, removeHs=False, sanitize=True)
        if mol is None:
            # V3000 molfiles sometimes need explicit flag
            mol = Chem.MolFromMolBlock(record, removeHs=False, sanitize=False)
            if mol is None:
                print(f"  ✗  RDKit could not parse {label} — skipping")
                continue
            try:
                Chem.SanitizeMol(mol)
            except Exception as e:
                print(f"  ✗  Sanitization failed for {label}: {e} — skipping")
                continue

        mol = Chem.RemoveHs(mol)   # clean slate before re-adding H

        conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None
        if conf and conf.Is3D():
            print(f"  →  {label}  already 3D — keeping + minimising")
        else:
            print(f"  →  {label}  generating 3D coords …")

        mol3d = make_3d(mol, label)
        mol3d.SetProp("_Name", label)
        mol3d.SetProp("Compound_ID", label)
        writer.write(mol3d)
        mols_3d.append((label, mol3d))
        print(f"       done  ({mol3d.GetNumAtoms()} atoms incl. H)")

    writer.close()
    print(f"\n3D SDF saved → {OUT_SDF}\n")

    if mols_3d:
        build_html(mols_3d)
    else:
        print("No molecules successfully processed — HTML not generated.")

if __name__ == "__main__":
    main()