# Virtual screening of small-molecule RIG-I agonists

## Overview

A computational pipeline for preparing small-molecule libraries, identifying potential binding pockets in RIG-I, and prioritizing candidate compounds through structure-based virtual screening.

## Workflow

1. **Compound preparation:** Extract candidate molecules from SDF libraries, generate 3D conformations using RDKit (ETKDGv3), and minimize structures with MMFF94.
2. **Binding-pocket identification:** Identify potential druggable pockets using DoGSiteScorer.
3. **Molecular docking:** Dock 14 candidate compounds to full-length RIG-I and its C-terminal domain (CTD) using AutoDock Vina.
4. **Candidate prioritization:** Compare predicted binding affinities across targets and docking modes to rank compounds for further investigation.

## Tools

Python, RDKit, Open Babel, DoGSiteScorer, AutoDock Vina.

## Key Outputs

* Prepared 3D compound structures
* Binding-pocket analysis results for full lengths RIGI
  The highest-scoring pocket overlaps the conserved nucleotide-binding site in the RIG-I helicase domain. This region contains Motif I (GCGKT; residues 267–271), including Lys270, and Motif II (DECH; residues 372–375), which are key components of the ATPase machinery (Jiang et al., 2011).
  <img width="1401" height="898" alt="Screenshot 2026-09-30 at 10 48 56 PM" src="https://github.com/user-attachments/assets/fb1058db-d7f0-4890-acab-fc5bb7c6d865" />
  
* Docking scores and candidate rankings
  
  <img width="603" height="314" alt="Screenshot 2026-09-29 at 10 42 31 PM" src="https://github.com/user-attachments/assets/09509390-2650-4d45-ba0c-e1148d2aa826" />

* Binding-pocket analysis results for CTD domain of RIGI
  Predicted binding pockets on the RIG-I C-Terminal domain. Pockets were identified using DoGSiteScorer (Volkamer et al., 2012) on the RIG-I CTD structure (PDB 2QFB). Pocket P0 had a volume of 417.86 Å³ and a drug score of 0.62, whereas pocket P1 had a volume of 537.57 Å³ and a drug score of 0.50. (A) Pocket P0 is lined by residues including His830, Phe853, Lys858, Lys861, Lys888, and Lys907 and overlaps the positively charged groove proposed as the 5′-triphosphate-binding site (Cui et al., 2008), P0 scores therefore indicate possible competition with 5′ppp RNA. (B) Pocket P1 is lined by residues 803–806, 820–824, 831, and 911–916; no specific functional role has been assigned to this site.

* Docking scores and candidate rankings
  <img width="932" height="490" alt="Screenshot 2026-09-30 at 9 21 18 PM" src="https://github.com/user-attachments/assets/4f740622-c96f-40a9-bdb1-dfd0f2e628a3" />

## Limitations

Docking scores and predicted binding poses are computational predictions and do not establish RIG-I agonist activity. Experimental validation is required to determine biological activity.
