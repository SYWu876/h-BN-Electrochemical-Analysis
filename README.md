# Structure-Guided Emergent Interfacial Electrochemistry in Liquid-Phase-Processed h-BN

This repository packages the data products, figure exports, and analysis scripts accompanying the manuscript **"Structure-Guided Emergent Interfacial Electrochemistry in Liquid-Phase-Processed h-BN"**.

The study treats liquid-phase-processed h-BN as a **structure-defined reference platform** for understanding how defect heterogeneity, partial restacking, and near-surface perturbation can generate a measurable but distinctly non-ideal interfacial electrochemical response in a nominally insulating layered solid. The manuscript links structural characterization (HRTEM/FFT, Raman, XPS) to electrochemical analysis across CV, GCD, and EIS, and further compares classical, continuous, and discrete quantum-assisted inference routes for the EIS model.

## Repository scope

This package is intended as a **GitHub-ready companion archive** for manuscript review, sharing, and later public release. It contains:

- raw experimental files for **TEM, Raman, XPS, CV, GCD, and EIS**
- processed CSV tables used to support the electrochemical analysis
- exported **main-text figures (Figures 1-8)** and **Supporting Information figures (Figures S1-S13)**
- lightweight Python scripts for rebuilding the included **EIS classical-fit and quantum-branch CSV outputs**
- a regenerated machine-readable file manifest and checksums for archive integrity

## Repository layout

```text
hBN_GitHub_Package_final/
├── README.md
├── CITATION.cff
├── requirements.txt
├── .gitignore
├── data/
│   ├── raw/
│   │   ├── CV/
│   │   ├── GCD/
│   │   ├── EIS/
│   │   ├── Raman/
│   │   ├── TEM/
│   │   └── XPS/
│   └── processed/
│       ├── CV/
│       ├── CV_segmentation/
│       ├── GCD/
│       └── EIS/
│           ├── classical_fit/
│           ├── quantum_branches/
│           └── qaoa_landscapes/
├── figures/
│   ├── main_text/
│   └── SI/
├── scripts/
│   ├── 01_classical_eis_fit.py
│   ├── 02_quantum_branch_comparison.py
│   └── 03_surrogate_qaoa_landscape.py
└── docs/
    ├── manifest.csv
    ├── checksums.sha256
    └── Note_S6_GitHub_v2.md
```

## What is in each folder?

### `data/raw/`
Original or source-level files used as inputs for analysis.

- `CV/`: scan-rate-resolved cyclic voltammetry tables and Dunn-partition source tables
- `GCD/`: full galvanostatic charge-discharge dataset and rate-specific preprocessing diagnostics
- `EIS/`: raw impedance spectrum used for the classical and quantum-branch comparison
- `Raman/`: raw Raman text export
- `TEM/`: source TEM image plus ROI-level enhancement, FFT, lattice-map, and maxima-overlay files
- `XPS/`: raw XPS spreadsheet export

### `data/processed/`
CSV tables prepared for plotting, fitting summaries, and manuscript-linked interpretation.

- `CV/`: descriptor tables, log(i)-log(v) outputs, and Dunn-separation products
- `CV_segmentation/`: Ising-type local-field labels, transition points, and Q-KPCA outputs
- `GCD/`: bounded-fit summaries, residual tables, and bootstrap outputs
- `EIS/classical_fit/`: final compact-circuit fit, fitted parameters, and residuals
- `EIS/quantum_branches/`: branch parameters, metrics, overlays, and deviations
- `EIS/qaoa_landscapes/`: surrogate slices and coarse/refined QAOA landscape tables

### `figures/`
PNG exports of manuscript figures.

- `main_text/`: Figures 1-8
- `SI/`: Figures S1-S13

### `scripts/`
Minimal Python scripts included for rebuilding the EIS data products distributed in `data/processed/EIS/`.

## Quick start

### 1. Create a Python environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Rebuild the EIS outputs

Run from the repository root:

```bash
python scripts/01_classical_eis_fit.py
python scripts/02_quantum_branch_comparison.py
python scripts/03_surrogate_qaoa_landscape.py
```

The scripts write CSV outputs into the corresponding folders under `data/processed/EIS/`.

## Manuscript-level analysis map

The study is organized around three electrochemical branches:

- **CV**: peak-current scaling, power-law `b`-value analysis, Dunn-type current separation, Ising-type kinetic segmentation, and Q-KPCA embedding
- **GCD**: preprocessing diagnostics, QAOA-based stable-window selection, bounded physics-informed fitting, and bootstrap/residual summaries
- **EIS**: compact classical anchor fit and comparison with continuous and discrete quantum-assisted inference under a shared complex-domain objective

The full file-to-analysis mapping is listed in `docs/manifest.csv`.

## Notes for public GitHub release

- macOS metadata files have been removed from this final package
- file integrity hashes are listed in `docs/checksums.sha256`
- before public release, you may additionally choose to add a project license and repository DOI/Zenodo record

## Suggested citation

If you use or redistribute this package, please cite the associated manuscript and reference this repository as its data-and-figure companion archive.

## Correspondence

**Sheng Yun Wu**  
Department of Physics, National Dong Hwa University, Hualien 97401, Taiwan  
Email: sywu@gms.ndhu.edu.tw
