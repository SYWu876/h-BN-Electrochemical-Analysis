# Structure-Guided Emergent Interfacial Electrochemistry in Liquid-Phase-Processed h-BN

This repository packages the data products and analysis scripts accompanying the manuscript **"Structure-Guided Emergent Interfacial Electrochemistry in Liquid-Phase-Processed h-BN"**.

The study treats liquid-phase-processed h-BN as a **structure-defined reference platform** for understanding how defect heterogeneity, partial restacking, and near-surface perturbation can generate a measurable but distinctly non-ideal interfacial electrochemical response in a nominally insulating layered solid. The manuscript links structural characterization (HRTEM/FFT, Raman, XPS) to electrochemical analysis across CV, GCD, and EIS, and further compares classical, continuous, and discrete quantum-assisted inference routes for the EIS model.

## Repository scope

This package is intended as a **GitHub-ready companion archive** for manuscript review, sharing, and later public release. It contains:

- raw experimental files for **TEM, Raman, XPS, CV, GCD, and EIS**
- processed CSV tables used to support structural, spectroscopic, and electrochemical analysis
- Python scripts for rebuilding selected **TEM, Raman, XPS, CV, GCD, and EIS** outputs
- repository metadata, license files, smoke tests, and CI configuration
- a regenerated machine-readable file manifest and checksums for archive integrity

## Repository layout

```text
hBN_GitHub_Package_final/
|-- README.md
|-- CITATION.cff
|-- LICENSE
|-- DATA_LICENSE.md
|-- NOTICE
|-- requirements.txt
|-- requirements-dev.txt
|-- data/
|   |-- raw/
|   |   |-- CV/
|   |   |-- GCD/
|   |   |-- EIS/
|   |   |-- Raman/
|   |   |-- TEM/
|   |   `-- XPS/
|   `-- processed/
|       |-- CV/
|       |-- CV_segmentation/
|       |-- GCD/
|       |-- EIS/
|       |-- Raman/
|       |-- TEM/
|       `-- XPS/
|-- scripts/
|   |-- cv/
|   |-- eis/
|   |   |-- 01_eis_classical_anchor_fit.py
|   |   |-- 02_eis_quantum_comparison_from_anchor.py
|   |   |-- 03_eis_surrogate_qaoa_landscape.py
|   |   `-- 04_eis_shared_objective_full_pipeline.py
|   |-- gcd/
|   |-- integrated/
|   |-- raman/
|   |-- tem/
|   `-- xps/
|-- tests/
`-- docs/
    |-- manifest.csv
    |-- checksums.sha256
    |-- Note_S7_Raman_GitHub.md
    |-- Note_S9_GitHub_v4.md
    |-- README_TEM_reproducibility.md
    `-- XPS_REPRODUCIBILITY_NOTES.md
```

## What is in each folder?

### `data/raw/`

Original or source-level files used as inputs for analysis.

- `CV/`: scan-rate-resolved cyclic voltammetry tables and Dunn-partition source tables
- `GCD/`: full galvanostatic charge-discharge dataset and rate-specific preprocessing diagnostics
- `EIS/`: raw impedance spectrum used for the classical and quantum-branch comparison
- `Raman/`: raw Raman text export
- `TEM/`: source TEM image plus ROI-level enhancement, FFT, lattice-map, and maxima-overlay files
- `XPS/`: raw XPS workbook exports, including an `.xlsx` fallback for reproducible script execution

### `data/processed/`

CSV tables prepared for plotting, fitting summaries, and manuscript-linked interpretation.

- `CV/`: descriptor tables, log(i)-log(v) outputs, Dunn-separation products, and figure-facing CV summaries
- `CV_segmentation/`: Ising-type local-field labels, transition points, and Q-KPCA outputs
- `GCD/`: bounded-fit summaries, residual tables, bootstrap outputs, diagnostics, and the final manuscript summary
- `integrated/`: conservative cross-domain descriptor tables for multi-modal evidence visualization
- `EIS/classical_fit/`: final compact-circuit fit, fitted parameters, and residuals
- `EIS/quantum_branches/`: branch parameters, metrics, overlays, and deviations
- `EIS/qaoa_landscapes/`: surrogate slices and coarse/refined QAOA landscape tables
- `Raman/`: baseline correction, peak fitting, and local spectral segmentation tables
- `TEM/`: patch-ensemble descriptor tables for the TEM reproducibility workflow
- `XPS/`: profile-fit summaries, adopted 13-peak table, and corrected descriptor matrix

### `scripts/`

Python scripts included for rebuilding selected analysis products. Domain-specific structural and electrochemical workflows are grouped under `scripts/cv/`, `scripts/eis/`, `scripts/gcd/`, `scripts/raman/`, `scripts/tem/`, and `scripts/xps/`. Cross-domain evidence tables are rebuilt from `scripts/integrated/`.

The first three EIS scripts in `scripts/eis/` are lightweight companion-archive rebuild helpers for selected processed EIS tables. In particular, `scripts/eis/02_eis_quantum_comparison_from_anchor.py` compares manuscript-linked continuous/discrete reference branch parameters against the classical EIS anchor, while `scripts/eis/03_eis_surrogate_qaoa_landscape.py` rebuilds the lightweight surrogate/QAOA CSV products.

For the full EIS shared-objective workflow, use `scripts/eis/04_eis_shared_objective_full_pipeline.py`. It starts from the raw EIS spectrum, performs the classical anchor fit, builds the local surrogate, evaluates the QAOA-compatible landscape, decodes the discrete branch, and writes the manuscript EIS surrogate slices, including `R1_Q1_slice.csv` and `Rs_alpha1_slice.csv`, to the selected output directory.

All EIS surrogate slices used in the manuscript can be regenerated from the full shared-objective pipeline.

## Quick start

### 1. Create a Python environment

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. Rebuild selected outputs

Run from the repository root:

```bash
python scripts/cv/00_run_all_cv_analysis.py
python scripts/gcd/00_run_all_gcd_analysis.py
python scripts/raman/01_raman_analysis_pipeline.py
python scripts/raman/02_plot_raman_ns14d_schematic.py
python scripts/tem/00_tem_patch_ensemble_analysis.py --repo-root . --source-image "data/raw/TEM/OneView 200kV 800kX 39972.jpg" --roi-csv data/raw/TEM/roi_boxes_template.csv
python scripts/xps/04_run_xps_pipeline.py
python scripts/integrated/00_build_cross_domain_evidence.py --repo-root .
python scripts/eis/01_eis_classical_anchor_fit.py
python scripts/eis/02_eis_quantum_comparison_from_anchor.py
python scripts/eis/03_eis_surrogate_qaoa_landscape.py
python scripts/eis/04_eis_shared_objective_full_pipeline.py --input data/raw/EIS/hbn_EIS_1.csv --output outputs/eis_shared_objective
```

The committed archive includes the manuscript-facing CSV tables. Generated figures and local rerun outputs are written to ignored output paths where possible, so they can be recreated locally without becoming part of the tracked archive.

## Manuscript-level analysis map

The study is organized around structural/spectroscopic inputs and three electrochemical branches:

- **TEM/Raman/XPS**: defect heterogeneity, lattice/patch descriptors, Raman line-shape analysis, and corrected XPS profile/descriptor tables
- **CV**: peak-current scaling, power-law `b`-value analysis, Dunn-type current separation, Ising-type kinetic segmentation, and Q-KPCA embedding
- **GCD**: preprocessing diagnostics, QAOA-based stable-window selection, bounded physics-informed fitting, final manuscript summary, and auto-scaled GCD summary regeneration
- **EIS**: compact classical anchor fit and comparison with continuous and discrete quantum-assisted inference under a shared complex-domain objective
- **Integrated evidence**: conservative descriptor matrix, heatmap table, and exploratory PCA projection assembled only from committed processed tables

The full file-to-analysis mapping is listed in `docs/manifest.csv`.

## Notes for public GitHub release

- macOS metadata files have been removed from this final package
- generated local figure/output folders are intentionally ignored
- file inventory and integrity hashes are listed in `docs/manifest.csv` and `docs/checksums.sha256`; these integrity files intentionally omit themselves
- repository licenses are included; before public release, you may additionally choose to add a repository DOI/Zenodo record

## License

Source code in `scripts/`, tests, and CI/configuration files are released under the MIT License; see `LICENSE`.

Data tables, image artifacts, and documentation are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0); see `DATA_LICENSE.md`.

The top-level `NOTICE` file records the path-level license scope for tools that need an explicit split-license summary.

## Suggested citation

If you use or redistribute this package, please cite the associated manuscript and reference this repository as its data-and-analysis companion archive.

## Correspondence

**Sheng Yun Wu**  
Department of Physics, National Dong Hwa University, Hualien 97401, Taiwan  
Email: sywu@gms.ndhu.edu.tw
