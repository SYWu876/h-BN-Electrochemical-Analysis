# Note S9. GitHub package organization, file inventory, and reproducibility scope

This package accompanies the manuscript **"Structure-Guided Emergent Interfacial Electrochemistry in Liquid-Phase-Processed h-BN"** and is organized as a GitHub-ready companion archive. It collects the raw measurement files, processed CSV outputs, repository metadata, licensing files, smoke tests, CI configuration, and analysis scripts needed to document the TEM, Raman, XPS, CV, GCD, and EIS workflows used in the manuscript and Supporting Information.

## Included content

The archive contains four main groups of materials:

1. **Raw data** in `data/raw/` for TEM, Raman, XPS, CV, GCD, and EIS.
2. **Processed data** in `data/processed/` for TEM descriptors, Raman line-shape tables, XPS profile/descriptor tables, CV descriptors, CV segmentation/Q-KPCA outputs, GCD bounded-fit summaries, classical/quantum EIS tables, and integrated cross-domain descriptor summaries.
3. **Repository support files** including citation metadata, license files, smoke tests, CI configuration, file inventory, and integrity hashes.
4. **Reusable scripts** in `scripts/` for rebuilding selected TEM, Raman, XPS, CV, GCD, and EIS outputs.

## Reproducibility scope

This package is designed primarily as a transparent numerical and graphical archive.

- The **TEM/Raman/XPS** workflows are represented by raw inputs, processed descriptor/profile tables, and domain-specific reproducibility scripts under `scripts/tem/`, `scripts/raman/`, and `scripts/xps/`.
- The **CV** branch is represented by raw and processed CSV tables that support rate-dependent scaling, `b`-value analysis, Dunn separation, Ising-type segmentation, and Q-KPCA embedding.
- The **GCD** branch is represented by raw traces, preprocessing diagnostics, bounded-fit summaries, bootstrap outputs, residual tables, and the corrected final manuscript summary used for the auto-scaled GCD summary workflow.
- The **EIS** branch includes both the raw spectrum and executable scripts under `scripts/eis/` that regenerate selected classical-fit, quantum-branch-comparison, surrogate-slice, and QAOA-landscape CSV products.
- The **QC circuit schematic** branch provides an interactive notebook under `scripts/QC Circuit/` for regenerating the four GCD/EIS quantum-circuit figures used to document the hybrid workflows.
- The **integrated evidence** branch assembles a conservative descriptor matrix, heatmap table, and exploratory PCA projection from committed processed outputs only. These tables are intended for multi-modal visualization and do not imply statistical significance or causal inference.

The first three EIS scripts are lightweight companion-archive rebuild helpers. The continuous and discrete branch vectors used by `scripts/eis/02_eis_quantum_comparison_from_anchor.py` are manuscript-linked reference branch parameters for comparison against the classical anchor, and `scripts/eis/03_eis_surrogate_qaoa_landscape.py` rebuilds selected lightweight surrogate/QAOA CSV products.

The full shared-objective EIS workflow is provided by `scripts/eis/04_eis_shared_objective_full_pipeline.py`. It starts from the raw EIS spectrum, performs the classical anchor fit, builds the local surrogate, evaluates the QAOA-compatible landscape, decodes the discrete branch, and exports the manuscript EIS surrogate slices, including `R1_Q1_slice.csv` and `Rs_alpha1_slice.csv`, to the selected output directory. All EIS surrogate slices used in the manuscript can be regenerated from the full shared-objective pipeline.

Accordingly, the present package is suitable for review sharing, manuscript archiving, and later migration to a public GitHub repository.

## File manifest

A regenerated manifest file is provided at `docs/manifest.csv`. Each row lists the package path, analysis family, file type, and manuscript-facing role of the corresponding artifact. Integrity hashes for the payload files are provided in `docs/checksums.sha256`; the manifest and checksum files intentionally omit themselves.

## Recommended next step for public release

For a public repository release, the present package can be used directly. Repository licenses are included; the remaining optional addition would be a permanent repository DOI.
