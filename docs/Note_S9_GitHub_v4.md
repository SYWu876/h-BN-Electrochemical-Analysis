# Note S9. GitHub package organization, file inventory, and reproducibility scope

This package accompanies the manuscript **"Structure-Guided Emergent Interfacial Electrochemistry in Liquid-Phase-Processed h-BN"** and is organized as a GitHub-ready companion archive. It collects the raw measurement files, processed CSV outputs, figure exports, and lightweight analysis scripts needed to document the structure-guided CV, GCD, and EIS workflow used in the manuscript and Supporting Information.

## Included content

The archive contains four main groups of materials:

1. **Raw data** in `data/raw/` for TEM, Raman, XPS, CV, GCD, and EIS.
2. **Processed data** in `data/processed/` for CV descriptors, CV segmentation/Q-KPCA outputs, GCD bounded-fit summaries, and classical/quantum EIS tables.
3. **Figure exports** in `figures/main_text/` and `figures/SI/` covering Figures 1-8 and Figures S1-S13.
4. **Reusable scripts** in `scripts/` for rebuilding the distributed EIS fitting and landscape tables.

## Reproducibility scope

This package is designed primarily as a transparent numerical and graphical archive.

- The **CV** branch is represented by raw and processed CSV tables that support rate-dependent scaling, `b`-value analysis, Dunn separation, Ising-type segmentation, and Q-KPCA embedding.
- The **GCD** branch is represented by raw traces, preprocessing diagnostics, bounded-fit summaries, bootstrap outputs, and residual tables supporting the stable-window workflow.
- The **EIS** branch includes both the raw spectrum and executable scripts that regenerate the classical-fit, quantum-branch-comparison, surrogate-slice, and QAOA-landscape CSV products.

Accordingly, the present package is suitable for review sharing, manuscript archiving, and later migration to a public GitHub repository.

## File manifest

A regenerated manifest file is provided at `docs/manifest.csv`. Each row lists the package path, analysis family, file type, and manuscript-facing role of the corresponding artifact. Integrity hashes for all tracked files are provided in `docs/checksums.sha256`.

## Recommended next step for public release

For a public repository release, the present package can be used directly. The only optional additions would be a project license and a permanent repository DOI.
