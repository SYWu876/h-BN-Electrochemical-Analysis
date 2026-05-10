# XPS reproducibility notes

This module reproduces the XPS analysis for the processed h-BN sample.

## Raw workbook format
The source file is supplied as legacy `.xls`. The helper functions try, in order:

1. direct `pandas.read_excel(...)`
2. conversion to `.xlsx` using LibreOffice in headless mode
3. reading the converted `.xlsx`

## Adopted 13-peak definition
The corrected ML analysis uses the peak centers and area fractions adopted in Table X:

- B1s: B-def, B–N, B-ox
- N1s: N-def, N–B, N-ox
- C1s: C–C, C–O, C=O, O–C=O
- O1s: O-lat, O-ads, H2O/O

## Descriptor set
For each fitted component, the ML feature vector is:

- Area fraction A_i
- pseudo-Voigt mixing parameter eta_i
- shift relative to dominant peak DeltaE_i
- FWHM Gamma_i

These descriptors are standardized before PCA and hierarchical clustering.
