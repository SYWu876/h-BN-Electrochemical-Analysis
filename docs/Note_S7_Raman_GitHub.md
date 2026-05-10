# Note S7. Machine-learning-assisted Raman line-shape analysis of the processed h-BN spectrum

This note documents the Raman workflow used to support Figure NS14. The analysis includes asymmetric least-squares baseline correction, pseudo-Voigt fitting of the dominant h-BN band near 1367 cm^-1, and unsupervised k-means segmentation of local spectral descriptors [I, dI/dω, d²I/dω²].

## Equations

### ALS baseline objective
min_z Σ_i w_i (y_i - z_i)^2 + λ Σ_i (Δ² z_i)^2

### Corrected intensity
I_corr(ω_i) = y_i - z_i

### Pseudo-Voigt profile
I_PV(ω) = A [η L(ω; ω0, γ) + (1-η) G(ω; ω0, σ)]

### Lorentzian component
L(ω; ω0, γ) = γ² / [(ω - ω0)² + γ²]

### Gaussian component
G(ω; ω0, σ) = exp(-(ω - ω0)² / (2σ²))

### FWHM
FWHM = ω_R - ω_L, where I_PV(ω_L) = I_PV(ω_R) = 0.5 I_PV,max

### Local descriptor vector
f_i = [I_corr(ω_i), dI/dω|_i, d²I/dω²|_i]

### Standardization
f_tilde = (f - μ) / s

### K-means objective
min Σ_k Σ_(x_i in C_k) ||x_i - μ_k||²

## Outputs

- Figure NS14a: raw Raman spectrum with ALS baseline and pseudo-Voigt fit
- Figure NS14b: baseline-corrected peak with fitted center and FWHM guide
- Figure NS14c: k-means segmentation of the corrected Raman band
- Figure NS14d: schematic link between Raman heterogeneity and CV/GCD/EIS interpretation
