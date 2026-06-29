# Note S7. Machine-learning-assisted Raman line-shape analysis of the processed h-BN spectrum

This note documents the Raman workflow used to support Figure NS14. The analysis includes asymmetric least-squares baseline correction, pseudo-Voigt fitting of the dominant h-BN band near 1367 cm^-1, and unsupervised k-means segmentation of local spectral descriptors [I, dI/dω, d²I/dω²].

## Equations

### ALS baseline objective
min_z Σ_i w_i (y_i - z_i)^2 + λ Σ_i (Δ² z_i)^2
$$\min_z \sum_i w_i (y_i - z_i)^2 + \lambda \sum_i (\Delta^2 z_i)^2$$

### Corrected intensity
I_corr(ω_i) = y_i - z_i
$$I_{\text{corr}}(\omega_i) = y_i - z_i$$

### Pseudo-Voigt profile
I_PV(ω) = A [η L(ω; ω0, γ) + (1-η) G(ω; ω0, σ)]
$$I_{\text{PV}}(\omega) = A [\eta L(\omega; \omega_0, \gamma) + (1-\eta) G(\omega; \omega_0, \sigma)]$$

### Lorentzian component
L(ω; ω0, γ) = γ² / [(ω - ω0)² + γ²]
$$L(\omega; \omega_0, \gamma) = \frac{\gamma^2}{(\omega - \omega_0)^2 + \gamma^2}$$

### Gaussian component
G(ω; ω0, σ) = exp(-(ω - ω0)² / (2σ²))
$$G(\omega; \omega_0, \sigma) = \exp\left(-\frac{(\omega - \omega_0)^2}{2\sigma^2}\right)$$

### FWHM
FWHM = ω_R - ω_L, where I_PV(ω_L) = I_PV(ω_R) = 0.5 I_PV,max
$$\text{FWHM} = \omega_R - \omega_L, \quad \text{where } I_{\text{PV}}(\omega_L) = I_{\text{PV}}(\omega_R) = 0.5 I_{\text{PV,max}}$$

### Local descriptor vector
f_i = [I_corr(ω_i), dI/dω|_i, d²I/dω²|_i]
$$\mathbf{f}_i = \left[I_{\text{corr}}(\omega_i), \left.\frac{dI}{d\omega}\right|_i, \left.\frac{d^2I}{d\omega^2}\right|_i\right]$$

### Standardization
f_tilde = (f - μ) / s
$$\tilde{\mathbf{f}} = \frac{\mathbf{f} - \boldsymbol{\mu}}{\mathbf{s}}$$

### K-means objective
min Σ_k Σ_(x_i in C_k) ||x_i - μ_k||²
$$\min \sum_k \sum_{\mathbf{x}_i \in C_k} \|\mathbf{x}_i - \boldsymbol{\mu}_k\|^2$$

## Outputs

- Figure NS14a: raw Raman spectrum with ALS baseline and pseudo-Voigt fit
- Figure NS14b: baseline-corrected peak with fitted center and FWHM guide
- Figure NS14c: k-means segmentation of the corrected Raman band
- Figure NS14d: schematic link between Raman heterogeneity and CV/GCD/EIS interpretation
