# Note S7. Machine-learning-assisted Raman line-shape analysis of the processed h-BN spectrum

This note documents the Raman workflow used to support Figure NS14. The analysis includes asymmetric least-squares baseline correction, pseudo-Voigt fitting of the dominant h-BN band near 1367 cm^-1, and unsupervised k-means segmentation of local spectral descriptors [I, dI/dω, d²I/dω²].

## Equations

### ALS baseline objective
$$\min_z \sum_i w_i (y_i - z_i)^2 + \lambda \sum_i (\Delta^2 z_i)^2$$

### Corrected intensity
$$I_{\text{corr}}(\omega_i) = y_i - z_i$$

### Pseudo-Voigt profile
$$I_{\text{PV}}(\omega) = A [\eta L(\omega; \omega_0, \gamma) + (1-\eta) G(\omega; \omega_0, \sigma)]$$

### Lorentzian component
$$L(\omega; \omega_0, \gamma) = \frac{\gamma^2}{(\omega - \omega_0)^2 + \gamma^2}$$

### Gaussian component
$$G(\omega; \omega_0, \sigma) = \exp\left(-\frac{(\omega - \omega_0)^2}{2\sigma^2}\right)$$

### FWHM
$$\text{FWHM} = \omega_R - \omega_L, \quad \text{where } I_{\text{PV}}(\omega_L) = I_{\text{PV}}(\omega_R) = 0.5 I_{\text{PV,max}}$$

### Local descriptor vector
$$\mathbf{f}_i = \left[I_{\text{corr}}(\omega_i), \left.\frac{dI}{d\omega}\right|_i, \left.\frac{d^2I}{d\omega^2}\right|_i\right]$$

### Standardization
$$\tilde{\mathbf{f}} = \frac{\mathbf{f} - \boldsymbol{\mu}}{\mathbf{s}}$$

### K-means objective
$$\min \sum_k \sum_{\mathbf{x}_i \in C_k} \|\mathbf{x}_i - \boldsymbol{\mu}_k\|^2$$

## Outputs

- Figure NS14a: raw Raman spectrum with ALS baseline and pseudo-Voigt fit
- Figure NS14b: baseline-corrected peak with fitted center and FWHM guide
- Figure NS14c: k-means segmentation of the corrected Raman band
- Figure NS14d: schematic link between Raman heterogeneity and CV/GCD/EIS interpretation
