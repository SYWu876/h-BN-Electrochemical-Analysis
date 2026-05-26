# h-BN TEM reproducibility package

This package contains the TEM-analysis part for the GitHub repository:

SYWu876/h-BN-Electrochemical-Analysis

## Included contents

- `scripts/tem/00_tem_patch_ensemble_analysis.py`
  Rebuilds the TEM descriptor workflow for the five-patch ensemble analysis.
- `data/raw/TEM/OneView 200kV 800kX 39972.jpg`
  Raw TEM image used for ROI cropping.
- `data/raw/TEM/roi_boxes_template.csv`
  Manuscript ROI crop boxes mapped onto the raw TEM image.
- `data/processed/TEM/`
  Manuscript reference descriptor tables and generated reconstruction figures.

## ROI selection provenance

The five manuscript ROIs were selected with the same logic used in the
interactive TEM analysis helper:

1. The TEM image is converted to grayscale and resized to a 512 x 512 working
   canvas while retaining the scale factor back to the original image.
2. A local texture map is computed as the sliding-window standard deviation,
   `sqrt(mean(I^2) - mean(I)^2)`, using a kernel of one quarter of the ROI
   window size.
3. Square candidate ROIs are scanned on a regular grid with a margin of half
   the ROI size and a step equal to the ROI size. The interactive helper used a
   32 x 32 pixel ROI on the 512 x 512 working canvas.
4. Candidate ROIs below the 40th percentile of texture score are removed. If
   too few candidates remain, the highest-texture windows are retained.
5. Each remaining candidate is analyzed by FFT after suppressing the central DC
   region, giving a dominant peak radius and orientation.
6. Candidates with similar spatial position, d-spacing, and FFT orientation are
   merged into lattice domains. Representative ROIs are then chosen by taking
   the highest-texture ROI from each primary domain and filling to at most five
   ROIs by texture score.

The repository stores the final manuscript ROI boxes in original-image pixel
coordinates in `roi_boxes_template.csv`. The analysis script uses these stored
coordinates as the reproducible archive input and writes an overlay image so the
selected regions can be checked against the raw TEM micrograph.

## Raw TEM image filename convention

This package uses the raw TEM image file:

`OneView 200kV 800kX 39972.jpg`

Expected location:
`data/raw/TEM/OneView 200kV 800kX 39972.jpg`

## Supported workflows

### 1. Run from pre-cropped ROI images
Place ROI images in:
`data/raw/TEM/rois/`

Example:
```bash
python scripts/tem/00_tem_patch_ensemble_analysis.py \
    --repo-root . \
    --roi-dir data/raw/TEM/rois
```

### 2. Run from the raw TEM image + ROI box file
The repository includes the manuscript ROI boxes in:
`data/raw/TEM/roi_boxes_template.csv`

Example:
```bash
python scripts/tem/00_tem_patch_ensemble_analysis.py \
    --repo-root . \
    --source-image "data/raw/TEM/OneView 200kV 800kX 39972.jpg" \
    --roi-csv data/raw/TEM/roi_boxes_template.csv
```

The default run uses the committed manuscript descriptor table,
`data/processed/TEM/TEM_descriptors.csv`, for the order-disorder map,
geometric patch weights, and descriptor heatmap. It also writes
`TEM_roi_selection_overlay.png/.pdf` so the ROI boxes can be visually checked
against the raw source image. Use `--recompute-descriptors` only for
exploratory checks of the current ROI crops; those estimates are not intended
to replace the manuscript reference table.

## Included example outputs in processed/TEM

- `TEM_descriptors.csv`
- `TEM_descriptors_normalized.csv`
- `TEM_order_disorder_map.png`
- `TEM_geometric_patch_weights.png`
- `TEM_descriptor_heatmap.png`

## Notes

- Add `matplotlib` to `requirements.txt` if it is not already listed.
- The workflow is intended for comparative structural-order analysis and five-patch ensemble representation, not definitive atomic-species identification.
