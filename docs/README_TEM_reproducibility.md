# h-BN TEM reproducibility package

This package contains the TEM-analysis part for the GitHub repository:

SYWu876/h-BN-Electrochemical-Analysis

## Included contents

- `scripts/tem/00_tem_patch_ensemble_analysis.py`
  Rebuilds the TEM descriptor workflow for the five-patch ensemble analysis.
- `data/raw/TEM/OneView 200kV 800kX 39972.jpg`
  Raw TEM image used for ROI cropping.
- `data/raw/TEM/roi_boxes_template.csv`
  Template for ROI cropping from the raw TEM image.
- `data/processed/TEM/`
  Included example processed descriptor tables and reconstruction figures.

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
Prepare ROI boxes in:
`data/raw/TEM/roi_boxes_template.csv`

Example:
```bash
python scripts/tem/00_tem_patch_ensemble_analysis.py \
    --repo-root . \
    --source-image "data/raw/TEM/OneView 200kV 800kX 39972.jpg" \
    --roi-csv data/raw/TEM/roi_boxes_template.csv
```

## Included example outputs in processed/TEM

- `TEM_descriptors.csv`
- `TEM_descriptors_normalized.csv`
- `TEM_order_disorder_map.png`
- `TEM_geometric_patch_weights.png`
- `TEM_descriptor_heatmap.png`

## Notes

- Add `matplotlib` to `requirements.txt` if it is not already listed.
- The workflow is intended for comparative structural-order analysis and five-patch ensemble representation, not definitive atomic-species identification.
