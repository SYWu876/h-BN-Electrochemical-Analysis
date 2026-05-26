#!/usr/bin/env python3
"""
Reproducible TEM ROI analysis for the h-BN manuscript.

This script documents the ROI-level descriptor workflow described in Note S1:
- high-pass preprocessing
- conservative contrast normalization and mild denoising
- windowed FFT magnitude
- band-pass lattice reconstruction
- descriptor extraction: P_FFT, A_off, R_bp, Δq_FFT
- composite LOI
- regularized geometric patch weights
- ensemble centroid in the LOI–Δq_FFT plane

It is designed to fit a GitHub repository layout such as:
    data/raw/TEM/
    data/processed/TEM/

Supported input modes
---------------------
1) Cropped ROI images already prepared under `data/raw/TEM/rois/`
   Example files:
       ROI-1.png, ROI-2.png, ROI-3.png, ROI-4.png, ROI-5.png

2) A single source TEM image plus ROI coordinates in a CSV file.
   The CSV must contain columns:
       roi, x, y, width, height
   where (x, y) is the upper-left corner in pixels.

ROI selection provenance
------------------------
The manuscript ROI boxes stored in `data/raw/TEM/roi_boxes_template.csv` come
from the interactive TEM helper's texture-guided selection workflow: grayscale
512 x 512 working canvas, 32 x 32 candidate windows, local standard-deviation
texture scoring, 40th-percentile candidate filtering, FFT-based
d-spacing/orientation extraction, lattice-domain grouping, and selection of up
to five representative high-texture ROIs across the primary domains. This script
uses the stored manuscript ROI coordinates as the reproducible archive input.

Outputs
-------
- data/processed/TEM/TEM_descriptors.csv
- data/processed/TEM/TEM_descriptors_normalized.csv
- data/processed/TEM/TEM_roi_selection_overlay.png/.pdf
- data/processed/TEM/TEM_order_disorder_map.png/.pdf
- data/processed/TEM/TEM_geometric_patch_weights.png/.pdf
- data/processed/TEM/TEM_descriptor_heatmap.png/.pdf
- ROI-level PNG exports: processed image, FFT magnitude, lattice map, maxima overlay

By default, the manuscript descriptor table already committed under
`data/processed/TEM/TEM_descriptors.csv` is used for the Figure 1f-h style
summary panels. Use `--recompute-descriptors` only for exploratory checks of
the current ROI crops; those image-processing estimates are not a replacement
for the manuscript reference descriptor table.

Minimal example
---------------
From the repository root:
    python scripts/tem/00_tem_patch_ensemble_analysis.py \
        --repo-root . \
        --roi-dir data/raw/TEM/rois

Or with a source image + ROI table:
python scripts/tem/00_tem_patch_ensemble_analysis.py \
    --repo-root . \
    --source-image "data/raw/TEM/OneView 200kV 800kX 39972.jpg" \
    --roi-csv data/raw/TEM/roi_boxes_template.csv
"""

from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import fft as sp_fft
from scipy import ndimage as ndi
from scipy.signal import find_peaks, savgol_filter
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle


# ---------------------------
# Configuration and datatypes
# ---------------------------

ROI_NAME_ORDER = ["ROI-1", "ROI-2", "ROI-3", "ROI-4", "ROI-5"]
IMAGE_EXTENSIONS = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp", ".webp"}
DESCRIPTOR_COLUMNS = [
    "ROI",
    "P_FFT",
    "A_off",
    "R_bp",
    "Delta_q_FFT",
    "LOI_0_100",
    "geometric_patch_weight_wi",
]
NORMALIZED_DESCRIPTOR_COLUMNS = [
    "ROI",
    "P_FFT_norm",
    "A_off_norm",
    "R_bp_norm",
    "Delta_q_FFT_norm_raw",
    "Delta_q_FFT_ordering_tendency_inverted",
    "LOI_0_100",
    "geometric_patch_weight_wi",
]


@dataclass
class ROIResult:
    roi: str
    raw: np.ndarray
    processed: np.ndarray
    fft_log: np.ndarray
    bandpass_recon: np.ndarray
    maxima_rc: np.ndarray
    p_fft: float
    a_off: float
    r_bp: float
    dq_fft: float
    peak_radius: float
    radial_profile_q: np.ndarray
    radial_profile_val: np.ndarray


# ---------------------------
# Utility functions
# ---------------------------

def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def copy_file_if_different(source: Path, destination: Path) -> None:
    if source.resolve() == destination.resolve():
        return
    ensure_dir(destination.parent)
    destination.write_bytes(source.read_bytes())


def portable_source_path(path: Path, repo_root: Optional[Path] = None) -> str:
    if repo_root is not None:
        try:
            return path.resolve().relative_to(repo_root.resolve()).as_posix()
        except ValueError:
            pass
    if path.is_absolute():
        return path.name
    return path.as_posix()


def load_image(path: Path) -> np.ndarray:
    """Load an image as float grayscale in [0, 1]."""
    suffix = path.suffix.lower()
    arr = None

    # Pillow is often available transitively with matplotlib; use it if present.
    try:
        from PIL import Image
        with Image.open(path) as im:
            im = im.convert("L")
            arr = np.asarray(im, dtype=np.float64)
    except Exception:
        pass

    if arr is None:
        try:
            import imageio.v3 as iio
            arr = iio.imread(path)
            if arr.ndim == 3:
                arr = arr[..., :3].mean(axis=2)
            arr = arr.astype(np.float64)
        except Exception:
            pass

    if arr is None:
        arr = plt.imread(path)
        if arr.ndim == 3:
            arr = arr[..., :3].mean(axis=2)
        arr = arr.astype(np.float64)

    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    arr_min, arr_max = float(arr.min()), float(arr.max())
    if arr_max > arr_min:
        arr = (arr - arr_min) / (arr_max - arr_min)
    else:
        arr = np.zeros_like(arr, dtype=np.float64)
    return arr


def minmax(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    xmin = np.min(x)
    xmax = np.max(x)
    if xmax <= xmin:
        return np.zeros_like(x)
    return (x - xmin) / (xmax - xmin)


def robust_rescale(x: np.ndarray, lo: float = 1.0, hi: float = 99.0) -> np.ndarray:
    p_lo, p_hi = np.percentile(x, [lo, hi])
    if p_hi <= p_lo:
        return np.clip(x, 0, 1)
    y = (x - p_lo) / (p_hi - p_lo)
    return np.clip(y, 0.0, 1.0)


# ---------------------------
# Image preprocessing and FFT
# ---------------------------

def preprocess_roi(img: np.ndarray, sigma_bg: float = 12.0, sigma_denoise: float = 0.8) -> np.ndarray:
    """
    Conservative preprocessing consistent with Note S1:
    - subtract Gaussian background (high-pass)
    - robust local/global rescaling
    - mild Gaussian denoising
    """
    bg = ndi.gaussian_filter(img, sigma=sigma_bg, mode="reflect")
    high_pass = img - bg

    # Conservative contrast enhancement by robust rescaling.
    proc = robust_rescale(high_pass)
    proc = ndi.gaussian_filter(proc, sigma=sigma_denoise, mode="reflect")
    proc = robust_rescale(proc)
    return proc


def hann2d(shape: Tuple[int, int]) -> np.ndarray:
    wy = np.hanning(shape[0])
    wx = np.hanning(shape[1])
    return np.outer(wy, wx)


def fft_log_magnitude(img: np.ndarray) -> np.ndarray:
    window = hann2d(img.shape)
    centered = (img - np.mean(img)) * window
    F = sp_fft.fftshift(sp_fft.fft2(centered))
    return np.log1p(np.abs(F))


def radial_profile(image: np.ndarray, center: Optional[Tuple[float, float]] = None) -> Tuple[np.ndarray, np.ndarray]:
    h, w = image.shape
    cy = (h - 1) / 2.0 if center is None else center[0]
    cx = (w - 1) / 2.0 if center is None else center[1]

    y, x = np.indices(image.shape)
    r = np.sqrt((y - cy) ** 2 + (x - cx) ** 2)
    r_int = np.floor(r).astype(int)

    max_r = min(h, w) // 2
    valid = r_int < max_r
    r_int = r_int[valid]
    vals = image[valid]

    sums = np.bincount(r_int, weights=vals, minlength=max_r)
    counts = np.bincount(r_int, minlength=max_r)
    counts[counts == 0] = 1
    prof = sums / counts
    q = np.arange(max_r, dtype=np.float64)
    return q, prof


def estimate_peak_and_width(q: np.ndarray, prof: np.ndarray, low_cut_frac: float = 0.08) -> Tuple[float, float, float]:
    """
    Returns peak_radius, p_fft, dq_fft.
    dq_fft is the FWHM in pixel-radius units of the dominant annular FFT peak.
    """
    prof_s = savgol_filter(prof, window_length=min(len(prof) // 2 * 2 - 1, 11) if len(prof) >= 11 else max(3, len(prof) // 2 * 2 - 1), polyorder=2, mode="interp")
    low_cut = max(2, int(low_cut_frac * len(q)))

    q_sub = q[low_cut:]
    prof_sub = prof_s[low_cut:]
    if len(q_sub) < 5:
        peak_idx = int(np.argmax(prof_s))
    else:
        peaks, _ = find_peaks(prof_sub, prominence=max(np.std(prof_sub) * 0.1, 1e-8))
        if len(peaks) == 0:
            peak_idx = int(low_cut + np.argmax(prof_sub))
        else:
            # Choose the most intense non-central ring.
            rel = peaks[np.argmax(prof_sub[peaks])]
            peak_idx = int(low_cut + rel)

    peak_val = float(prof_s[peak_idx])

    # Background from the non-peak annular profile using a robust median.
    mask_bg = np.ones_like(prof_s, dtype=bool)
    left_bg = max(0, peak_idx - 5)
    right_bg = min(len(prof_s), peak_idx + 6)
    mask_bg[left_bg:right_bg] = False
    bg_vals = prof_s[mask_bg & (q >= low_cut)]
    bg_val = float(np.median(bg_vals)) if bg_vals.size else max(float(np.median(prof_s)), 1e-8)
    p_fft = peak_val / max(bg_val, 1e-8)

    # FWHM around the dominant peak.
    half = bg_val + 0.5 * (peak_val - bg_val)
    left = peak_idx
    while left > 0 and prof_s[left] >= half:
        left -= 1
    right = peak_idx
    while right < len(prof_s) - 1 and prof_s[right] >= half:
        right += 1

    # Linear interpolation for sub-bin width.
    def interp_x(i1: int, i2: int, target: float) -> float:
        x1, y1 = q[i1], prof_s[i1]
        x2, y2 = q[i2], prof_s[i2]
        if abs(y2 - y1) < 1e-12:
            return float(x1)
        return float(x1 + (target - y1) * (x2 - x1) / (y2 - y1))

    if left < peak_idx:
        x_left = interp_x(left, min(left + 1, len(q) - 1), half)
    else:
        x_left = float(q[left])
    if right > peak_idx:
        x_right = interp_x(max(right - 1, 0), right, half)
    else:
        x_right = float(q[right])

    dq_fft = max(x_right - x_left, 1e-6)
    return float(q[peak_idx]), float(p_fft), float(dq_fft)


def build_bandpass_mask(shape: Tuple[int, int], radius: float, sigma_r: float) -> np.ndarray:
    h, w = shape
    cy = (h - 1) / 2.0
    cx = (w - 1) / 2.0
    y, x = np.indices(shape)
    r = np.sqrt((y - cy) ** 2 + (x - cx) ** 2)
    mask = np.exp(-0.5 * ((r - radius) / max(sigma_r, 1e-6)) ** 2)
    return mask


def reconstruct_bandpass(img: np.ndarray, peak_radius: float, dq_fft: float) -> np.ndarray:
    window = hann2d(img.shape)
    centered = (img - np.mean(img)) * window
    F = sp_fft.fftshift(sp_fft.fft2(centered))
    sigma_r = max(dq_fft / 2.355, 1.5)  # convert approximate FWHM to sigma
    mask = build_bandpass_mask(img.shape, radius=peak_radius, sigma_r=sigma_r)
    recon = np.real(sp_fft.ifft2(sp_fft.ifftshift(F * mask)))
    recon = robust_rescale(recon)
    return recon


def autocorr_offcenter_amplitude(img: np.ndarray, exclusion_radius: int = 3) -> float:
    F = sp_fft.fft2(img - np.mean(img))
    ac = np.real(sp_fft.fftshift(sp_fft.ifft2(np.abs(F) ** 2)))
    ac = np.maximum(ac, 0.0)
    h, w = ac.shape
    cy = (h - 1) // 2
    cx = (w - 1) // 2
    y, x = np.indices(ac.shape)
    rr = np.sqrt((y - cy) ** 2 + (x - cx) ** 2)
    mask = rr > exclusion_radius
    denom = max(float(ac[cy, cx]), 1e-12)
    return float(np.max(ac[mask]) / denom)


def bandpass_rms(img: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(img))))


def detect_local_maxima(img: np.ndarray, min_distance: int = 5, threshold_rel: float = 0.55, max_points: int = 80) -> np.ndarray:
    """Simple local-max detection using scipy only."""
    work = np.asarray(img, dtype=np.float64)
    threshold = threshold_rel * float(np.max(work))
    maxf = ndi.maximum_filter(work, size=max(3, min_distance * 2 + 1), mode="reflect")
    peaks = (work == maxf) & (work >= threshold)
    coords = np.argwhere(peaks)
    if coords.size == 0:
        return np.empty((0, 2), dtype=int)
    values = work[coords[:, 0], coords[:, 1]]
    order = np.argsort(values)[::-1]
    coords = coords[order]

    selected: List[np.ndarray] = []
    for rc in coords:
        if len(selected) >= max_points:
            break
        if all(np.linalg.norm(rc - s) >= min_distance for s in selected):
            selected.append(rc)
    return np.array(selected, dtype=int) if selected else np.empty((0, 2), dtype=int)


# ---------------------------
# ROI loading
# ---------------------------

def discover_roi_files(roi_dir: Path) -> Dict[str, Path]:
    if not roi_dir.exists():
        raise FileNotFoundError(f"ROI directory not found: {roi_dir}")

    mapping: Dict[str, Path] = {}
    for path in sorted(roi_dir.iterdir()):
        if path.suffix.lower() not in IMAGE_EXTENSIONS:
            continue
        stem = path.stem.strip()
        upper = stem.upper().replace("_", "-")
        for roi in ROI_NAME_ORDER:
            if upper == roi.upper() or upper.endswith(roi.upper()) or roi.upper() in upper:
                mapping[roi] = path
                break

    missing = [roi for roi in ROI_NAME_ORDER if roi not in mapping]
    if missing:
        raise FileNotFoundError(
            f"Missing ROI image files for: {missing}. Expected names containing ROI-1 ... ROI-5 in {roi_dir}."
        )
    return mapping


def read_roi_box_table(roi_csv: Path, image_shape: Tuple[int, int]) -> pd.DataFrame:
    df = pd.read_csv(roi_csv)
    if "roi" not in df.columns and "roi_name" in df.columns:
        df = df.rename(columns={"roi_name": "roi"})
    required = {"roi", "x", "y", "width", "height"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"ROI CSV missing required columns: {sorted(missing)}")

    height, width = image_shape
    rows = []
    seen_rois = set()
    seen_boxes = set()
    for _, row in df.iterrows():
        roi = str(row["roi"]).strip()
        if roi not in ROI_NAME_ORDER:
            raise ValueError(f"Unexpected ROI label in ROI CSV: {roi}")
        if roi in seen_rois:
            raise ValueError(f"Duplicate ROI label in ROI CSV: {roi}")
        seen_rois.add(roi)

        x = int(row["x"])
        y = int(row["y"])
        w = int(row["width"])
        h = int(row["height"])
        if x < 0 or y < 0 or w <= 0 or h <= 0:
            raise ValueError(f"Invalid ROI box for {roi}: x={x}, y={y}, width={w}, height={h}")
        if x + w > width or y + h > height:
            raise ValueError(
                f"ROI box for {roi} exceeds source image bounds: "
                f"x={x}, y={y}, width={w}, height={h}, image={width}x{height}"
            )

        box = (x, y, w, h)
        if box in seen_boxes:
            raise ValueError(f"Duplicate ROI crop box in ROI CSV: {box}")
        seen_boxes.add(box)
        rows.append({"roi": roi, "x": x, "y": y, "width": w, "height": h})

    missing_order = [roi for roi in ROI_NAME_ORDER if roi not in seen_rois]
    if missing_order:
        raise ValueError(f"ROI CSV did not provide all expected ROIs: {missing_order}")

    ordered = pd.DataFrame(rows).set_index("roi").loc[ROI_NAME_ORDER].reset_index()
    return ordered


def crop_rois_from_source(source_image: Path, roi_csv: Path) -> Dict[str, np.ndarray]:
    img = load_image(source_image)
    df = read_roi_box_table(roi_csv, image_shape=img.shape[:2])

    rois: Dict[str, np.ndarray] = {}
    for _, row in df.iterrows():
        roi = str(row["roi"]).strip()
        x = int(row["x"])
        y = int(row["y"])
        w = int(row["width"])
        h = int(row["height"])
        crop = img[y:y+h, x:x+w]
        if crop.size == 0:
            raise ValueError(f"Empty crop for {roi}: x={x}, y={y}, width={w}, height={h}")
        rois[roi] = crop

    return rois


def load_rois(repo_root: Path, roi_dir: Optional[Path], source_image: Optional[Path], roi_csv: Optional[Path]) -> Dict[str, np.ndarray]:
    if roi_dir is not None:
        files = discover_roi_files(repo_root / roi_dir if not roi_dir.is_absolute() else roi_dir)
        return {roi: load_image(path) for roi, path in files.items()}

    if source_image is not None and roi_csv is not None:
        src = repo_root / source_image if not source_image.is_absolute() else source_image
        csv = repo_root / roi_csv if not roi_csv.is_absolute() else roi_csv
        return crop_rois_from_source(src, csv)

    default_dir = repo_root / "data" / "raw" / "TEM" / "rois"
    if default_dir.exists():
        files = discover_roi_files(default_dir)
        return {roi: load_image(path) for roi, path in files.items()}

    raise ValueError(
        "No TEM input mode found. Provide --roi-dir or both --source-image and --roi-csv."
    )


# ---------------------------
# Analysis pipeline
# ---------------------------

def analyze_single_roi(roi_name: str, img: np.ndarray) -> ROIResult:
    proc = preprocess_roi(img)
    fft_log = fft_log_magnitude(proc)
    q, prof = radial_profile(fft_log)
    peak_radius, p_fft, dq_fft = estimate_peak_and_width(q, prof)
    recon = reconstruct_bandpass(proc, peak_radius=peak_radius, dq_fft=dq_fft)
    a_off = autocorr_offcenter_amplitude(recon)
    r_bp = bandpass_rms(recon)
    maxima = detect_local_maxima(recon)

    return ROIResult(
        roi=roi_name,
        raw=img,
        processed=proc,
        fft_log=fft_log,
        bandpass_recon=recon,
        maxima_rc=maxima,
        p_fft=p_fft,
        a_off=a_off,
        r_bp=r_bp,
        dq_fft=dq_fft,
        peak_radius=peak_radius,
        radial_profile_q=q,
        radial_profile_val=prof,
    )


def results_to_tables(results: Sequence[ROIResult], geometric_eps: float) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, float]]:
    if not np.isfinite(geometric_eps) or geometric_eps <= 0:
        raise ValueError("geometric_eps must be a finite value greater than 0")

    df = pd.DataFrame(
        {
            "ROI": [r.roi for r in results],
            "P_FFT": [r.p_fft for r in results],
            "A_off": [r.a_off for r in results],
            "R_bp": [r.r_bp for r in results],
            "Delta_q_FFT": [r.dq_fft for r in results],
        }
    )

    Pn = minmax(df["P_FFT"].to_numpy())
    An = minmax(df["A_off"].to_numpy())
    Rn = minmax(df["R_bp"].to_numpy())
    delta_q_norm_raw = minmax(df["Delta_q_FFT"].to_numpy())
    Cn = 1.0 - delta_q_norm_raw

    loi = 100.0 * (Pn + An + Rn + Cn) / 4.0
    g = ((Pn + geometric_eps) * (An + geometric_eps) * (Rn + geometric_eps) * (Cn + geometric_eps)) ** 0.25
    weight_sum = np.sum(g)
    if not np.isfinite(weight_sum) or weight_sum <= 0:
        raise ValueError("geometric patch weight denominator must be positive and finite")
    wi = g / weight_sum

    df["LOI_0_100"] = loi
    df["geometric_patch_weight_wi"] = wi

    rank_order = df.sort_values("LOI_0_100", ascending=False)["ROI"].tolist()
    ranking_str = " > ".join(rank_order)

    centroid_loi = float(np.sum(df["LOI_0_100"].to_numpy() * wi))
    centroid_dq = float(np.sum(df["Delta_q_FFT"].to_numpy() * wi))

    df_norm = pd.DataFrame(
        {
            "ROI": df["ROI"],
            "P_FFT_norm": Pn,
            "A_off_norm": An,
            "R_bp_norm": Rn,
            "Delta_q_FFT_norm_raw": delta_q_norm_raw,
            "Delta_q_FFT_ordering_tendency_inverted": Cn,
            "LOI_0_100": loi,
            "geometric_patch_weight_wi": wi,
        }
    )

    metadata = {
        "geometric_eps": float(geometric_eps),
        "ensemble_centroid_LOI": centroid_loi,
        "ensemble_centroid_Delta_q_FFT": centroid_dq,
        "ordering_rank": ranking_str,
    }
    return df, df_norm, metadata


def reference_descriptors_to_tables(
    descriptor_csv: Path,
    geometric_eps: float,
    normalized_csv: Optional[Path] = None,
    repo_root: Optional[Path] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, Dict[str, object]]:
    if not np.isfinite(geometric_eps) or geometric_eps <= 0:
        raise ValueError("geometric_eps must be a finite value greater than 0")
    if not descriptor_csv.exists():
        raise FileNotFoundError(f"Reference TEM descriptor table not found: {descriptor_csv}")

    df = pd.read_csv(descriptor_csv)
    missing = [column for column in DESCRIPTOR_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"Reference TEM descriptor table missing columns: {missing}")

    df = df[DESCRIPTOR_COLUMNS].copy()
    df["ROI"] = df["ROI"].astype(str).str.strip()
    if set(df["ROI"]) != set(ROI_NAME_ORDER):
        raise ValueError(f"Reference TEM descriptor table must contain exactly {ROI_NAME_ORDER}")
    df = df.set_index("ROI").loc[ROI_NAME_ORDER].reset_index()

    wi = df["geometric_patch_weight_wi"].to_numpy(dtype=float)
    weight_sum = float(np.sum(wi))
    if not np.isfinite(weight_sum) or weight_sum <= 0:
        raise ValueError("Reference TEM geometric patch weights must sum to a positive finite value")
    wi_norm = wi / weight_sum

    if normalized_csv is not None and normalized_csv.exists():
        df_norm = pd.read_csv(normalized_csv)
        missing_norm = [column for column in NORMALIZED_DESCRIPTOR_COLUMNS if column not in df_norm.columns]
        if missing_norm:
            raise ValueError(f"Reference TEM normalized descriptor table missing columns: {missing_norm}")
        df_norm = df_norm[NORMALIZED_DESCRIPTOR_COLUMNS].copy()
        df_norm["ROI"] = df_norm["ROI"].astype(str).str.strip()
        if set(df_norm["ROI"]) != set(ROI_NAME_ORDER):
            raise ValueError(f"Reference TEM normalized descriptor table must contain exactly {ROI_NAME_ORDER}")
        df_norm = df_norm.set_index("ROI").loc[ROI_NAME_ORDER].reset_index()
    else:
        Pn = minmax(df["P_FFT"].to_numpy())
        An = minmax(df["A_off"].to_numpy())
        Rn = minmax(df["R_bp"].to_numpy())
        delta_q_norm_raw = minmax(df["Delta_q_FFT"].to_numpy())
        Cn = 1.0 - delta_q_norm_raw
        df_norm = pd.DataFrame(
            {
                "ROI": df["ROI"],
                "P_FFT_norm": Pn,
                "A_off_norm": An,
                "R_bp_norm": Rn,
                "Delta_q_FFT_norm_raw": delta_q_norm_raw,
                "Delta_q_FFT_ordering_tendency_inverted": Cn,
                "LOI_0_100": df["LOI_0_100"],
                "geometric_patch_weight_wi": df["geometric_patch_weight_wi"],
            }
        )

    metadata = {
        "geometric_eps": float(geometric_eps),
        "descriptor_source": portable_source_path(descriptor_csv, repo_root=repo_root),
        "ensemble_centroid_LOI": float(np.sum(df["LOI_0_100"].to_numpy(dtype=float) * wi_norm)),
        "ensemble_centroid_Delta_q_FFT": float(np.sum(df["Delta_q_FFT"].to_numpy(dtype=float) * wi_norm)),
        "ordering_rank": " > ".join(df.sort_values("LOI_0_100", ascending=False)["ROI"].tolist()),
    }
    return df, df_norm, metadata


# ---------------------------
# Plotting
# ---------------------------

def set_plot_style(font_scale: float = 1.0) -> None:
    plt.rcParams.update(
        {
            "font.size": 12 * font_scale,
            "axes.labelsize": 14 * font_scale,
            "axes.titlesize": 14 * font_scale,
            "xtick.labelsize": 12 * font_scale,
            "ytick.labelsize": 12 * font_scale,
            "axes.linewidth": 1.3,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", top=True, right=True, width=1.2, length=5)
    for spine in ax.spines.values():
        spine.set_linewidth(1.3)


def save_roi_artifacts(result: ROIResult, out_dir: Path) -> None:
    roi_dir = out_dir / result.roi
    ensure_dir(roi_dir)

    def save_img(img: np.ndarray, filename: str, cmap: str = "gray") -> None:
        fig, ax = plt.subplots(figsize=(3.0, 3.0), dpi=250)
        ax.imshow(img, cmap=cmap)
        ax.set_axis_off()
        fig.tight_layout(pad=0)
        fig.savefig(roi_dir / filename, bbox_inches="tight", pad_inches=0)
        plt.close(fig)

    save_img(result.raw, "raw.png")
    save_img(result.processed, "processed.png")
    save_img(result.fft_log, "fft_log.png")
    save_img(result.bandpass_recon, "bandpass_reconstruction.png")

    fig, ax = plt.subplots(figsize=(3.0, 3.0), dpi=250)
    ax.imshow(result.bandpass_recon, cmap="gray")
    if len(result.maxima_rc):
        ax.scatter(result.maxima_rc[:, 1], result.maxima_rc[:, 0], s=10, marker="o", facecolors="none", edgecolors="tab:red", linewidths=0.6)
    ax.set_axis_off()
    fig.tight_layout(pad=0)
    fig.savefig(roi_dir / "bandpass_maxima_overlay.png", bbox_inches="tight", pad_inches=0)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.0, 3.0), dpi=250)
    ax.plot(result.radial_profile_q, result.radial_profile_val)
    ax.set_xlabel("q (pixel-radius)")
    ax.set_ylabel("Radial FFT intensity")
    ax.set_title(result.roi)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(roi_dir / "radial_fft_profile.png", bbox_inches="tight")
    plt.close(fig)


def save_roi_selection_overlay(source_image: Path, roi_csv: Path, out_dir: Path) -> None:
    img = load_image(source_image)
    boxes = read_roi_box_table(roi_csv, image_shape=img.shape[:2])

    set_plot_style(1.0)
    fig, ax = plt.subplots(figsize=(7.2, 7.2), dpi=300)
    ax.imshow(img, cmap="gray")
    for _, row in boxes.iterrows():
        x = int(row["x"])
        y = int(row["y"])
        w = int(row["width"])
        h = int(row["height"])
        roi = str(row["roi"])
        ax.add_patch(Rectangle((x, y), w, h, fill=False, edgecolor="#f2c400", linewidth=2.2))
        ax.text(
            x,
            max(0, y - 24),
            roi,
            ha="left",
            va="bottom",
            fontsize=9,
            fontweight="bold",
            color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5},
        )
    ax.set_axis_off()
    fig.tight_layout(pad=0)
    fig.savefig(out_dir / "TEM_roi_selection_overlay.png", bbox_inches="tight", pad_inches=0)
    fig.savefig(out_dir / "TEM_roi_selection_overlay.pdf", bbox_inches="tight", pad_inches=0)
    plt.close(fig)



def plot_figure_1g(df: pd.DataFrame, metadata: Dict[str, object], out_dir: Path) -> None:
    set_plot_style(1.15)
    fig, ax = plt.subplots(figsize=(6.0, 5.2), dpi=300)

    x = df["LOI_0_100"].to_numpy()
    y = df["Delta_q_FFT"].to_numpy()
    ax.scatter(x, y, s=70)

    cx = float(metadata["ensemble_centroid_LOI"])
    cy = float(metadata["ensemble_centroid_Delta_q_FFT"])
    ax.scatter([cx], [cy], marker="*", s=180)

    median_x = float(np.median(x))
    median_y = float(np.median(y))
    ax.axvline(median_x, linestyle="--", linewidth=1.0)
    ax.axhline(median_y, linestyle="--", linewidth=1.0)

    ax.set_xlabel("Local Ordering Index, LOI (0–100)")
    ax.set_ylabel(r"$\Delta q_{\mathrm{FFT}}$ (a.u.; peak width)")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "TEM_order_disorder_map.png", bbox_inches="tight")
    fig.savefig(out_dir / "TEM_order_disorder_map.pdf", bbox_inches="tight")
    plt.close(fig)



def plot_figure_1h(df: pd.DataFrame, out_dir: Path) -> None:
    set_plot_style(1.15)
    fig, ax = plt.subplots(figsize=(6.0, 5.2), dpi=300)
    bars = ax.bar(df["ROI"], df["geometric_patch_weight_wi"])
    ax.set_xlabel("TEM ROI")
    ax.set_ylabel(r"Geometric patch weight, $w_i$")
    for bar, val in zip(bars, df["geometric_patch_weight_wi"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01, f"{val:.3f}", ha="center", va="bottom", fontsize=10)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "TEM_geometric_patch_weights.png", bbox_inches="tight")
    fig.savefig(out_dir / "TEM_geometric_patch_weights.pdf", bbox_inches="tight")
    plt.close(fig)



def plot_figure_1i(df: pd.DataFrame, out_dir: Path) -> None:
    set_plot_style(1.10)

    matrix_raw = np.column_stack(
        [
            df["P_FFT"].to_numpy(),
            df["A_off"].to_numpy(),
            df["R_bp"].to_numpy(),
            df["Delta_q_FFT"].to_numpy(),
            df["LOI_0_100"].to_numpy(),
        ]
    )

    matrix_norm = np.column_stack(
        [
            minmax(df["P_FFT"].to_numpy()),
            minmax(df["A_off"].to_numpy()),
            minmax(df["R_bp"].to_numpy()),
            1.0 - minmax(df["Delta_q_FFT"].to_numpy()),
            minmax(df["LOI_0_100"].to_numpy()),
        ]
    )

    fig, ax = plt.subplots(figsize=(7.4, 4.8), dpi=300)
    im = ax.imshow(matrix_norm, aspect="auto", cmap="viridis", norm=Normalize(vmin=0.0, vmax=1.0))

    ax.set_xticks(range(5))
    ax.set_xticklabels([r"P$_{\mathrm{FFT}}$", r"A$_{\mathrm{off}}$", r"R$_{\mathrm{bp}}$", r"$\Delta q_{\mathrm{FFT}}$", "LOI"])
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["ROI"])

    for i in range(matrix_raw.shape[0]):
        for j in range(matrix_raw.shape[1]):
            val = matrix_raw[i, j]
            text = f"{val:.3f}" if j < 4 else f"{val:.1f}"
            ax.text(j, i, text, ha="center", va="center", fontsize=9, color="white" if matrix_norm[i, j] < 0.35 else "black")

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Normalized ordering tendency")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(out_dir / "TEM_descriptor_heatmap.png", bbox_inches="tight")
    fig.savefig(out_dir / "TEM_descriptor_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)


# ---------------------------
# Main CLI
# ---------------------------

def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Rebuild TEM ROI descriptors and ensemble panels for the h-BN repository.")
    parser.add_argument("--repo-root", type=Path, default=Path("."), help="Repository root directory.")
    parser.add_argument("--roi-dir", type=Path, default=None, help="Directory containing cropped ROI images.")
    parser.add_argument("--source-image", type=Path, default=None, help="Path to a single source TEM image.")
    parser.add_argument("--roi-csv", type=Path, default=None, help="CSV containing ROI crop boxes for the source image.")
    parser.add_argument(
        "--descriptor-csv",
        type=Path,
        default=Path("data/processed/TEM/TEM_descriptors.csv"),
        help="Manuscript reference descriptor table used for summary panels unless --recompute-descriptors is set.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("data/processed/TEM"), help="Output directory relative to repo root.")
    parser.add_argument("--geometric-eps", type=float, default=0.05, help="Regularization floor used in the geometric ensemble score.")
    parser.add_argument(
        "--recompute-descriptors",
        action="store_true",
        help="Recompute exploratory descriptors from the current ROI crops instead of using the manuscript reference table.",
    )
    parser.add_argument("--save-json-summary", action="store_true", help="Also save a JSON summary of the centroid and rank order.")
    return parser



def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    if not math.isfinite(args.geometric_eps) or args.geometric_eps <= 0:
        parser.error("--geometric-eps must be a positive finite float")

    repo_root = args.repo_root.resolve()
    out_dir = args.output_dir if args.output_dir.is_absolute() else repo_root / args.output_dir
    ensure_dir(out_dir)
    descriptor_csv = args.descriptor_csv if args.descriptor_csv.is_absolute() else repo_root / args.descriptor_csv
    normalized_descriptor_csv = descriptor_csv.with_name("TEM_descriptors_normalized.csv")

    rois = load_rois(repo_root=repo_root, roi_dir=args.roi_dir, source_image=args.source_image, roi_csv=args.roi_csv)
    ordered_rois = {roi: rois[roi] for roi in ROI_NAME_ORDER}
    if args.source_image is not None and args.roi_csv is not None:
        source_image = repo_root / args.source_image if not args.source_image.is_absolute() else args.source_image
        roi_csv = repo_root / args.roi_csv if not args.roi_csv.is_absolute() else args.roi_csv
        save_roi_selection_overlay(source_image, roi_csv, out_dir)

    results: List[ROIResult] = []
    for roi_name, img in ordered_rois.items():
        result = analyze_single_roi(roi_name, img)
        results.append(result)
        save_roi_artifacts(result, out_dir)

    if args.recompute_descriptors:
        df, df_norm, metadata = results_to_tables(results, geometric_eps=args.geometric_eps)
        metadata["descriptor_source"] = "recomputed_from_current_roi_crops"
    else:
        df, df_norm, metadata = reference_descriptors_to_tables(
            descriptor_csv,
            geometric_eps=args.geometric_eps,
            normalized_csv=normalized_descriptor_csv,
            repo_root=repo_root,
        )

    descriptor_out = out_dir / "TEM_descriptors.csv"
    normalized_descriptor_out = out_dir / "TEM_descriptors_normalized.csv"
    if args.recompute_descriptors:
        df.to_csv(descriptor_out, index=False)
        df_norm.to_csv(normalized_descriptor_out, index=False)
    else:
        copy_file_if_different(descriptor_csv, descriptor_out)
        if normalized_descriptor_csv.exists():
            copy_file_if_different(normalized_descriptor_csv, normalized_descriptor_out)
        else:
            df_norm.to_csv(normalized_descriptor_out, index=False)

    if args.save_json_summary:
        with open(out_dir / "TEM_ensemble_summary.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)

    plot_figure_1g(df, metadata, out_dir)
    plot_figure_1h(df, out_dir)
    plot_figure_1i(df, out_dir)

    print("TEM analysis completed.")
    print(f"Output directory: {out_dir}")
    print("Generated files:")
    for name in [
        "TEM_descriptors.csv",
        "TEM_descriptors_normalized.csv",
        "TEM_roi_selection_overlay.png",
        "TEM_order_disorder_map.png",
        "TEM_geometric_patch_weights.png",
        "TEM_descriptor_heatmap.png",
    ]:
        print(f"  - {out_dir / name}")
    print(f"Descriptor source: {metadata.get('descriptor_source', 'recomputed_from_current_roi_crops')}")
    print(f"Ensemble centroid: LOI = {metadata['ensemble_centroid_LOI']:.3f}, Delta_q_FFT = {metadata['ensemble_centroid_Delta_q_FFT']:.3f}")
    print(f"Ordering rank: {metadata['ordering_rank']}")


if __name__ == "__main__":
    main()
