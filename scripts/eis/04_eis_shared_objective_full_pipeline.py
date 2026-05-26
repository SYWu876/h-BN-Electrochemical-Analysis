#!/usr/bin/env python3
"""
04_eis_shared_objective_full_pipeline.py

End-to-end, shared-objective EIS pipeline for the h-BN manuscript repository.

This script starts from a raw EIS CSV file, performs a classical complex-domain
least-squares fit of the compact non-ideal circuit,

    Z(w) = Rs + j*w*L + (Rct || CPE1) + CPE2,

builds a local quadratic surrogate around the classical anchor in latent
parameter space (log coordinates for positive parameters and logit coordinates
for alpha values), solves the continuous surrogate branch, discretizes that
local trust region using 3 bits per parameter, maps the surrogate to a
QUBO/Ising-style binary objective, evaluates an exact p=1 statevector QAOA
angle landscape over the discretized trust-region lattice, decodes the best
sampled bitstring, and exports spectra, tables,
surrogate slices, landscape tables, and decoded results. The exported slice
columns named normalized_x/y are local latent displacements kept for
compatibility with the release package. A fast diagnostic proxy is available
through --fast-proxy-qaoa, but the default qaoa_p1_landscape.csv/png outputs are
exact p=1 statevector calculations.

The script is intentionally self-contained and transparent. It is intended for
reproducibility and review, not for claiming quantum advantage.

Example
-------
python scripts/eis/04_eis_shared_objective_full_pipeline.py \
    --input data/raw/EIS/hbn_EIS_1.csv \
    --output outputs/eis_shared_objective \
    --trust-delta 0.08 \
    --bits-per-parameter 3 \
    --qaoa-grid 9

Dependencies: numpy, pandas, scipy, matplotlib
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize
import matplotlib.pyplot as plt


PARAM_NAMES = ["Rs", "L", "Rct", "Q1", "alpha1", "Q2", "alpha2"]
LOG_PARAMS = {"Rs", "L", "Rct", "Q1", "Q2"}
ALPHA_PARAMS = {"alpha1", "alpha2"}


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


def bits_per_parameter(value: str) -> int:
    parsed = positive_int(value)
    if parsed > 3:
        raise argparse.ArgumentTypeError("must be between 1 and 3 for this lightweight release pipeline")
    return parsed


def read_eis_csv(path: Path) -> pd.DataFrame:
    """Read raw EIS CSV and standardize key columns."""
    df = pd.read_csv(path)
    required = ["Frequency (Hz)", "Z' (Ohms)", "-Z\" (Ohms)"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required EIS columns: {missing}")

    out = pd.DataFrame(
        {
            "frequency_Hz": pd.to_numeric(df["Frequency (Hz)"], errors="coerce"),
            "Zreal_ohm": pd.to_numeric(df["Z' (Ohms)"], errors="coerce"),
            "minus_Zimag_ohm": pd.to_numeric(df["-Z\" (Ohms)"], errors="coerce"),
        }
    ).dropna()
    out = out[out["frequency_Hz"] > 0].copy()
    out = out.sort_values("frequency_Hz", ascending=False).reset_index(drop=True)
    return out


def cpe_impedance(Q: float, alpha: float, omega: np.ndarray) -> np.ndarray:
    """Constant phase element impedance: Z_CPE = 1/(Q*(jw)^alpha)."""
    return 1.0 / (Q * (1j * omega) ** alpha)


def circuit_impedance(params: Sequence[float], freq_hz: np.ndarray) -> np.ndarray:
    """Compact EIS model: Rs + jwL + (Rct || CPE1) + CPE2."""
    Rs, L, Rct, Q1, alpha1, Q2, alpha2 = params
    omega = 2.0 * np.pi * freq_hz
    Z_L = 1j * omega * L
    Z_Q1 = cpe_impedance(Q1, alpha1, omega)
    Z_parallel = 1.0 / (1.0 / Rct + 1.0 / Z_Q1)
    Z_Q2 = cpe_impedance(Q2, alpha2, omega)
    return Rs + Z_L + Z_parallel + Z_Q2


def physical_to_latent(params: Sequence[float]) -> np.ndarray:
    """Map physical parameters to unconstrained latent coordinates."""
    if len(params) != len(PARAM_NAMES):
        raise ValueError(f"Expected {len(PARAM_NAMES)} physical parameters, got {len(params)}")
    x = []
    for name, value in zip(PARAM_NAMES, params):
        if name in LOG_PARAMS:
            x.append(np.log(max(float(value), 1e-300)))
        elif name in ALPHA_PARAMS:
            v = float(np.clip(value, 1e-9, 1 - 1e-9))
            x.append(np.log(v / (1 - v)))
        else:
            raise ValueError(name)
    return np.array(x, dtype=float)


def latent_to_physical(x: Sequence[float]) -> np.ndarray:
    """Map unconstrained latent coordinates to physical parameters."""
    if len(x) != len(PARAM_NAMES):
        raise ValueError(f"Expected {len(PARAM_NAMES)} latent coordinates, got {len(x)}")
    out = []
    for name, value in zip(PARAM_NAMES, x):
        if name in LOG_PARAMS:
            out.append(float(np.exp(value)))
        elif name in ALPHA_PARAMS:
            # Stable sigmoid
            if value >= 0:
                z = np.exp(-value)
                out.append(float(1.0 / (1.0 + z)))
            else:
                z = np.exp(value)
                out.append(float(z / (1.0 + z)))
        else:
            raise ValueError(name)
    return np.array(out, dtype=float)


def complex_residual(
    x: Sequence[float], freq_hz: np.ndarray, z_exp: np.ndarray, weighted: bool = False
) -> np.ndarray:
    """Stacked real/imaginary residuals for least-squares fitting.

    The default unweighted residual matches the classical anchor script and
    gives the Nyquist overlay used for visual fit assessment.  The weighted
    form is retained only as a relative-error diagnostic.
    """
    params = latent_to_physical(x)
    z_model = circuit_impedance(params, freq_hz)
    if weighted:
        scale = np.maximum(np.abs(z_exp), 1.0)
    else:
        scale = np.ones_like(np.abs(z_exp))
    return np.r_[(z_model.real - z_exp.real) / scale, (z_model.imag - z_exp.imag) / scale]


def loss_from_latent(
    x: Sequence[float], freq_hz: np.ndarray, z_exp: np.ndarray, weighted: bool = False
) -> float:
    r = complex_residual(x, freq_hz, z_exp, weighted=weighted)
    return float(np.sum(r * r))


def fit_classical_anchor(freq_hz: np.ndarray, z_exp: np.ndarray) -> Tuple[np.ndarray, Dict[str, float]]:
    """Fit the compact EIS model from several physically plausible starts."""
    initial_guesses = [
        [1.3, 1e-7, 34.0, 0.05, 0.87, 0.04, 0.91],
        [1.33, 1.83e-7, 34.0, 0.00331, 0.874, 0.00231, 0.915],
        [1.0, 1e-7, 100.0, 0.005, 0.75, 0.003, 0.95],
        [1.0, 1e-6, 40.0, 0.02, 0.80, 0.02, 0.90],
        [2.0, 1e-7, 80.0, 0.01, 0.80, 0.005, 0.90],
    ]

    best = None
    for guess in initial_guesses:
        x0 = physical_to_latent(guess)
        result = least_squares(
            complex_residual,
            x0,
            args=(freq_hz, z_exp, False),
            method="trf",
            max_nfev=30000,
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
        )
        if best is None or result.cost < best.cost:
            best = result

    if best is None:
        raise RuntimeError("Classical fitting failed to start.")

    x_anchor = np.asarray(best.x, dtype=float)
    params = latent_to_physical(x_anchor)
    stats = {
        "least_squares_cost": float(best.cost),
        "unweighted_SSE": loss_from_latent(x_anchor, freq_hz, z_exp, weighted=False),
        "relative_weighted_SSE": loss_from_latent(x_anchor, freq_hz, z_exp, weighted=True),
        "nfev": int(best.nfev),
        "optimality": float(best.optimality),
        "success": bool(best.success),
    }
    return params, stats


def local_displacement_to_latent(displacement: np.ndarray, x_anchor: np.ndarray) -> np.ndarray:
    """Map a local latent-space displacement to absolute latent coordinates."""
    return x_anchor + np.asarray(displacement, dtype=float)


def build_quadratic_surrogate(
    x_anchor: np.ndarray,
    freq_hz: np.ndarray,
    z_exp: np.ndarray,
    fd_step: float,
) -> Tuple[float, np.ndarray, np.ndarray]:
    """Finite-difference quadratic surrogate in local latent coordinates y.

    E(y) = c + g^T y + 0.5 y^T H y
    where y = x - x_anchor.
    """
    n = len(x_anchor)
    c = loss_from_latent(x_anchor, freq_hz, z_exp)
    g = np.zeros(n)
    H = np.zeros((n, n))

    eye = np.eye(n)
    for i in range(n):
        fp = loss_from_latent(x_anchor + fd_step * eye[i], freq_hz, z_exp)
        fm = loss_from_latent(x_anchor - fd_step * eye[i], freq_hz, z_exp)
        g[i] = (fp - fm) / (2.0 * fd_step)
        H[i, i] = (fp - 2.0 * c + fm) / (fd_step * fd_step)

    for i in range(n):
        for j in range(i + 1, n):
            fpp = loss_from_latent(x_anchor + fd_step * eye[i] + fd_step * eye[j], freq_hz, z_exp)
            fpm = loss_from_latent(x_anchor + fd_step * eye[i] - fd_step * eye[j], freq_hz, z_exp)
            fmp = loss_from_latent(x_anchor - fd_step * eye[i] + fd_step * eye[j], freq_hz, z_exp)
            fmm = loss_from_latent(x_anchor - fd_step * eye[i] - fd_step * eye[j], freq_hz, z_exp)
            Hij = (fpp - fpm - fmp + fmm) / (4.0 * fd_step * fd_step)
            H[i, j] = H[j, i] = Hij

    return c, g, H


def surrogate_energy(y: np.ndarray, c: float, g: np.ndarray, H: np.ndarray) -> np.ndarray:
    """Evaluate E(y)=c+g*y+0.5*y'H*y for one or many y vectors."""
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        return float(c + g @ y + 0.5 * y @ H @ y)
    return c + y @ g + 0.5 * np.einsum("bi,ij,bj->b", y, H, y)


def levels_for_bits(bits: int, trust_delta: float) -> np.ndarray:
    return np.linspace(-trust_delta, trust_delta, 2**bits)


def enumerate_lattice_minimum(
    c: float, g: np.ndarray, H: np.ndarray, levels: np.ndarray, chunk_size: int = 200000
) -> Tuple[np.ndarray, float, int]:
    """Enumerate the 8^7 local lattice efficiently and return best indices."""
    n = len(g)
    total = len(levels) ** n
    best_E = np.inf
    best_indices = None
    best_linear_index = -1

    # chunk over Cartesian product to avoid large permanent arrays
    batch_idx = []
    start_index = 0
    for linear_index, idx in enumerate(itertools.product(range(len(levels)), repeat=n)):
        batch_idx.append(idx)
        if len(batch_idx) >= chunk_size:
            arr_idx = np.asarray(batch_idx, dtype=int)
            y = levels[arr_idx]
            E = surrogate_energy(y, c, g, H)
            k = int(np.argmin(E))
            if E[k] < best_E:
                best_E = float(E[k])
                best_indices = arr_idx[k].copy()
                best_linear_index = start_index + k
            start_index += len(batch_idx)
            batch_idx = []
    if batch_idx:
        arr_idx = np.asarray(batch_idx, dtype=int)
        y = levels[arr_idx]
        E = surrogate_energy(y, c, g, H)
        k = int(np.argmin(E))
        if E[k] < best_E:
            best_E = float(E[k])
            best_indices = arr_idx[k].copy()
            best_linear_index = start_index + k

    if best_indices is None:
        raise RuntimeError("Lattice enumeration failed.")
    return best_indices, best_E, best_linear_index


def index_to_bits(index: int, bits: int) -> List[int]:
    return [(index >> k) & 1 for k in range(bits)]


def bits_to_index(bits_list: Sequence[int]) -> int:
    return int(sum(int(b) << k for k, b in enumerate(bits_list)))


def lattice_indices_to_bitstring(indices: Sequence[int], bits: int) -> str:
    # Human-readable: per parameter, most significant bit first.
    pieces = []
    for idx in indices:
        b = index_to_bits(int(idx), bits)
        pieces.append("".join(str(x) for x in b[::-1]))
    return " ".join(pieces)


def decode_state_level_indices(states: np.ndarray, n_params: int, bits: int) -> np.ndarray:
    """Decode integer bitstring states to lattice-level indices in bounded chunks."""
    decoded = np.zeros((states.size, n_params), dtype=np.int16)
    for p in range(n_params):
        val = np.zeros(states.size, dtype=np.int16)
        for b in range(bits):
            bit_pos = p * bits + b
            val += (((states >> bit_pos) & 1).astype(np.int16) << b)
        decoded[:, p] = val
    return decoded


def diagonal_surrogate_energies(
    c: float, g: np.ndarray, H: np.ndarray, levels: np.ndarray, bits: int, chunk_size: int = 200000
) -> np.ndarray:
    """Compute diagonal energies for every bitstring state with bounded peak memory."""
    n_params = len(g)
    nbits = n_params * bits
    nstates = 1 << nbits
    energies = np.empty(nstates, dtype=np.float64)

    for start in range(0, nstates, chunk_size):
        stop = min(start + chunk_size, nstates)
        states = np.arange(start, stop, dtype=np.uint64)
        decoded = decode_state_level_indices(states, n_params, bits)
        y = levels[decoded]
        energies[start:stop] = surrogate_energy(y, c, g, H)

    return energies


def apply_rx_all_qubits(state: np.ndarray, beta: float, nbits: int) -> np.ndarray:
    """Apply product exp(-i beta X) to every qubit using vectorized blocks."""
    c = math.cos(beta)
    s = -1j * math.sin(beta)
    psi = state.copy()
    for q in range(nbits):
        step = 1 << q
        paired = psi.reshape(-1, 2, step)
        a = paired[:, 0, :].copy()
        b = paired[:, 1, :].copy()
        paired[:, 0, :] = c * a + s * b
        paired[:, 1, :] = s * a + c * b
    return psi


def qaoa_p1_probabilities(diag_E: np.ndarray, gamma: float, beta: float) -> np.ndarray:
    """Exact statevector probabilities for p=1 QAOA on a diagonal cost Hamiltonian."""
    nstates = diag_E.size
    nbits = int(round(math.log2(nstates)))
    E_centered = diag_E - np.mean(diag_E)
    scale = np.std(E_centered)
    if scale > 0:
        E_phase = E_centered / scale
    else:
        E_phase = E_centered
    psi = np.exp(-1j * gamma * E_phase).astype(np.complex128) / math.sqrt(nstates)
    psi = apply_rx_all_qubits(psi, beta, nbits)
    prob = np.abs(psi) ** 2
    return prob / prob.sum()


def qaoa_p1_expectation(diag_E: np.ndarray, gamma: float, beta: float) -> float:
    """Exact statevector expectation for p=1 QAOA on a diagonal cost Hamiltonian."""
    prob = qaoa_p1_probabilities(diag_E, gamma, beta)
    return float(np.dot(prob, diag_E))


def scaled_proxy_energies(diag_E: np.ndarray) -> Tuple[np.ndarray, np.ndarray, float]:
    E = np.asarray(diag_E, dtype=float)
    Esc = E - E.min()
    scale = np.std(Esc)
    if scale > 0:
        Esc = Esc / scale
    return E, Esc, float(E.mean())


def proxy_cost_distribution(scaled_E: np.ndarray, gamma: float) -> np.ndarray:
    # stable softmin selected by gamma
    logits = -float(gamma) * scaled_E
    logits -= logits.max()
    p_cost = np.exp(logits)
    p_cost /= p_cost.sum()
    return p_cost


def proxy_qaoa_probabilities(diag_E: np.ndarray, gamma: float, beta: float) -> np.ndarray:
    """Fast diagnostic proxy distribution over the binary lattice.

    Exact statevector p=1 QAOA for 21 qubits is computationally heavy in a
    lightweight repository script. This optional proxy preserves the same
    two-angle table interface and diagonal QUBO/surrogate objective, but it is
    not a physical QAOA expectation landscape: gamma controls softmin cost
    selectivity and beta controls mixing with the uniform distribution. Use it
    only when --fast-proxy-qaoa is explicitly requested.
    """
    E, scaled_E, _ = scaled_proxy_energies(diag_E)
    p_cost = proxy_cost_distribution(scaled_E, gamma)
    mix = math.sin(float(beta)) ** 2
    p = (1.0 - mix) * p_cost + mix / E.size
    return p / p.sum()


def evaluate_qaoa_landscape(
    diag_E: np.ndarray, grid: int, gamma_max: float = math.pi, beta_max: float = math.pi / 2, exact: bool = True
) -> pd.DataFrame:
    gammas = np.linspace(0.0, gamma_max, grid)
    betas = np.linspace(0.0, beta_max, grid)
    rows = []
    if exact:
        for gamma in gammas:
            for beta in betas:
                expE = qaoa_p1_expectation(diag_E, float(gamma), float(beta))
                rows.append(
                    {
                        "gamma": gamma,
                        "beta": beta,
                        "expected_surrogate_energy": expE,
                        "mode": "exact_statevector_p1",
                    }
                )
    else:
        E, scaled_E, uniform_expectation = scaled_proxy_energies(diag_E)
        for gamma in gammas:
            p_cost = proxy_cost_distribution(scaled_E, float(gamma))
            cost_expectation = float(np.dot(p_cost, E))
            for beta in betas:
                mix = math.sin(float(beta)) ** 2
                expE = (1.0 - mix) * cost_expectation + mix * uniform_expectation
                rows.append(
                    {
                        "gamma": gamma,
                        "beta": beta,
                        "expected_surrogate_energy": expE,
                        "mode": "fast_diagnostic_proxy_not_qaoa",
                    }
                )
    return pd.DataFrame(rows)


def sample_best_from_qaoa(diag_E: np.ndarray, gamma: float, beta: float, shots: int, seed: int = 123, exact: bool = False) -> Tuple[int, pd.DataFrame]:
    """Sample the selected exact/proxy QAOA distribution and return best observed state."""
    rng = np.random.default_rng(seed)
    nstates = diag_E.size
    if exact:
        prob = qaoa_p1_probabilities(diag_E, gamma, beta)
    else:
        prob = proxy_qaoa_probabilities(diag_E, gamma, beta)
    sampled = rng.choice(nstates, size=shots, replace=True, p=prob)
    counts = pd.Series(sampled).value_counts().rename_axis("state_index").reset_index(name="counts")
    counts["surrogate_energy"] = diag_E[counts["state_index"].to_numpy(dtype=int)]
    counts = counts.sort_values(["surrogate_energy", "counts"], ascending=[True, False]).reset_index(drop=True)
    best_state = int(counts.loc[0, "state_index"])
    return best_state, counts


def continuous_surrogate_branch(
    c: float,
    g: np.ndarray,
    H: np.ndarray,
    x_anchor: np.ndarray,
    trust_delta: float,
) -> Tuple[np.ndarray, np.ndarray, float, Dict[str, object]]:
    """Find the continuous minimum of the local surrogate in the trust region."""

    def objective(y: np.ndarray) -> float:
        return float(surrogate_energy(y, c, g, H))

    bounds = [(-float(trust_delta), float(trust_delta)) for _ in range(len(PARAM_NAMES))]
    result = minimize(
        objective,
        np.zeros(len(PARAM_NAMES), dtype=float),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 2000, "ftol": 1e-14},
    )
    y_continuous = np.asarray(result.x, dtype=float)
    params_continuous = latent_to_physical(local_displacement_to_latent(y_continuous, x_anchor))
    stats = {
        "success": bool(result.success),
        "surrogate_energy": float(result.fun),
        "nit": int(getattr(result, "nit", 0)),
    }
    return params_continuous, y_continuous, float(result.fun), stats


def state_index_to_lattice_indices(state_index: int, n_params: int, bits: int) -> np.ndarray:
    indices = []
    for p in range(n_params):
        idx = 0
        for b in range(bits):
            bit_pos = p * bits + b
            idx += ((state_index >> bit_pos) & 1) << b
        indices.append(idx)
    return np.asarray(indices, dtype=int)


def export_spectra(
    out_dir: Path,
    df: pd.DataFrame,
    classical_params: np.ndarray,
    continuous_params: np.ndarray,
    discrete_params: np.ndarray,
) -> None:
    f = df["frequency_Hz"].to_numpy(float)
    z_exp = df["Zreal_ohm"].to_numpy(float) - 1j * df["minus_Zimag_ohm"].to_numpy(float)
    z_classical = circuit_impedance(classical_params, f)
    z_continuous = circuit_impedance(continuous_params, f)
    z_discrete = circuit_impedance(discrete_params, f)
    out = pd.DataFrame(
        {
            "frequency_Hz": f,
            "Zreal_exp_ohm": z_exp.real,
            "minus_Zimag_exp_ohm": -z_exp.imag,
            "Zreal_classical_ohm": z_classical.real,
            "minus_Zimag_classical_ohm": -z_classical.imag,
            "Zreal_continuous_ohm": z_continuous.real,
            "minus_Zimag_continuous_ohm": -z_continuous.imag,
            "Zreal_discrete_ohm": z_discrete.real,
            "minus_Zimag_discrete_ohm": -z_discrete.imag,
            "absZ_exp_ohm": np.abs(z_exp),
            "absZ_classical_ohm": np.abs(z_classical),
            "absZ_continuous_ohm": np.abs(z_continuous),
            "absZ_discrete_ohm": np.abs(z_discrete),
            "phase_exp_deg": np.angle(z_exp, deg=True),
            "phase_classical_deg": np.angle(z_classical, deg=True),
            "phase_continuous_deg": np.angle(z_continuous, deg=True),
            "phase_discrete_deg": np.angle(z_discrete, deg=True),
        }
    )
    out.to_csv(out_dir / "eis_fitted_spectra_classical_discrete.csv", index=False)
    out.to_csv(out_dir / "eis_fitted_spectra_classical_continuous_discrete.csv", index=False)


def export_surrogate_slices(out_dir: Path, c: float, g: np.ndarray, H: np.ndarray, trust_delta: float) -> None:
    """Export manuscript-facing surrogate slices and standalone support files.

    The combined table `surrogate_2d_slices.csv` contains all regenerated
    two-parameter surrogate slices.  For reviewer-facing traceability, this
    function also exports the two manuscript-supporting slice files
    `R1_Q1_slice.csv` and `Rs_alpha1_slice.csv` directly from the same
    regenerated surrogate surface, so these files are no longer archived-only
    artifacts.

    Internally, the equivalent-circuit resistance is named `Rct`; the manuscript
    and legacy figure files may label the same parameter as `R1`.  Therefore,
    `R1_Q1_slice.csv` is regenerated from the internal (`Rct`, `Q1`) pair and
    exported with manuscript-facing alias labels (`R1`, `Q1`).
    """
    # (internal_x, internal_y, manuscript_x, manuscript_y, optional_standalone_file)
    pairs = [
        ("Rs", "Rct", "Rs", "Rct", None),
        ("Q1", "alpha1", "Q1", "alpha1", None),
        ("Q2", "alpha2", "Q2", "alpha2", None),
        ("Rct", "Q1", "R1", "Q1", "R1_Q1_slice.csv"),
        ("Rs", "alpha1", "Rs", "alpha1", "Rs_alpha1_slice.csv"),
    ]
    grid = np.linspace(-trust_delta, trust_delta, 81)
    all_rows = []
    standalone_rows = {}

    for internal_x, internal_y, label_x, label_y, standalone_name in pairs:
        ix, iy = PARAM_NAMES.index(internal_x), PARAM_NAMES.index(internal_y)
        rows = []
        for vx in grid:
            for vy in grid:
                y = np.zeros(len(PARAM_NAMES))
                y[ix] = vx
                y[iy] = vy
                row = {
                    "slice": f"{label_x}_vs_{label_y}",
                    "param_x": label_x,
                    "param_y": label_y,
                    "internal_param_x": internal_x,
                    "internal_param_y": internal_y,
                    "normalized_x": vx,
                    "normalized_y": vy,
                    "surrogate_energy": surrogate_energy(y, c, g, H),
                }
                rows.append(row)
                all_rows.append(row)
        if standalone_name is not None:
            standalone_rows[standalone_name] = rows

    pd.DataFrame(all_rows).to_csv(out_dir / "surrogate_2d_slices.csv", index=False)
    for standalone_name, rows in standalone_rows.items():
        pd.DataFrame(rows).to_csv(out_dir / standalone_name, index=False)


def make_plots(out_dir: Path) -> None:
    spectra_path = out_dir / "eis_fitted_spectra_classical_continuous_discrete.csv"
    if not spectra_path.exists():
        spectra_path = out_dir / "eis_fitted_spectra_classical_discrete.csv"
    spec = pd.read_csv(spectra_path)
    plt.figure(figsize=(5.2, 4.6))
    plt.plot(spec["Zreal_exp_ohm"], spec["minus_Zimag_exp_ohm"], "o", ms=4, label="Raw EIS")
    plt.plot(spec["Zreal_classical_ohm"], spec["minus_Zimag_classical_ohm"], "-", lw=2, label="Classical fit")
    if {"Zreal_continuous_ohm", "minus_Zimag_continuous_ohm"}.issubset(spec.columns):
        plt.plot(
            spec["Zreal_continuous_ohm"],
            spec["minus_Zimag_continuous_ohm"],
            "-.",
            lw=2,
            label="Continuous fit",
        )
    plt.plot(spec["Zreal_discrete_ohm"], spec["minus_Zimag_discrete_ohm"], "--", lw=2, label="Decoded discrete")
    plt.xlabel("Z' (Ω)")
    plt.ylabel("-Z'' (Ω)")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(out_dir / "nyquist_classical_discrete.png", dpi=300)
    plt.savefig(out_dir / "nyquist_classical_continuous_discrete.png", dpi=300)
    plt.close()

    landscape_path = out_dir / "qaoa_p1_landscape.csv"
    if landscape_path.exists():
        land = pd.read_csv(landscape_path)
        value_col = "expected_surrogate_energy" if "expected_surrogate_energy" in land.columns else "qaoa_expected_surrogate_energy"
        piv = land.pivot(index="beta", columns="gamma", values=value_col)
        plt.figure(figsize=(5.4, 4.6))
        plt.imshow(
            piv.values,
            origin="lower",
            aspect="auto",
            extent=[land["gamma"].min(), land["gamma"].max(), land["beta"].min(), land["beta"].max()],
        )
        mode = str(land["mode"].iloc[0]) if "mode" in land.columns and len(land) else "qaoa"
        mode_label = {
            "exact_statevector_p1": "exact p=1 statevector QAOA",
            "fast_diagnostic_proxy_not_qaoa": "fast diagnostic proxy, not QAOA",
        }.get(mode, mode.replace("_", " "))
        plt.xlabel("gamma")
        plt.ylabel("beta")
        plt.title(mode_label)
        plt.colorbar(label="QAOA expected surrogate energy")
        plt.tight_layout()
        plt.savefig(out_dir / "qaoa_p1_landscape.png", dpi=300)
        plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Shared-objective h-BN EIS classical + surrogate + QAOA pipeline")
    parser.add_argument("--input", required=True, type=Path, help="Raw EIS CSV file")
    parser.add_argument("--output", required=True, type=Path, help="Output directory")
    parser.add_argument("--trust-delta", type=positive_float, default=0.08, help="Trust-region half-width in latent coordinates")
    parser.add_argument("--fd-step", type=positive_float, default=0.015, help="Finite-difference step for surrogate Hessian")
    parser.add_argument("--bits-per-parameter", type=bits_per_parameter, default=3, help="Binary discretization per parameter, from 1 to 3 in this release script")
    parser.add_argument("--qaoa-grid", type=positive_int, default=9, help="Number of gamma and beta grid points")
    parser.add_argument("--shots", type=positive_int, default=4096, help="Shot count for sampling the best QAOA bitstring")
    parser.add_argument("--skip-qaoa", action="store_true", help="Skip QAOA landscape and use the exact lattice minimum only")
    parser.add_argument("--exact-statevector-qaoa", action="store_true", help="Compatibility flag; exact p=1 statevector QAOA is now the default.")
    parser.add_argument("--fast-proxy-qaoa", action="store_true", help="Use the old fast diagnostic proxy instead of exact p=1 QAOA. The proxy is smooth and is not a physical QAOA expectation landscape.")
    parser.add_argument("--seed", type=int, default=123, help="Random seed for QAOA shot sampling")
    args = parser.parse_args()

    if args.exact_statevector_qaoa and args.fast_proxy_qaoa:
        parser.error("--exact-statevector-qaoa and --fast-proxy-qaoa cannot be used together")

    if not args.input.is_file():
        parser.error(f"--input does not exist or is not a file: {args.input}")

    out_dir = args.output
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_eis_csv(args.input)
    df.to_csv(out_dir / "eis_raw_standardized.csv", index=False)
    freq = df["frequency_Hz"].to_numpy(float)
    z_exp = df["Zreal_ohm"].to_numpy(float) - 1j * df["minus_Zimag_ohm"].to_numpy(float)

    classical_params, fit_stats = fit_classical_anchor(freq, z_exp)
    x_anchor = physical_to_latent(classical_params)

    c, g, H = build_quadratic_surrogate(x_anchor, freq, z_exp, args.fd_step)

    # Export classical parameters and surrogate coefficients.
    param_df = pd.DataFrame({"parameter": PARAM_NAMES, "classical_anchor_value": classical_params, "latent_anchor": x_anchor})
    param_df.to_csv(out_dir / "classical_anchor_parameters.csv", index=False)
    pd.DataFrame(H, index=PARAM_NAMES, columns=PARAM_NAMES).to_csv(out_dir / "surrogate_hessian.csv")
    pd.DataFrame({"parameter": PARAM_NAMES, "linear_gradient": g}).to_csv(out_dir / "surrogate_gradient.csv", index=False)
    with open(out_dir / "fit_and_surrogate_metadata.json", "w", encoding="utf-8") as fh:
        json.dump(
            {
                "input_file": str(args.input),
                "circuit": "Rs + j*w*L + (Rct || CPE1) + CPE2",
                "parameter_order": PARAM_NAMES,
                "fit_stats": fit_stats,
                "surrogate_constant": c,
                "fd_step": args.fd_step,
                "trust_delta": args.trust_delta,
                "bits_per_parameter": args.bits_per_parameter,
                "shots": args.shots,
            },
            fh,
            indent=2,
        )

    levels = levels_for_bits(args.bits_per_parameter, args.trust_delta)
    continuous_params, continuous_y, continuous_surrogate_E, continuous_stats = continuous_surrogate_branch(
        c,
        g,
        H,
        x_anchor,
        args.trust_delta,
    )
    best_indices, best_surrogate_E, best_linear_index = enumerate_lattice_minimum(c, g, H, levels)
    y_lattice = levels[best_indices]
    lattice_params = latent_to_physical(local_displacement_to_latent(y_lattice, x_anchor))

    # Optional p=1 QAOA landscape.
    qaoa_best_state = None
    qaoa_counts = None
    if not args.skip_qaoa:
        diag_E = diagonal_surrogate_energies(c, g, H, levels, args.bits_per_parameter)
        use_exact_qaoa = not args.fast_proxy_qaoa
        landscape = evaluate_qaoa_landscape(diag_E, args.qaoa_grid, exact=use_exact_qaoa)
        landscape.to_csv(out_dir / "qaoa_p1_landscape.csv", index=False)
        k = int(landscape["expected_surrogate_energy"].idxmin())
        gamma_best = float(landscape.loc[k, "gamma"])
        beta_best = float(landscape.loc[k, "beta"])
        qaoa_best_state, qaoa_counts = sample_best_from_qaoa(
            diag_E, gamma_best, beta_best, args.shots, args.seed, exact=use_exact_qaoa
        )
        qaoa_counts.to_csv(out_dir / "qaoa_sampled_bitstrings.csv", index=False)
        qaoa_indices = state_index_to_lattice_indices(qaoa_best_state, len(PARAM_NAMES), args.bits_per_parameter)
        y_qaoa = levels[qaoa_indices]
        qaoa_params = latent_to_physical(local_displacement_to_latent(y_qaoa, x_anchor))
        qaoa_surrogate_E = float(diag_E[qaoa_best_state])
        discrete_source = "p1_qaoa_sampled_best" if use_exact_qaoa else "diagnostic_proxy_sampled_best"
        discrete_indices = qaoa_indices
        discrete_y = y_qaoa
        discrete_params = qaoa_params
        discrete_surrogate_E = qaoa_surrogate_E
        qaoa_summary = {
            "gamma_best": gamma_best,
            "beta_best": beta_best,
            "qaoa_best_state_index": int(qaoa_best_state),
            "qaoa_landscape_mode": "exact_statevector_p1" if use_exact_qaoa else "fast_diagnostic_proxy_not_qaoa",
        }
    else:
        discrete_source = "exact_lattice_minimum_skip_qaoa"
        discrete_indices = best_indices
        discrete_y = y_lattice
        discrete_params = lattice_params
        discrete_surrogate_E = best_surrogate_E
        qaoa_summary = {}

    true_classical_loss = loss_from_latent(x_anchor, freq, z_exp)
    true_continuous_loss = loss_from_latent(physical_to_latent(continuous_params), freq, z_exp)
    true_discrete_loss = loss_from_latent(physical_to_latent(discrete_params), freq, z_exp)

    decoded_df = pd.DataFrame(
        {
            "parameter": PARAM_NAMES,
            "classical_anchor": classical_params,
            "continuous_surrogate_branch": continuous_params,
            "local_y_continuous": continuous_y,
            "decoded_discrete": discrete_params,
            "local_y_discrete": discrete_y,
            "lattice_level_index": discrete_indices,
        }
    )
    decoded_df.to_csv(out_dir / "decoded_discrete_parameters.csv", index=False)

    comparison = pd.DataFrame(
        [
            {
                "branch": "classical_anchor",
                "unweighted_true_SSE": true_classical_loss,
                "relative_weighted_SSE": loss_from_latent(x_anchor, freq, z_exp, weighted=True),
                "surrogate_energy": c,
            },
            {
                "branch": discrete_source,
                "unweighted_true_SSE": true_discrete_loss,
                "relative_weighted_SSE": loss_from_latent(physical_to_latent(discrete_params), freq, z_exp, weighted=True),
                "surrogate_energy": discrete_surrogate_E,
            },
            {
                "branch": "continuous_surrogate_minimum",
                "unweighted_true_SSE": true_continuous_loss,
                "relative_weighted_SSE": loss_from_latent(physical_to_latent(continuous_params), freq, z_exp, weighted=True),
                "surrogate_energy": continuous_surrogate_E,
            },
            {
                "branch": "exact_lattice_minimum",
                "unweighted_true_SSE": loss_from_latent(physical_to_latent(lattice_params), freq, z_exp),
                "relative_weighted_SSE": loss_from_latent(physical_to_latent(lattice_params), freq, z_exp, weighted=True),
                "surrogate_energy": best_surrogate_E,
            },
        ]
    )
    comparison.to_csv(out_dir / "branch_loss_comparison.csv", index=False)

    bitstring = lattice_indices_to_bitstring(discrete_indices, args.bits_per_parameter)
    with open(out_dir / "decoded_bitstring_summary.json", "w", encoding="utf-8") as fh:
        json.dump(
            {
                "discrete_source": discrete_source,
                "parameter_order": PARAM_NAMES,
                "bits_per_parameter": args.bits_per_parameter,
                "bitstring_by_parameter_msb_first": bitstring,
                "lattice_level_indices": [int(x) for x in discrete_indices],
                "qaoa_summary": qaoa_summary,
                "continuous_branch": {
                    "source": "continuous_surrogate_minimum",
                    "local_y_continuous": [float(x) for x in continuous_y],
                    "optimization_stats": continuous_stats,
                },
                "loss_definition": "unweighted complex SSE over real and imaginary residuals",
                "true_classical_loss": true_classical_loss,
                "true_continuous_loss": true_continuous_loss,
                "true_discrete_loss": true_discrete_loss,
                "surrogate_continuous_energy": continuous_surrogate_E,
                "surrogate_discrete_energy": discrete_surrogate_E,
            },
            fh,
            indent=2,
        )

    export_spectra(out_dir, df, classical_params, continuous_params, discrete_params)
    export_surrogate_slices(out_dir, c, g, H, args.trust_delta)
    make_plots(out_dir)

    print("Shared-objective EIS pipeline completed.")
    print(f"Output directory: {out_dir}")
    print(f"Classical unweighted SSE: {true_classical_loss:.6g}")
    print(f"Discrete unweighted SSE:  {true_discrete_loss:.6g}")
    print(f"Discrete bitstring: {bitstring}")


if __name__ == "__main__":
    main()
