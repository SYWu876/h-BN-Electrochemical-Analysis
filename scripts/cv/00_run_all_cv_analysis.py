
from pathlib import Path
import subprocess, sys

script_dir = Path(__file__).resolve().parent
ordered = [
    "01_cv_overlay.py",
    "02_current_heatmap.py",
    "03_peak_vs_sqrt_v.py",
    "04_b_value_profile.py",
    "05_dunn_separation.py",
    "06_ising_qkpca.py",
    "07_descriptor_summary.py",
]
for name in ordered:
    print(f"Running {name} ...")
    subprocess.run([sys.executable, str(script_dir / name)], check=True)
print("All CV analysis scripts completed.")
