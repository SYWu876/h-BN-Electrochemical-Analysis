
from pathlib import Path
import subprocess, sys

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts" / "xps"
scripts = [
    SCRIPTS / "01_make_table_x.py",
    SCRIPTS / "02_profile_fit_ns15_abcd.py",
    SCRIPTS / "03_ml_analysis_ns15_ef.py",
]
for script in scripts:
    print(f"Running {script.name} ...")
    subprocess.run([sys.executable, str(script)], check=True)
print("XPS reproducibility pipeline completed.")
