
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts" / "gcd"

def run(script_name: str) -> None:
    cmd = [sys.executable, str(SCRIPTS / script_name)]
    print("Running", " ".join(cmd))
    subprocess.run(cmd, check=True)

def main():
    run("01_preprocess_and_select_windows.py")
    run("02_fit_and_summarize.py")
    print("All hBN GCD analysis steps completed.")

if __name__ == "__main__":
    main()
