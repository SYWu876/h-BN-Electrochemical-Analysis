
from pathlib import Path
from utils_xps import table_x_dataframe, ensure_dir

out_dir = Path(__file__).resolve().parents[2] / "data" / "processed" / "XPS"
ensure_dir(out_dir)
table_x_dataframe().to_csv(out_dir / "XPS_13_peak_assignments.csv", index=False)
print("Saved:", out_dir / "XPS_13_peak_assignments.csv")
