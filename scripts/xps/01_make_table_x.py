
from pathlib import Path
from utils_xps import table_x_dataframe, ensure_dir

out_dir = Path(__file__).resolve().parents[2] / "data" / "processed" / "XPS"
ensure_dir(out_dir)
table_x_dataframe().to_csv(out_dir / "Table_X_13_peaks.csv", index=False)
print("Saved:", out_dir / "Table_X_13_peaks.csv")
