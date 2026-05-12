# QC Circuit notebooks

This folder contains the interactive notebook used to regenerate the GCD and EIS quantum-circuit schematic figures for the manuscript companion archive.

## Notebook

- `generate_gcd_eis_qc_circuits.ipynb`

The notebook builds four circuit figures:

- GCD discharge-window QAOA circuit
- GCD bounded hardware-efficient ansatz circuit
- EIS bounded-decoder hardware-efficient ansatz circuit
- EIS 21-qubit QUBO/QAOA circuit for 3-bit encoding of seven EIS parameters

The notebook keeps bounded decoding, objective evaluation, and optimizer updates as classical helper steps outside the quantum gate sequence. Run the cells from top to bottom in Jupyter to regenerate the figures interactively.

Install the top-level repository requirements before opening the notebook:

```bash
python -m pip install -r requirements.txt
```
