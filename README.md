# Coarse PEEC Solver

A hybrid Python + C++ PEEC workflow for interconnect parasitic extraction and circuit-level simulation.

## What This Repository Contains

- `Integrated_PEEC_Solver/` (Python)
  - Builds PEEC model files from `Node/Branch` geometry
  - Extracts `L/P` matrices and writes `B2N/PORT`
  - Can launch the C++ solver directly
- `Quasi_Static_Solver/` (C++)
  - Quasi-static PEEC solver with MNA assembly
  - Supports frequency-domain (`freq`) and time-domain (`time`) solving

## Quick Start

1. Prepare input files under `Integrated_PEEC_Solver/input/` (currently `run.py` reads from `input/68/`).
2. Configure paths and units in `Integrated_PEEC_Solver/config.py`.
3. Run:

```powershell
cd Integrated_PEEC_Solver
python run.py
```

Results are written to `Integrated_PEEC_Solver/workspace/`.

## Documentation

Detailed usage, algorithm flow, and model customization instructions:

- `Integrated_PEEC_Solver/README.md`
