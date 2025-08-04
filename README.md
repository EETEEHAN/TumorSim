# TumorSim: Agent-Based Tumor Growth Simulator

TumorSim is an agent-based modeling framework for simulating cancer cell growth on a 2D grid. It tracks individual cells, their mutations, proliferation behavior, and interactions with neighboring cells. The simulator is designed to be biologically plausible and has been calibrated using real tumor growth data.

---

## Project Structure

```
TumorSim/
├── TumorSimV8.py              # Core simulation logic (Environment, Tumor, Cell classes)
├── calibrate.py               # Script to fit parameters using real tumor growth data
├── plot_final.py              # Generates final plot comparing calibrated simulation to real data
├── compare_real_data.py       # Earlier version of comparison script
├── Sweep.py                   # Parameter sweep script
├── analyze_time_series.py     # Analyzes and visualizes simulation time series
├── Generate_time_series.py    # Runs simulations and exports tumor growth over time
├── Visualize_sweep.py         # Plots summary metrics from parameter sweep
├── tumorgrowth.xlsx           # Real tumor growth data (Excel format)
├── timeseries_output/         # Folder with exported CSVs of time series simulations
├── parameter_sweep_results.csv
├── tumor_growth_comparison.png
├── final_cell_heatmap.png
├── tumor_growth_animation.gif
└── README.md                  # This file
```

---

## Simulation Overview

Each simulation begins with a single cancer cell on a 2D grid. At each step:

- Cells grow and accumulate mutations
- Cancer cells evaluate whether to divide based on proliferation chance and local crowding
- New daughter cells may inherit mutations or acquire new ones
- Tumor expands spatially and tracks cell-level metrics

---

## Biological Calibration

### Growth Calibration

Each simulation step is interpreted as one full cell division cycle (∼24 hours). This assumption allows tumor doubling times to be estimated and aligned with empirical data.

For example:
- A tumor progressing from 1 to ~128 cells in ~7–10 steps implies a doubling time of 3–5 days
- These growth dynamics fall within known ranges for breast, colorectal, and small cell lung cancers

**References**:  
- Spratt et al., *Radiology*, 1996  
- Friberg & Mattson, *Anticancer Research*, 1997

---

### Mutation Burden Calibration

According to TCGA and COSMIC:
- Most solid tumors carry ~10–100 mutations per megabase
- This translates to ~300–3000 coding mutations per tumor
- In our simulation, tumors with ~200 cells typically show 2–5 mutations per cell, aligning with these values

**References**:  
- Lawrence et al., *Nature*, 2013  
- COSMIC database

---

## Parameter Fitting

To identify realistic biological parameters, we used a calibration script (`calibrate.py`) that:

- Loads real tumor growth data from `tumorgrowth.xlsx`
- Simulates tumors using a parameterized agent-based model
- Interpolates both real and simulated growth curves to 100 timepoints
- Minimizes mean squared error (MSE) between normalized growth curves

**Best-Fit Parameters Identified**:
```
Mutation rate:     0.0100
Proliferation:     0.3000
Aggressiveness:    1.2000
Final MSE:         ~0.0074
```

The final comparison plot (in `tumor_fit_vs_real.png`) shows close alignment between real tumor data and the calibrated simulation.

---

## Getting Started

**Dependencies**:
- Python 3.11+
- numpy, pandas, matplotlib, seaborn, scikit-learn, openpyxl

**To run the simulation**:
```bash
python TumorSimV8.py
```

**To calibrate the parameters**:
```bash
python calibrate.py
```

**To generate the final plot**:
```bash
python plot_final.py
```

---

## Acknowledgments

This project was developed as part of a computational biology module focused on modeling tumor dynamics and parameter calibration using real-world data.

