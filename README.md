# TumorSim: A Biologically Calibrated Tumor Growth Simulator

TumorSim is an agent-based simulation of tumor growth written in Python. It models the spatial expansion, mutation accumulation, and proliferative behavior of cancer cells on a discrete 2D grid. This project was calibrated using real biological data from tumor growth studies to match observed in vitro expansion dynamics, making it suitable for hypothesis generation and educational demonstrations.

## Features

- Agent-based simulation of individual cancer cells
- Customizable parameters: mutation rate, proliferation chance, aggressiveness
- Spatial pressure and crowding constraints
- Mutation type tracking and visualization
- Time series logging and growth curve generation
- Parameter sweep and optimization tools
- Calibration using real tumor growth data via ML-based curve fitting

---

## Simulation Architecture

### Environment & Tumor Classes

The simulation is based on two main classes:

- `Environment`: defines the 2D grid, manages spatial occupancy, and tracks cell pressure
- `Tumor`: maintains a list of all `Cell` and `Cancer_Cell` agents, handles timestep logic, cell growth, mutation inheritance, and division

Each `Cancer_Cell` object contains parameters for:

- `mutation_rate`: probability of gaining a mutation at each timestep
- `proliferation_chance`: base probability of cell division
- `aggressiveness`: modifier that scales proliferation based on local pressure

### Growth and Division Dynamics

At each timestep:

1. All cells grow in a randomized order.
2. For cancer cells, local pressure is computed based on neighboring occupancy.
3. If a cell decides to divide (stochastically based on pressure and parameters), a daughter cell is spawned.
4. Mutation inheritance is modeled with additional random mutations.

---

## Calibration to Real Biological Data

### Growth Curve Alignment

We calibrated simulation parameters using a dataset of tumor size measurements over time, extracted from the `tumorgrowth.xlsx` file (`AnalysisData` sheet). To enable direct comparison, both simulated and real tumor sizes were:

- Normalized to [0, 1] range (to compare growth shapes, not absolute size)
- Mapped to normalized time (0 to 1), since simulation steps and experimental days may differ

This allows for a time-independent, shape-based alignment of tumor growth curves.

### Parameter Fitting via Machine Learning

We used a machine learning-inspired optimization strategy to identify parameters that best match the real tumor growth curve. Specifically:

- **Objective**: Minimize the mean squared error (MSE) between real and simulated tumor growth trajectories.
- **Approach**: Used `scipy.optimize.minimize` (L-BFGS-B algorithm) to search over:
  - Mutation rate ∈ [0.001, 0.05]
  - Proliferation chance ∈ [0.1, 0.8]
  - Aggressiveness ∈ [0.5, 2.5]
- **Data Handling**:
  - Interpolated both real and simulated growth curves to 100 evenly spaced timepoints
  - Compared trajectories using `sklearn.metrics.mean_squared_error`

Best-fit parameters obtained:
```
Mutation Rate:    0.0100
Proliferation:    0.3000
Aggressiveness:   1.2000
Final MSE:        ~0.007
```

These parameters were then used in the final model validation.

---

## Visualizations

- `tumor_growth_comparison.png`: Comparison of normalized real and simulated tumor growth over time
- `final_cell_heatmap.png`: Final spatial distribution of cells
- `final_mutation_heatmap.png`: Spatial mutation map of the tumor
- `tumor_growth_timeseries.png`: Raw time series data of simulated growth
- `tumor_growth_animation.gif`: Animated growth progression of the tumor

---

## Project Structure

```
TumorSim/
├── TumorSimV8.py                # Main simulation logic
├── calibrate.py                 # MSE-based parameter optimization
├── plot_final.py                # Final visualization script (calibrated vs. real)
├── Sweep.py                     # Parameter sweep over growth settings
├── analyze_time_series.py       # Plot and summary stats from multiple runs
├── Generate_time_series.py      # Batch generator for time series
├── compare_real_data.py         # Overlay real and simulated data
├── tumorgrowth.xlsx             # Experimental tumor size data
├── timeseries_output/           # CSVs from simulation runs
├── tumor_growth.json            # Optional export of tumor state
├── *.png, *.gif                 # Visual outputs
└── README.md                    # Project documentation (this file)
```

---

## Dependencies

- Python ≥ 3.8
- NumPy
- pandas
- matplotlib
- seaborn
- scipy
- scikit-learn

Install with:

```bash
pip install -r requirements.txt
```

---

## Running the Simulation

To run a calibrated simulation and generate plots:

```bash
python plot_final.py
```

To perform parameter fitting against real data:

```bash
python calibrate.py
```

To generate multiple time series and analyze trends:

```bash
python Generate_time_series.py
python analyze_time_series.py
```

---

## References

- Spratt et al., *Radiology*, 1996
- Friberg & Mattson, *Anticancer Research*, 1997
- Lawrence et al., *Nature*, 2013
- COSMIC: Catalogue Of Somatic Mutations In Cancer
