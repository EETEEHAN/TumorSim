import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import sys

# Import TumorSimV8 module
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from TumorSimV8 import Environment, Tumor, Cancer_Cell

# --- Step 1: Load real tumor growth data ---
real_df = pd.read_excel("tumorgrowth.xlsx", sheet_name="AnalysisData")
real_df = real_df.groupby("Day")["Size"].mean().reset_index()
real_df["TumorSize"] = real_df["Size"] / real_df["Size"].max()
real_df["Day"] = real_df["Day"] / real_df["Day"].max()
real_df["Source"] = "Real Data"

# --- Step 2: Run calibrated simulation ---
def run_calibrated_simulation(mutation_rate, proliferation_chance, aggressiveness, steps=100):
    env = Environment(width=20, height=20)
    env.initialize_grid()
    tumor = Tumor(env)
    tumor.seed_initial_cancer(Cancer_Cell(
        position=(10, 10),
        mutation_rate=mutation_rate,
        proliferation_chance=proliferation_chance,
        aggressiveness=aggressiveness
    ))

    for _ in range(steps):
        tumor.step()

    sim_df = pd.DataFrame(tumor.history)
    sim_df["TumorSize"] = sim_df["cancer_cell_count"] / sim_df["cancer_cell_count"].max()
    sim_df["Day"] = np.linspace(0, 1, len(sim_df))
    sim_df["Source"] = "Simulated"
    return sim_df[["Day", "TumorSize", "Source"]]

# Use best-fit parameters
sim_df = run_calibrated_simulation(mutation_rate=0.01, proliferation_chance=0.3, aggressiveness=1.2)

# --- Step 3: Combine and plot ---
combined_df = pd.concat([real_df[["Day", "TumorSize", "Source"]], sim_df])

plt.figure(figsize=(10, 6))
sns.lineplot(data=combined_df, x="Day", y="TumorSize", hue="Source", marker="o")
plt.title("Tumor Growth: Calibrated Simulation vs Real Data (Normalized)")
plt.xlabel("Normalized Time")
plt.ylabel("Normalized Tumor Size")
plt.tight_layout()
plt.savefig("tumor_fit_vs_real.png")
plt.show()

