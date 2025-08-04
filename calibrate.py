import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from sklearn.metrics import mean_squared_error
from scipy.optimize import minimize
import sys
import os

# --- Dynamically import TumorSimV8 ---
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from TumorSimV8 import Environment, Tumor, Cancer_Cell


def run_simulation(mutation_rate, proliferation, aggressiveness, steps=100):
    env = Environment(width=20, height=20)
    env.initialize_grid()
    tumor = Tumor(env)
    tumor.seed_initial_cancer(Cancer_Cell(
        position=(10, 10),
        mutation_rate=mutation_rate,
        proliferation_chance=proliferation,
        aggressiveness=aggressiveness
    ))

    for _ in range(steps):
        tumor.step()

    df = pd.DataFrame(tumor.history)
    df["TumorSize"] = df["cancer_cell_count"] / df["cancer_cell_count"].max()  # normalize
    df["Day"] = np.linspace(0, 1, len(df))  # normalized time
    return df[["Day", "TumorSize"]]


def load_real_data():
    df = pd.read_excel("tumorgrowth.xlsx", sheet_name="AnalysisData")
    df = df.groupby("Day")["Size"].mean().reset_index()
    df["TumorSize"] = df["Size"] / df["Size"].max()  # normalize size
    df["Day"] = df["Day"] / df["Day"].max()  # normalize time
    return df[["Day", "TumorSize"]]


def compute_mse(params):
    mutation_rate, proliferation, aggressiveness = params

    try:
        sim_df = run_simulation(mutation_rate, proliferation, aggressiveness)
        real_df = load_real_data()

        # Interpolate both to 100 timepoints for fair comparison
        common_x = np.linspace(0, 1, 100)
        sim_interp = interp1d(sim_df["Day"], sim_df["TumorSize"], fill_value="extrapolate")(common_x)
        real_interp = interp1d(real_df["Day"], real_df["TumorSize"], fill_value="extrapolate")(common_x)

        mse = mean_squared_error(real_interp, sim_interp)
        print(f"Tested {params} → MSE = {mse:.5f}")
        return mse

    except Exception as e:
        print(f"Error at params {params}: {e}")
        return np.inf


if __name__ == "__main__":
    initial_guess = [0.01, 0.3, 1.2]  # mutation_rate, proliferation, aggressiveness
    bounds = [(0.001, 0.05), (0.1, 0.8), (0.5, 2.5)]

    result = minimize(compute_mse, initial_guess, method="L-BFGS-B", bounds=bounds)

    print("\n--- Best Fit Parameters ---")
    print(f"Mutation Rate:    {result.x[0]:.4f}")
    print(f"Proliferation:    {result.x[1]:.4f}")
    print(f"Aggressiveness:   {result.x[2]:.4f}")
    print(f"Final MSE:        {result.fun:.6f}")

