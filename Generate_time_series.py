import pandas as pd
import itertools
import os
from TumorSimV8 import Environment, Tumor, Cell, Cancer_Cell

# Parameters
mutation_rates = [0.01, 0.015, 0.02]
proliferation_chances = [0.3, 0.35, 0.4]
n_steps = 120
grid_size = 20

# Output directory
os.makedirs("timeseries_output", exist_ok=True)

# Parameter sweep loop
for m_rate, p_chance in itertools.product(mutation_rates, proliferation_chances):
    print(f"Running: mutation_rate={m_rate}, proliferation_chance={p_chance}")
    
    env = Environment(width=grid_size, height=grid_size)
    env.initialize_grid()
    tumor = Tumor(env)
    
    # Seed the tumor with a custom Cancer_Cell
    tumor.seed_initial_cancer(
        Cancer_Cell(position=(grid_size // 2, grid_size // 2),
                    mutation_rate=m_rate,
                    proliferation_chance=p_chance)
    )

    # Simulate
    for _ in range(n_steps):
        tumor.step()

    # Save history as CSV
    df = pd.DataFrame(tumor.history)
    df["Mutation Rate"] = m_rate
    df["Proliferation Chance"] = p_chance
    filename = f"timeseries_output/timeseries_m{m_rate}_p{p_chance}.csv"
    df.to_csv(filename, index=False)
