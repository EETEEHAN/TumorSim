import csv
import random
import numpy as np
from TumorSimV8 import Environment, Tumor, Cancer_Cell

# --- Sweep settings ---
mutation_rates = [0.005,0.01,0.02]
prolifs = [0.2,0.3,0.4]
REPEATS = 5
STEPS = 120
WIDTH = 20
HEIGHT = 20

# --- Output CSV ---
with open("parameter_sweep_results.csv", mode="w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Mutation Rate", "Proliferation Chance", "Avg Final Cell Count", "Avg Mutation Count", "Avg Age"])

    for m_rate in mutation_rates:
        for p_chance in prolifs:
            total_cells = 0
            total_mut = 0
            total_age = 0

            for i in range(REPEATS):
                random.seed(42 + i + int(m_rate * 1000) + int(p_chance * 1000))
                np.random.seed(42 + i + int(m_rate * 1000) + int(p_chance * 1000))

                env = Environment(WIDTH, HEIGHT)
                env.initialize_grid()

                tumor = Tumor(env)
                tumor.seed_initial_cancer(Cancer_Cell(
                    position=(WIDTH // 2, HEIGHT // 2),
                    mutation_rate=m_rate,
                    proliferation_chance=p_chance
                ))

                for _ in range(STEPS):
                    tumor.step()

                # Collect stats
                total_cells += len(tumor.cells)
                avg_mut = sum(cell.mutation_count for cell in tumor.cells) / len(tumor.cells)
                avg_age = sum(cell.age for cell in tumor.cells) / len(tumor.cells)
                total_mut += avg_mut
                total_age += avg_age

            # Compute averages
            avg_cells = total_cells / REPEATS
            avg_mut = total_mut / REPEATS
            avg_age = total_age / REPEATS

            writer.writerow([m_rate, p_chance, round(avg_cells), round(avg_mut, 2), round(avg_age, 2)])
