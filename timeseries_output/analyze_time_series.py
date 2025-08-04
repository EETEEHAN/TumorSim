import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Step 1: Load and concatenate all CSVs
directory = "timeseries_output"
all_dfs = []

for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        df = pd.read_csv(os.path.join(directory, filename))
        all_dfs.append(df)

if not all_dfs:
    raise ValueError("No CSV files found in timeseries_output/")

full_df = pd.concat(all_dfs, ignore_index=True)

# Step 2: Clean and round parameters for easier plotting
full_df["Mutation Rate"] = full_df["Mutation Rate"].round(3)
full_df["Proliferation Chance"] = full_df["Proliferation Chance"].round(2)

# Step 3: Time Series Plot – Tumor Cell Count Over Time
plt.figure(figsize=(10, 6))
for (m_rate, p_chance), group in full_df.groupby(["Mutation Rate", "Proliferation Chance"]):
    label = f"m={m_rate}, p={p_chance}"
    plt.plot(group["step"], group["cancer_cell_count"], label=label)

plt.xlabel("Step")
plt.ylabel("Cancer Cell Count")
plt.title("Tumor Growth Over Time by Parameter Combination")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
plt.tight_layout()
plt.savefig("tumor_growth_timeseries.png")
plt.close()

# Step 4: Summary Heatmaps (final values only)

# Efficiently get final step for each parameter combo
final_steps = (
    full_df.loc[full_df.groupby(["Mutation Rate", "Proliferation Chance"])["step"].idxmax()]
    .reset_index(drop=True)
)

# Heatmap 1: Final Cell Count
pivot_cells = final_steps.pivot(
    index="Mutation Rate", 
    columns="Proliferation Chance", 
    values="cancer_cell_count"
)
plt.figure(figsize=(8, 6))
sns.heatmap(pivot_cells, annot=True, fmt=".0f", cmap="YlGnBu", cbar_kws={'label': 'Final Cell Count'})
plt.title("Final Tumor Size Across Parameters")
plt.savefig("final_cell_heatmap.png")
plt.close()

# Heatmap 2: Final Avg Mutation Count
pivot_mut = final_steps.pivot(
    index="Mutation Rate", 
    columns="Proliferation Chance", 
    values="average_mutations"
)
plt.figure(figsize=(8, 6))
sns.heatmap(pivot_mut, annot=True, fmt=".2f", cmap="OrRd", cbar_kws={'label': 'Avg Mutations'})
plt.title("Final Average Mutation Count Across Parameters")
plt.savefig("final_mutation_heatmap.png")
plt.close()

print("All plots generated and saved.")
