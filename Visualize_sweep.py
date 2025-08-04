import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv("parameter_sweep_results.csv")
df["Mutation Rate"] = df["Mutation Rate"].round(3)
df["Proliferation Chance"] = df["Proliferation Chance"].round(2)

sns.set(style="whitegrid")

# --- Heatmap ---
pivot = df.pivot_table(
    index="Mutation Rate",
    columns="Proliferation Chance",
    values="Avg Final Cell Count"
)

plt.figure(figsize=(8, 6))
sns.heatmap(pivot, annot=True, fmt=".0f", cmap="YlOrRd", cbar_kws={'label': 'Final Cancer Cell Count'})
plt.title("Final Tumor Size Across Parameter Sweep")
plt.xlabel("Proliferation Chance")
plt.ylabel("Mutation Rate")
plt.tight_layout()
plt.savefig("final_cell_heatmap.png")
plt.show()

# --- Line plots: Avg mutations and age ---
plt.figure(figsize=(14, 6))

# Plot 1 – Average Mutation Count
plt.subplot(1, 2, 1)
sns.lineplot(
    data=df,
    x="Proliferation Chance",
    y="Avg Mutation Count",
    hue="Mutation Rate",
    marker="o"
)
plt.title("Avg Mutations per Cell")
plt.xlabel("Proliferation Chance")
plt.ylabel("Average Mutations")

# Plot 2 – Average Age
plt.subplot(1, 2, 2)
sns.lineplot(
    data=df,
    x="Proliferation Chance",
    y="Avg Age",
    hue="Mutation Rate",
    marker="o",
    legend=False
)
plt.title("Avg Cancer Cell Age")
plt.xlabel("Proliferation Chance")
plt.ylabel("Average Age")

plt.suptitle("Tumor Characteristics Across Parameter Sweep", fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.88)
plt.savefig("tumor_metrics_summary.png")
plt.show()
