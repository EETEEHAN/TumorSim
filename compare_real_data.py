import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# --- Load real data ---
real_df = pd.read_excel("tumorgrowth.xlsx", sheet_name="AnalysisData")
real_df = real_df.groupby("Day")["Size"].mean().reset_index()
real_df.rename(columns={"Size": "TumorSize"}, inplace=True)
real_df["Source"] = "Real Data"

# --- Load simulated data ---
sim_dir = "timeseries_output"
sim_all = []

for file in os.listdir(sim_dir):
    if file.endswith(".csv"):
        df = pd.read_csv(os.path.join(sim_dir, file))
        sim_all.append(df)

sim_df = pd.concat(sim_all, ignore_index=True)

# --- Simulated tumor summary ---
sim_summary = sim_df.groupby("step")["cancer_cell_count"].mean().reset_index()
sim_summary.rename(columns={"step": "Step", "cancer_cell_count": "TumorSize"}, inplace=True)
sim_summary["Source"] = "Simulated"

# --- Normalize tumor size ---
real_df["TumorSize"] = real_df["TumorSize"] / real_df["TumorSize"].max()
sim_summary["TumorSize"] = sim_summary["TumorSize"] / sim_summary["TumorSize"].max()

# --- Normalize time (Step → Day) ---
real_df["TimeNorm"] = real_df["Day"] / real_df["Day"].max()
sim_summary["TimeNorm"] = sim_summary["Step"] / sim_summary["Step"].max()

# --- Combine and plot ---
combined = pd.concat([real_df[["TimeNorm", "TumorSize", "Source"]],
                      sim_summary[["TimeNorm", "TumorSize", "Source"]]])

plt.figure(figsize=(10, 6))
sns.lineplot(data=combined, x="TimeNorm", y="TumorSize", hue="Source", marker="o")
plt.title("Tumor Growth: Real Data vs Simulation (Normalized Time and Size)")
plt.xlabel("Normalized Time")
plt.ylabel("Normalized Tumor Size")
plt.tight_layout()
plt.savefig("tumor_growth_comparison_normalized_time.png")
plt.show()
