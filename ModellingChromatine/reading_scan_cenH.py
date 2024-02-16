import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Assuming you have a CSV file with the data, replace 'your_csv_file.csv' with your actual CSV file name.
csv_filename = 'ScanCenH_Proba_0.8_polymerasecount_0_F_10_newpolyproba_0.01.csv'
df = pd.read_csv(csv_filename)
length_chro = 198

# Define the number of simulation steps and frames to show
simulation_steps = 100000
frames_to_show = simulation_steps // 200

# Extract and reshape 'A in gene' column
As = df['A in gene'].to_numpy()
As_reshaped = As.reshape(frames_to_show, len(As) // frames_to_show)
meanA, stdA = As_reshaped.mean(axis=1), As_reshaped.std(axis=1)

# Create a scatter plot using seaborn
g = sns.relplot(
    data=df,
    x="Burst Frequency", y="CenHsize",
    hue=meanA, size=1,  # Use a constant size for all points
    palette='viridis',
    height=6,  # Set the height of the plot
)

# Set log scale for both x and y axes
g.set(xscale="log", yscale="log")

# Add colorbar
cbar = plt.colorbar(g.ax.collections[0])
cbar.set_label('Mean A in gene')

# Set axis labels and title
plt.xlabel('Burst Frequency (log scale)')
plt.ylabel('CenHsize (log scale)')
plt.title('Scatter Plot of Burst Frequency vs. CenHsize')

# Show the plot
plt.show()
