import pandas as pd
import matplotlib.pyplot as plt

# Assuming you have a CSV file with the data, replace 'your_csv_file.csv' with your actual CSV file name.
csv_filename = 'counting_lists_dataframe_polymerasecount_0_alpha_0.9871794871794872_F_77_addingpolyprobaintercept_0_addingpolyprobaslope_0.csv'
df = pd.read_csv(csv_filename)

# Function to update the plot in each animation frame
def update(frame):
    plt.clf()  # Clear the previous frame

    # Subplot for Polymerase Count
    plt.subplot(2, 1, 1)
    plt.plot(df['Time Steps'][:frame+1], df['Polymerase Count'][:frame+1],label='Polymerase Count')
    plt.title('Polymerase Count')
    plt.xlabel('Time Steps')
    plt.ylabel('Count')
    plt.legend()

    # Subplot for Histone Counts
    plt.subplot(2, 1, 2)

    plt.plot(df['Time Steps'][:frame+1], df['Acetylated Histone Count'][:frame+1], label='Acetylated Histone Count', color = "b")
    plt.plot(df['Time Steps'][:frame+1], df['Methylated Histone Count'][:frame+1], label='Methylated Histone Count', color = "r")
    plt.plot(df['Time Steps'][:frame+1], df['Unmodified Histone Count'][:frame+1], label='Unmodified Histone Count', color = "y")
    plt.title('Histone Counts')
    plt.xlabel('Time Steps')
    plt.ylabel('Count')
    plt.legend()

# Create a static plot without chromatine visualization
fig = plt.figure(figsize=(12, 8))
update(len(df) - 1)  # Use the last frame to show the final state
plt.tight_layout()
plt.show()
