import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Read the CSV file
file_path = 'simulation_results.csv'  # Replace with the actual file path
df = pd.read_csv(file_path)

# Calculate mean and standard deviation for each transition
mean_transitions = df.groupby('Transition')['Transitions'].mean()
std_transitions = df.groupby('Transition')['Transitions'].std()

# Plot the results
plt.figure(figsize=(12, 8))
sns.barplot(x=mean_transitions.index, y=mean_transitions, yerr=std_transitions, capsize=5)
plt.title('Mean Transitions with Error Bars')
plt.xlabel('Transition')
plt.ylabel('Mean Transitions')
plt.show()
