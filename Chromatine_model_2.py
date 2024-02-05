import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Chromatine:
    def __init__(self, histones_count):
        # Initialize chromatine with histones
        self.histones = np.full(histones_count, 'M', dtype='U1')

        # Randomly set approximately 10% of histones to 'U' (unmodified)
        num_unmodified = int(0.1 * histones_count)
        unmodified_positions = np.random.choice(histones_count, size=num_unmodified, replace=False)
        self.histones[unmodified_positions] = 'U'

    def regenerate_histones(self, positions):
        # Replace histones at specified unmodified positions with new histones ('A')
        unmodified_positions = np.where(self.histones == 'U')[0]
        selected_positions = np.random.choice(unmodified_positions, size=len(positions), replace=False)
        self.histones[selected_positions] = 'A'

    def add_polymerases(self, count, existing_polymerase_positions, adding_position):
        # Add a specified number of new polymerases at non-overlapping positions
        for _ in range(count):
            new_position = adding_position
            while new_position in existing_polymerase_positions:
                new_position += 1
            existing_polymerase_positions.append(new_position)

    def adding_poly_proba(self, adding_position):
        # Linear function of the local density of histones
        # Let's calculate local density as the count of active histones in the vicinity of the polymerase
        vicinity_size = 5  # Adjust this size as needed
        start_index = max(0, adding_position - vicinity_size)
        end_index = min(len(self.histones), adding_position + vicinity_size + 1)

        local_density = np.sum(np.isin(self.histones[start_index:end_index], ['M', 'A']))

        # Linear function parameters (you may adjust these)
        slope = 1e-4  # Adjust the slope to control the influence of local density
        intercept = 1e-1 # Adjust the intercept to control the baseline probability

        # Calculate the probability of adding a new polymerase
        probability = slope * local_density + intercept

        return probability

    def change_next_histones(self, position):
        # Simulate the influence of neighbors on the next histones
        if 0 <= position < len(self.histones) - 1:
            current_histone = self.histones[position]
            next_histone = self.histones[position + 1]

            # Transition U A  -> U U 
            if current_histone == 'U' and next_histone == 'A':
                # If the next histone is Acetylated, change the next histone to Unmodified
                self.histones[position + 1] = 'U'

            # Transition U M -> U U 
            elif current_histone == 'U' and next_histone == 'M':
                # If the next histone is Methylated, change the next histone to Unmodified
                self.histones[position + 1] = 'U'

            # Transition A U  -> A A 
            elif current_histone == 'A' and next_histone == 'U':
                # If the next histone is Unmodified, change the next histone to Acetylated
                self.histones[position + 1] = 'A'

            # Transition A M  -> A A 
            elif current_histone == 'A' and next_histone == 'M':
                # If the next histone is Methylated, change the next histone to Acetylated
                self.histones[position + 1] = 'A'

            # Transition M A  -> M M 
            elif current_histone == 'M' and next_histone == 'A':
                # If the next histone is Acetylated, change the next histone to Methylated
                self.histones[position + 1] = 'M'

            # Transition M U  -> M M 
            elif current_histone == 'M' and next_histone == 'U':
                # If the next histone is Unmodified, change the next histone to Methylated
                self.histones[position + 1] = 'M'

class Polymerase:
    def __init__(self, chromatine, position=10, temperature=1.0):
        # Initialize polymerase with a reference to the chromatine, a starting position, and a temperature for movement
        self.chromatine = chromatine
        self.position = position
        self.temperature = temperature

    def delete(self):
        polymerases.remove(self)

    def move(self, chromatine):
        # Define two possible states for movement (left and right)
        states = [0, 1]

        next_positions = [self.position + state for state in states]

        probabilities = [1/2, 30]

        # Normalize probabilities to sum to 1
        total_prob = np.sum(probabilities)
        normalized_probabilities = probabilities / total_prob

        # Choose the next position based on normalized probabilities
        self.position = np.random.choice(next_positions, p=normalized_probabilities)

        # Bounding conditions
        if self.position == len(chromatine.histones):
            self.delete()

    def change_histones(self, chromatine):
        # Simulate the histone change process by polymerase
        if 0 <= self.position < len(chromatine.histones) and chromatine.histones[self.position] == 'M':
            chromatine.histones[self.position] = 'A'

def visualize_chromatine(histones, polymerase_positions=None):
    # Display chromatine state as a bar chart
    plt.bar(range(len(histones)), [1] * len(histones),
            color=[('gray' if hist == 'U' else 'blue' if hist == 'M' else 'green' if hist == 'A' else 'red') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")

    if polymerase_positions:
        # Draw arrows indicating polymerase positions
        for position in polymerase_positions:
            plt.arrow(position, -0.2, 0, -0.1, color='red', head_width=0.5, head_length=0.1)

# Function to update the plot in each animation frame
def update(frame):
    deleted_positions = []  # Keep track of deleted positions for regeneration
    polymerase_positions = []  # Clear the polymerase_positions list

    for polymerase in polymerases:
        polymerase.move(chromatine)
        polymerase.change_histones(chromatine)
        deleted_positions.append(polymerase.position)
        polymerase_positions.append(polymerase.position)  # Append the current position

    # Randomly add new polymerase at the beginning of the chromatine with a certain probability
    if np.random.random() < chromatine.adding_poly_proba(adding_position):  # Adjust the probability as needed
        # Add new polymerases with non-overlapping random positions
        chromatine.add_polymerases(1, existing_polymerase_positions, adding_position)
        new_polymerase_positions = existing_polymerase_positions[-1:]
        new_polymerases = [Polymerase(chromatine, position=pos, temperature=1.0) for pos in new_polymerase_positions]
        polymerases.extend(new_polymerases)

    # Regenerate histones at unmodified positions
    if np.random.random() < 0.4:
        chromatine.regenerate_histones(deleted_positions)

    # Change the next histones based on the influence of first neighbors
    for position in existing_polymerase_positions:
        chromatine.change_next_histones(position)

    # Update the number of polymerases and active histones lists
    polymerase_count_over_time.append(len(polymerases))
    active_histone_count = np.sum(np.isin(chromatine.histones, ['M', 'A']))
    acetylated_histone_count = np.sum(chromatine.histones == 'A')
    methylated_histone_count = np.sum(chromatine.histones == 'M')
    unmodified_histone_count = np.sum(chromatine.histones == 'U')
    active_histone_count_over_time.append(active_histone_count)
    unmodified_histone_count_over_time.append(unmodified_histone_count)
    acetylated_histone_count_over_time.append(acetylated_histone_count)
    methylated_histone_count_over_time.append(methylated_histone_count)
    

    # Only extend the lists if the frame is greater than the current length of the lists
    if frame + 1 > len(polymerase_count_over_time):
        polymerase_count_over_time.append(0)
        acetylated_histone_count_over_time.append(0)
        methylated_histone_count_over_time.append(0)
        active_histone_count_over_time.append(0)

    # Clear the previous frame after updating the data
    axs[0, 0].clear()
    axs[0, 1].clear()
    axs[1, 1].clear()

    axs[0, 0].plot(range(1, frame + 2), polymerase_count_over_time[:frame + 1], marker='o')
    axs[0, 0].set_title('Number of polymerases Over Time')
    axs[0, 0].set_xlabel('Time Steps')
    axs[0, 0].set_ylabel('Number of polymerases')

    axs[0, 1].plot(range(1, frame + 2), active_histone_count_over_time[:frame + 1], marker='o', color='red', label='Active Histones')
    axs[0, 1].plot(range(1, frame + 2), acetylated_histone_count_over_time[:frame + 1], marker='o', color='green', label="Acetylated Histones")
    axs[0, 1].plot(range(1, frame + 2), methylated_histone_count_over_time[:frame + 1], marker='o', color='blue', label="Methylated Histones")
    axs[0, 1].plot(range(1, frame + 2), unmodified_histone_count_over_time[:frame + 1], marker='o', color='gray', label="Unmodified Histones")
    axs[0, 1].set_title('Number of Active Histones Over Time')
    axs[0, 1].set_xlabel('Time Steps')
    axs[0, 1].set_ylabel('Number of Histones')
    axs[0, 1].legend()

    # Visualize chromatine structure with arrows indicating polymerase positions
    visualize_chromatine(chromatine.histones, polymerase_positions=polymerase_positions)

# Parameters for simulation
chromatine_size = 50
polymerase_count = 0
simulation_steps = 100
adding_position = 25

# Initialize chromatine and polymerases with a specified temperature
chromatine = Chromatine(chromatine_size)
polymerases = [Polymerase(chromatine, temperature=1.0) for _ in range(polymerase_count)]

# Track existing polymerase positions using a list to avoid duplicates
existing_polymerase_positions = [polymerase.position for polymerase in polymerases]

# Lists to store the number of polymerases and active histones over time
polymerase_count_over_time = []
active_histone_count_over_time = []
acetylated_histone_count_over_time = []
methylated_histone_count_over_time = []
unmodified_histone_count_over_time = []

# Create an animated plot
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()
