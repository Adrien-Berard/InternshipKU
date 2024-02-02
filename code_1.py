import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

class Chromatine:
    def __init__(self, histones_count):
        # Initialize chromatine with histones
        self.histones = ['M' for _ in range(histones_count)]

    def regenerate_histones(self, positions):
        # Replace histones at specified positions with new histones ('A')
        for pos in positions:
            if 0 <= pos < len(self.histones):
                self.histones[pos] = 'A'


    def add_polymerases(self, count, existing_polymerase_positions):
        # Add a specified number of new polymerases at non-overlapping positions
        for _ in range(count):
            new_position = 10
            while new_position in existing_polymerase_positions:
                new_position += 1
            existing_polymerase_positions.append(new_position)

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
        total_prob = sum(probabilities)
        normalized_probabilities = [prob / total_prob for prob in probabilities]

        # Choose the next position based on normalized probabilities
        self.position = random.choices(next_positions, weights=normalized_probabilities, k=1)[0]

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
            color=[('gray' if hist == ' ' else 'blue' if hist == 'M' else 'green' if hist == 'A' else 'red') for hist in histones])
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
    if random.random() < 0.05:  # Adjust the probability as needed
        # Add new polymerases with non-overlapping random positions
        chromatine.add_polymerases(1, existing_polymerase_positions)
        new_polymerase_positions = existing_polymerase_positions[-1:]
        new_polymerases = [Polymerase(chromatine, position=pos, temperature=1.0) for pos in new_polymerase_positions]
        polymerases.extend(new_polymerases)

    # Regenerate histones at deleted positions
    if random.random() < 0.4:
        chromatine.regenerate_histones(deleted_positions)

    # Update the number of polymerases and active histones lists
    polymerase_count_over_time.append(len(polymerases))
    active_histone_count = sum(1 for histone in chromatine.histones if histone in {'M', 'A'})
    acetylated_histone_count = sum(1 for histone in chromatine.histones if histone in {'A'})
    methylated_histone_count = sum(1 for histone in chromatine.histones if histone in {'M'})
    active_histone_count_over_time.append(active_histone_count)
    acetylated_histone_count_over_time.append(acetylated_histone_count)
    methylated_histone_count_over_time.append(methylated_histone_count)
    
    # Pad polymerase_count_over_time with zeros if needed
    polymerase_count_over_time.extend([0] * (frame + 2 - len(polymerase_count_over_time)))
    acetylated_histone_count_over_time.extend([0] * (frame + 2 - len(acetylated_histone_count_over_time)))
    methylated_histone_count_over_time.extend([0] * (frame + 2 - len(methylated_histone_count_over_time)))
    active_histone_count_over_time.extend([0] * (frame + 2 - len(active_histone_count_over_time)))

    # Clear the previous frame after updating the data
    axs[0, 0].clear()
    axs[0, 1].clear()
    axs[1, 1].clear()

    axs[0, 0].plot(range(1, frame + 2), polymerase_count_over_time[:frame + 1], marker='o')
    axs[0, 0].set_title('Number of polymerases Over Time')
    axs[0, 0].set_xlabel('Time Steps')
    axs[0, 0].set_ylabel('Number of polymerases')

    axs[0, 1].plot(range(1, frame + 2), active_histone_count_over_time[:frame + 1], marker='o', color='red', label ='Active Histones')
    axs[0, 1].plot(range(1, frame + 2), acetylated_histone_count_over_time[:frame + 1], marker='o', color='green', label = "Acetylated Histones")
    axs[0, 1].plot(range(1, frame + 2), methylated_histone_count_over_time[:frame + 1], marker='o', color='blue', label = "Methylated Histones")
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

# Create an animated plot
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()