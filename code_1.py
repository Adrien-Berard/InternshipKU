import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

class Chromatine:
    def __init__(self, histones_count):
        # Initialize chromatine with histones
        self.histones = ['H' for _ in range(histones_count)]

    def regenerate_histones(self, positions):
        # Replace histones at specified positions with new histones ('N')
        for pos in positions:
            self.histones[pos] = 'N'

    def add_rnases(self, count, existing_rnase_positions):
        # Add a specified number of new RNases at the beginning without overlapping existing RNase positions
        for _ in range(count):
            new_position = 0
            while new_position in existing_rnase_positions:
                new_position += 1
            existing_rnase_positions.append(new_position)

class RNase:
    def __init__(self, chromatine, position=0, temperature=1.0):
        # Initialize RNase with a reference to the chromatine, a starting position, and a temperature for movement
        self.chromatine = chromatine
        self.position = position
        self.temperature = temperature

    def calculate_attractive_energy(self, chromatine):
        # Calculate attractive energy based on histones in front and back of the RNase
        attractive_force = -10  # You may adjust this as needed
        max_distance = 20  # Maximum distance to consider (2 points in front and 2 points behind)

        # Calculate energy for histones in front of the RNase
        front_energy = sum([attractive_force / max(1, abs(self.position - i)) for i in range(self.position + 1, self.position + 1 + max_distance)])

        # Calculate energy for histones behind the RNase
        back_energy = sum([attractive_force / max(1, abs(self.position - i)) for i in range(self.position - 1, self.position - 1 - max_distance, -1)])

        # Max to avoid dividing by zero and physically consider the histone on the current position as inactive
        return front_energy, back_energy

    def calculate_boltzmann_probability(self, energy_difference):
        # Calculate Boltzmann probability
        k_B = 1.0  # Boltzmann constant (you can adjust this as needed)
        return math.exp(-energy_difference / (k_B * self.temperature))

    def move(self, chromatine):
        # Define two possible states for movement (left and right)
        states = [-1, 1]

        # Calculate energies for the current and potential next positions
        current_energy = 0
        next_positions = [self.position + state for state in states]
        next_energies = list(self.calculate_attractive_energy(chromatine))

        # Calculate probabilities using Boltzmann distribution
        temperature = self.temperature
        probabilities = [math.exp(-energy / temperature) for energy in next_energies]

        # Normalize probabilities to sum to 1
        total_prob = sum(probabilities)
        normalized_probabilities = [prob / total_prob for prob in probabilities]

        # Choose the next position based on normalized probabilities
        self.position = random.choices(next_positions, weights=normalized_probabilities, k=1)[0]

        # Boundering conditions
        self.position = max(0, min(self.position, len(chromatine.histones) - 1))

    def cleave_histone(self, chromatine):
        # Simulate the histone cleavage process by RNase
        if chromatine.histones[self.position] in {'H', 'N'}:
            chromatine.histones[self.position] = ' '

def visualize_chromatine(histones, rnase_positions=None):
    # Display chromatine state as a bar chart
    plt.bar(range(len(histones)), [1] * len(histones),
            color=[('gray' if hist == ' ' else 'blue' if hist == 'H' else 'green' if hist == 'N' else 'red') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")

    if rnase_positions:
        # Draw arrows indicating RNase positions
        for position in rnase_positions:
            plt.arrow(position, -0.2, 0, -0.1, color='red', head_width=0.5, head_length=0.1)

# Function to update the plot in each animation frame
def update(frame):
    deleted_positions = []  # Keep track of deleted positions for regeneration
    for rnase in rnases:
        rnase.move(chromatine)
        rnase.cleave_histone(chromatine)
        deleted_positions.append(rnase.position)

    # Randomly add new RNase at the beginning of the chromatine with a certain probability
    if random.random() < 0.01:  # Adjust the probability as needed
        # Add new RNases with random positions
        new_rnase_positions = random.sample(range(1, len(chromatine.histones)), k=1)  # Avoid position 0 to prevent overlap
        chromatine.add_rnases(1, new_rnase_positions)
        existing_rnase_positions.extend(new_rnase_positions)
        rnases.append(RNase(chromatine, temperature=1.0))

    # Regenerate histones at deleted positions
    if random.random() < 0.1:
        chromatine.regenerate_histones(deleted_positions)

    # Update the number of RNases and active histones lists
    rnase_count_over_time.append(len(rnases))
    active_histone_count = sum(1 for histone in chromatine.histones if histone in {'H', 'N'})
    active_histone_count_over_time.append(active_histone_count)

    # Clear the previous frame after updating the data
    axs[0, 0].clear()
    axs[0, 1].clear()

    # Plot the number of RNases and active histones over time
    axs[0, 0].clear()
    axs[0, 1].clear()

    axs[0, 0].plot(range(1, frame + 2), rnase_count_over_time[:frame + 1], marker='o')
    axs[0, 0].set_title('Number of RNases Over Time')
    axs[0, 0].set_xlabel('Time Steps')
    axs[0, 0].set_ylabel('Number of RNases')

    axs[0, 1].plot(range(1, frame + 2), active_histone_count_over_time[:frame + 1], marker='o', color='green')
    axs[0, 1].set_title('Number of Active Histones Over Time')
    axs[0, 1].set_xlabel('Time Steps')
    axs[0, 1].set_ylabel('Number of Active Histones')

    # Visualize chromatine structure with arrows indicating RNase positions
    visualize_chromatine(chromatine.histones, rnase_positions=[rnase.position for rnase in rnases])

# Parameters for simulation
chromatine_size = 50
rnase_count = 2
simulation_steps = 100

# Initialize chromatine and RNases with a specified temperature
chromatine = Chromatine(chromatine_size)
rnases = [RNase(chromatine, temperature=1.0) for _ in range(rnase_count)]

# Track existing RNase positions using a list to avoid duplicates
existing_rnase_positions = [rnase.position for rnase in rnases]

# Lists to store the number of RNases and active histones over time
rnase_count_over_time = []
active_histone_count_over_time = []

# Create an animated plot
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()
