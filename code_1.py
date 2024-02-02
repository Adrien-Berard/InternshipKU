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

    def add_polymerases(self, count, existing_polymerase_positions):
        # Add a specified number of new polymerases at non-overlapping positions
        for _ in range(count):
            new_position = 0
            while new_position in existing_polymerase_positions:
                new_position = random.randint(0, len(self.histones) - 1)
            existing_polymerase_positions.append(new_position)

class Polymerase:
    def __init__(self, chromatine, position=0, temperature=1.0):
        # Initialize polymerase with a reference to the chromatine, a starting position, and a temperature for movement
        self.chromatine = chromatine
        self.position = position
        self.temperature = temperature

    def calculate_attractive_energy(self, chromatine):
        # Calculate attractive energy based on histones in front and back of the polymerase
        attractive_force = 10  # You may adjust this as needed
        max_distance = 20  # Maximum distance to consider (20 points in front and 20 points behind)

        # Calculate energy for histones in front of the polymerase
        front_energy = sum([attractive_force / max(1, abs(self.position - i)) for i in range(self.position + 1, self.position + 1 + max_distance)])

        # Calculate energy for histones behind the polymerase
        back_energy = sum([attractive_force / max(1, abs(self.position - i)) for i in range(self.position - 1, self.position - 1 - max_distance, -1)])

        # Max to avoid dividing by zero and physically consider the histone on the current position as inactive
        return front_energy, back_energy

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
        # Simulate the histone cleavage process by polymerase
        if chromatine.histones[self.position] in {'H', 'N'}:
            chromatine.histones[self.position] = ' '

def visualize_chromatine(histones, polymerase_positions=None):
    # Display chromatine state as a bar chart
    plt.bar(range(len(histones)), [1] * len(histones),
            color=[('gray' if hist == ' ' else 'blue' if hist == 'H' else 'green' if hist == 'N' else 'red') for hist in histones])
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
        polymerase.cleave_histone(chromatine)
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
    active_histone_count = sum(1 for histone in chromatine.histones if histone in {'H', 'N'})
    active_histone_count_over_time.append(active_histone_count)

    # Clear the previous frame after updating the data
    axs[0, 0].clear()
    axs[0, 1].clear()
    axs[1, 1].clear()

    axs[0, 0].plot(range(1, frame + 2), polymerase_count_over_time[:frame + 1], marker='o')
    axs[0, 0].set_title('Number of polymerases Over Time')
    axs[0, 0].set_xlabel('Time Steps')
    axs[0, 0].set_ylabel('Number of polymerases')

    axs[0, 1].plot(range(1, frame + 2), active_histone_count_over_time[:frame + 1], marker='o', color='green')
    axs[0, 1].set_title('Number of Active Histones Over Time')
    axs[0, 1].set_xlabel('Time Steps')
    axs[0, 1].set_ylabel('Number of Active Histones')

    # Visualize chromatine structure with arrows indicating polymerase positions
    visualize_chromatine(chromatine.histones, polymerase_positions=polymerase_positions)


# Parameters for simulation
chromatine_size = 50
polymerase_count = 2
simulation_steps = 100

# Initialize chromatine and polymerases with a specified temperature
chromatine = Chromatine(chromatine_size)
polymerases = [Polymerase(chromatine, temperature=1.0) for _ in range(polymerase_count)]

# Track existing polymerase positions using a list to avoid duplicates
existing_polymerase_positions = [polymerase.position for polymerase in polymerases]

# Lists to store the number of polymerases and active histones over time
polymerase_count_over_time = []
active_histone_count_over_time = []

# Create an animated plot
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()
