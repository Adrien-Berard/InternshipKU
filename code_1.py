
import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

class Chromatine:
    def __init__(self, histones_count):
        self.histones = ['H' for _ in range(histones_count)]
        self.new_histones = set()

    def regenerate_histones(self, positions):
        for pos in positions:
            if pos not in self.new_histones:
                self.new_histones.add(pos)
                self.histones[pos] = 'N'

    def add_rnases(self, count):
        for _ in range(count):
            self.histones.insert(0, 'H')

class RNase:
    def __init__(self, chromatine, position=0, temperature=1.0):
        self.chromatine = chromatine
        self.position = position
        self.temperature = temperature

    def calculate_repulsive_energy(self, other_rnase):
        repulsive_force = 1.0  # You may adjust this as needed
        distance = abs(self.position - other_rnase.position)
        return repulsive_force / distance

    def calculate_boltzmann_probability(self, energy_difference):
        k_B = 1.0  # Boltzmann constant (you can adjust this as needed)
        return math.exp(-energy_difference / (k_B * self.temperature))

    def move(self, other_rnases):
        # Define two possible states for movement (left and right)
        states = [-1, 1]

        # Calculate energies for the current and potential next positions
        current_energy = 0
        next_positions = [self.position + state for state in states]
        next_energies = [self.calculate_repulsive_energy(other_rnase) for other_rnase in other_rnases]

        # Calculate probabilities using Boltzmann distribution
        probabilities = [self.calculate_boltzmann_probability(energy - current_energy) for energy in next_energies]

        # Choose the next position based on probabilities
        self.position = random.choices(next_positions, probabilities)[0]

        # Boundering conditions
        self.position = max(0, min(self.position, len(self.chromatine.histones) - 1))

    def cleave_histone(self):
        if self.chromatine.histones[self.position] in {'H', 'N'}:
            self.chromatine.histones[self.position] = ' '

def visualize_chromatine(histones, rnase_positions=None):
    plt.bar(range(len(histones)), [1] * len(histones),
            color=[('gray' if hist == ' ' else 'blue' if hist == 'H' else 'green' if hist == 'N' else 'red') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")

    if rnase_positions:
        for position in rnase_positions:
            plt.arrow(position, -0.2, 0, -0.1, color='red', head_width=0.5, head_length=0.1)

# Function to update the plot in each animation frame
def update(frame):
    deleted_positions = []  # Keep track of deleted positions for regeneration
    for i, rnase in enumerate(rnases):
        other_rnases = rnases[:i] + rnases[i+1:]  # Exclude the current RNase from interaction
        rnase.move(other_rnases)
        rnase.cleave_histone()
        deleted_positions.append(rnase.position)

    if frame % 10 == 0:
        chromatine.add_rnases(len(rnases))  # Add new RNases every 10 steps

    plt.clf()  # Clear the previous frame

    # Visualize chromatine structure with arrows indicating RNase positions
    visualize_chromatine(chromatine.histones, rnase_positions=[rnase.position for rnase in rnases])

    # Regenerate histones at deleted positions
    chromatine.regenerate_histones(deleted_positions)

# Parameters for simulation
chromatine_size = 50
rnase_count = 2
simulation_steps = 100

# Initialize chromatine and RNases with a specified temperature
chromatine = Chromatine(chromatine_size)
rnases = [RNase(chromatine, temperature=1.0) for _ in range(rnase_count)]

# Create an animated plot
fig, ax = plt.subplots()
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()
