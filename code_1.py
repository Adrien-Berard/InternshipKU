import random
import matplotlib.pyplot as plt
from IPython.display import display, clear_output
import time

class Chromatine:
    def __init__(self, histones_count):
        # Initialize chromatine with histones
        self.histones = ['H' for _ in range(histones_count)]

class RNase:
    def __init__(self, chromatine, position=0):
        # Initialize RNase with a reference to the chromatine and a starting position
        self.chromatine = chromatine
        self.position = position

    def move(self):
        # Stochastic movement of RNase
        movement = random.choice([-1, 1])
        self.position += movement

        # Boundering conditions
        self.position = max(0, min(self.position, len(self.chromatine.histones) - 1))

    def cleave_histone(self):
        # Simulate the histone cleavage process by RNase
        if self.chromatine.histones[self.position] == 'H':
            self.chromatine.histones[self.position] = ' '

def visualize_chromatine(histones):
    # Display chromatine state as a bar chart
    plt.bar(range(len(histones)), [1] * len(histones), color=[('gray' if hist == ' ' else 'blue') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")
    plt.show()

def animate_simulation(chromatine, rnases, steps):
    for step in range(steps):
        for rnase in rnases:
            rnase.move()
            rnase.cleave_histone()

        # Visualize the chromatine state
        visualize_chromatine(chromatine.histones)

        # Pause to create an animation effect
        time.sleep(0.5)

        # Clear the output for the next frame
        clear_output(wait=True)

# Parameters for simulation
chromatine_size = 50
rnase_count = 1
simulation_steps = 20  # Reduced steps for demonstration purposes

# Initialize chromatine and RNase
chromatine = Chromatine(chromatine_size)
rnases = [RNase(chromatine) for _ in range(rnase_count)]

# Execute the simulation with animation
animate_simulation(chromatine, rnases, simulation_steps)

# Display the final state of chromatine
visualize_chromatine(chromatine.histones)