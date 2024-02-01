import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Chromatine:
    def __init__(self, histones_count):
        self.histones = ['H' for _ in range(histones_count)]
        self.new_histones = set()

    def regenerate_histones(self, positions):
        for pos in positions:
            if pos not in self.new_histones:
                self.new_histones.add(pos)
                self.histones[pos] = 'N'

class RNase:
    def __init__(self, chromatine, position=0):
        self.chromatine = chromatine
        self.position = position

    def move(self):
        movement = random.choice([-1, 1])
        self.position += movement
        self.position = max(0, min(self.position, len(self.chromatine.histones) - 1))

    def cleave_histone(self):
        if self.chromatine.histones[self.position] in {'H', 'N'}:
            self.chromatine.histones[self.position] = ' '

def visualize_chromatine(histones, rnase_position=None):
    plt.bar(range(len(histones)), [1] * len(histones),
            color=[('gray' if hist == ' ' else 'blue' if hist == 'H' else 'green' if hist == 'N' else 'red') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")

    if rnase_position is not None:
        plt.arrow(rnase_position, 1.2, 0, 0.1, color='red', head_width=0.5, head_length=0.1)

# Function to update the plot in each animation frame
def update(frame):
    deleted_positions = []  # Keep track of deleted positions for regeneration
    for rnase in rnases:
        rnase.move()
        rnase.cleave_histone()
        deleted_positions.append(rnase.position)

    plt.clf()  # Clear the previous frame

    # Visualize chromatine structure with arrow indicating RNase position
    visualize_chromatine(chromatine.histones, rnase_position=rnases[0].position)

    # Regenerate histones at deleted positions
    chromatine.regenerate_histones(deleted_positions)

# Parameters for simulation
chromatine_size = 50
rnase_count = 1
simulation_steps = 20

# Initialize chromatine and RNase
chromatine = Chromatine(chromatine_size)
rnases = [RNase(chromatine) for _ in range(rnase_count)]

# Create an animated plot
fig, ax = plt.subplots()
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()

# Print a message when the simulation is done
print("Simulation is done!")

# Visualize the final state of the chromatine
plt.figure()
visualize_chromatine(chromatine.histones)
plt.title("Final Chromatine Structure")
plt.show()
