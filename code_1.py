import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Chromatine:
    def __init__(self, histones_count):
        self.histones = ['H' for _ in range(histones_count)]

class RNase:
    def __init__(self, chromatine, position=0):
        self.chromatine = chromatine
        self.position = position

    def move(self):
        movement = random.choice([-1, 1])
        self.position += movement
        self.position = max(0, min(self.position, len(self.chromatine.histones) - 1))

    def cleave_histone(self):
        if self.chromatine.histones[self.position] == 'H':
            self.chromatine.histones[self.position] = ' '

def visualize_chromatine(histones):
    plt.bar(range(len(histones)), [1] * len(histones),
            color=[('gray' if hist == ' ' else 'blue') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")

# Function to update the plot in each animation frame
def update(frame):
    for rnase in rnases:
        rnase.move()
        rnase.cleave_histone()

    plt.clf()  # Clear the previous frame
    visualize_chromatine(chromatine.histones)

# Parameters for simulation
chromatine_size = 50
rnase_count = 10
simulation_steps = 20

# Initialize chromatine and RNase
chromatine = Chromatine(chromatine_size)
rnases = [RNase(chromatine) for _ in range(rnase_count)]

# Create an animated plot
fig, ax = plt.subplots()
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()
