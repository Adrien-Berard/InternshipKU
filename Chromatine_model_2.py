import random
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

class Chromatine:
    def __init__(self, histones_count):
        # Initialize chromatine with histones
        self.histones = np.array(['M' for _ in range(histones_count)])

    def regenerate_histones(self, positions):
        # Replace histones at specified positions with new histones ('A')
        self.histones[positions] = 'A'

    def add_polymerases(self, count, existing_polymerase_positions, adding_position):
        # Add a specified number of new polymerases at non-overlapping positions
        new_positions = np.arange(adding_position, adding_position + count)
        existing_polymerase_positions = np.concatenate([existing_polymerase_positions, new_positions])

    def adding_poly_proba(self, adding_position):
        # Linear function of the local density of histones
        vicinity_size = 5
        start_index = max(0, adding_position - vicinity_size)
        end_index = min(self.histones.shape[0], adding_position + vicinity_size + 1)

        local_density = np.sum(np.isin(self.histones[start_index:end_index], ['M', 'A']))

        slope = 1e-5
        intercept = 1e-2
        probability = slope * local_density + intercept

        return probability

    def change_next_histones(self, position):
        if 0 <= position < self.histones.shape[0] - 1:
            current_histone = self.histones[position]
            next_histone = self.histones[position + 1]

            if current_histone == 'A' and next_histone == 'M':
                self.histones[position + 1] = 'A'

            if position >= 2:
                two_back_histone = self.histones[position - 2]
                if two_back_histone == 'M':
                    self.histones[position + 1] = 'A'

            if position - 1 >= 0:
                one_up_histone = self.histones[position - 1]
                if one_up_histone == 'M':
                    self.histones[position + 1] = 'A'

            if current_histone == 'M' and next_histone == 'A':
                self.histones[position + 1] = 'M'

            if position >= 2:
                two_back_histone = self.histones[position - 2]
                if two_back_histone == 'A':
                    self.histones[position + 1] = 'M'

            if position - 1 >= 0:
                one_up_histone = self.histones[position - 1]
                if one_up_histone == 'A':
                    self.histones[position + 1] = 'M'

class Polymerase:
    def __init__(self, chromatine, position=10, temperature=1.0):
        self.chromatine = chromatine
        self.position = position
        self.temperature = temperature

    def delete(self):
        polymerases.remove(self)

    def move(self, chromatine):
        states = [0, 1]
        next_positions = self.position + states
        probabilities = [1/2, 30]
        total_prob = sum(probabilities)
        normalized_probabilities = [prob / total_prob for prob in probabilities]
        self.position = np.random.choice(next_positions, p=normalized_probabilities)

        if self.position == chromatine.histones.shape[0]:
            self.delete()

    def change_histones(self, chromatine):
        if 0 <= self.position < chromatine.histones.shape[0] and chromatine.histones[self.position] == 'M':
            chromatine.histones[self.position] = 'A'

def visualize_chromatine(histones, polymerase_positions=None):
    plt.bar(range(histones.shape[0]), [1] * histones.shape[0],
            color=[('gray' if hist == ' ' else 'blue' if hist == 'M' else 'green' if hist == 'A' else 'red') for hist in histones])
    plt.title("Chromatine Structure")
    plt.xlabel("Histone Position")
    plt.ylabel("Histone")

    if polymerase_positions:
        for position in polymerase_positions:
            plt.arrow(position, -0.2, 0, -0.1, color='red', head_width=0.5, head_length=0.1)

def update(frame):
    deleted_positions = np.array('')
    polymerase_positions = np.array('')

    for polymerase in polymerases:
        polymerase.move(chromatine)
        polymerase.change_histones(chromatine)
        deleted_positions = np.append(deleted_positions, polymerase.position)
        polymerase_positions = np.append(polymerase_positions, polymerase.position)

    if random.random() < chromatine.adding_poly_proba(adding_position):
        chromatine.add_polymerases(1, existing_polymerase_positions, adding_position)
        new_polymerase_positions = existing_polymerase_positions[-1:]
        new_polymerases = np.array([Polymerase(chromatine, position=pos, temperature=1.0) for pos in new_polymerase_positions])
        polymerases = np.concatenate([polymerases, new_polymerases])

    chromatine.regenerate_histones(deleted_positions)

    for position in existing_polymerase_positions:
        chromatine.change_next_histones(position)

    polymerase_count_over_time = np.append(polymerase_count_over_time, polymerases.shape[0])
    active_histone_count = np.sum(np.isin(chromatine.histones, ['M', 'A']))
    acetylated_histone_count = np.sum(chromatine.histones == 'A')
    methylated_histone_count = np.sum(chromatine.histones == 'M')
    active_histone_count_over_time = np.append(active_histone_count_over_time, active_histone_count)
    acetylated_histone_count_over_time = np.append(acetylated_histone_count_over_time, acetylated_histone_count)
    methylated_histone_count_over_time = np.append(methylated_histone_count_over_time, methylated_histone_count)

    polymerase_count_over_time = np.pad(polymerase_count_over_time, (0, frame + 2 - polymerase_count_over_time.shape[0]), 'constant')
    acetylated_histone_count_over_time = np.pad(acetylated_histone_count_over_time, (0, frame + 1 - acetylated_histone_count_over_time.shape[0]), 'constant')
    methylated_histone_count_over_time = np.pad(methylated_histone_count_over_time, (0, frame + 1 - methylated_histone_count_over_time.shape[0]), 'constant')
    active_histone_count_over_time = np.pad(active_histone_count_over_time, (0, frame + 1 - active_histone_count_over_time.shape[0]), 'constant')

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
    axs[0, 1].set_title('Number of Active Histones Over Time')
    axs[0, 1].set_xlabel('Time Steps')
    axs[0, 1].set_ylabel('Number of Histones')
    axs[0, 1].legend()

    visualize_chromatine(chromatine.histones, polymerase_positions=polymerase_positions)

chromatine_size = 50
polymerase_count = 0
simulation_steps = 100
adding_position = 25

chromatine = Chromatine(chromatine_size)
polymerases = np.array([Polymerase(chromatine, temperature=1.0) for _ in range(polymerase_count)])
existing_polymerase_positions = np.array([polymerase.position for polymerase in polymerases])

polymerase_count_over_time = np.array('')
active_histone_count_over_time = np.array('')
acetylated_histone_count_over_time = np.array('')
methylated_histone_count_over_time = np.array('')

fig, axs = plt.subplots(2, 2, figsize=(12, 8))
ani = FuncAnimation(fig, update, frames=simulation_steps, repeat=False)
plt.show()
