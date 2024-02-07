import random
import numpy as np
import pandas as pd

class Chromatine:
    def __init__(self, histones_count):
        # Initialize chromatine with histones
        self.histones = np.full(histones_count, 'M', dtype='U1')

        # Randomly set approximately 30% of histones to 'U' (unmodified)
        num_unmodified = int(0.3 * histones_count)
        unmodified_positions = np.random.choice(histones_count, size=num_unmodified, replace=False)
        self.histones[unmodified_positions] = 'U'

    def regenerate_histones(self, positions_to_delete):
        unmodified_positions = np.where(self.histones == 'U')[0]

        # Ensure that the size of the sample is not greater than the size of the population
        if len(unmodified_positions) >= len(positions_to_delete):
            selected_positions = np.random.choice(unmodified_positions, size=len(positions_to_delete), replace=False)
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
        slope = 1e-5  # Adjust the slope to control the influence of local density
        intercept = 1e-1 # Adjust the intercept to control the baseline probability

        # Calculate the probability of adding a new polymerase
        probability = slope * local_density + intercept

        return probability

    def change_next_histones(self, position, p_recruitment, p_change, nth_neighbor=1, vicinity_size=5):
        # Simulate the influence of neighbors on the next histones
        if 1 <= position < len(self.histones) - 1:
            current_histone = self.histones[position]

            # Calculate the influence of vicinity on the recruitment probability
            adjusted_p_recruitment = p_recruitment / (1e-1 * nth_neighbor)  # Add 1 to avoid division by zero
            if adjusted_p_recruitment > 1:
                adjusted_p_recruitment = 1

            # Probabilistically recruit an enzyme
            if np.random.random() < adjusted_p_recruitment:
                # Apply changes with probability p_change
                if np.random.random() < p_change:
                    nth_position = position + nth_neighbor
                    if nth_position < len(self.histones):
                        nth_histone = self.histones[nth_position]

                        if current_histone == 'U' and nth_histone == 'A':
                            self.histones[nth_position] = 'U'
                            # Implement changes based on the nth neighbor
                            transition_key = f"{current_histone}_and_{nth_histone}_to_{current_histone}_and_{current_histone}"
                            transitions_dict[transition_key] = transitions_dict.get(transition_key, 0) + 1
                        elif current_histone == 'U' and nth_histone == 'M':
                            self.histones[nth_position] = 'U'
                            # Implement changes based on the nth neighbor
                            transition_key = f"{current_histone}_and_{nth_histone}_to_{current_histone}_and_{current_histone}"
                            transitions_dict[transition_key] = transitions_dict.get(transition_key, 0) + 1
                        elif current_histone == 'A' and nth_histone == 'U':
                            self.histones[nth_position] = 'A'
                            # Implement changes based on the nth neighbor
                            transition_key = f"{current_histone}_and_{nth_histone}_to_{current_histone}_and_{current_histone}"
                            transitions_dict[transition_key] = transitions_dict.get(transition_key, 0) + 1
                        elif current_histone == 'A' and nth_histone == 'M':
                            self.histones[nth_position] = 'A'
                            # Implement changes based on the nth neighbor
                            transition_key = f"{current_histone}_and_{nth_histone}_to_{current_histone}_and_{current_histone}"
                            transitions_dict[transition_key] = transitions_dict.get(transition_key, 0) + 1
                        elif current_histone == 'M' and nth_histone == 'A':
                            self.histones[nth_position] = 'M'
                            # Implement changes based on the nth neighbor
                            transition_key = f"{current_histone}_and_{nth_histone}_to_{current_histone}_and_{current_histone}"
                            transitions_dict[transition_key] = transitions_dict.get(transition_key, 0) + 1
                        elif current_histone == 'M' and nth_histone == 'U':
                            self.histones[nth_position] = 'M'
                            # Implement changes based on the nth neighbor
                            transition_key = f"{current_histone}_and_{nth_histone}_to_{current_histone}_and_{current_histone}"
                            transitions_dict[transition_key] = transitions_dict.get(transition_key, 0) + 1

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

        next_position = self.position + 1  # Check only one position ahead

        # Check if there is another polymerase in the next position
        if next_position in existing_polymerase_positions:
            # Do not move if there is another polymerase 1 place after
            return

        probabilities = [1/2, 1/2]

        # Normalize probabilities to sum to 1
        total_prob = np.sum(probabilities)
        normalized_probabilities = probabilities / total_prob

        # Choose the next position based on normalized probabilities
        self.position = np.random.choice([self.position, next_position], p=normalized_probabilities)

        # Bounding conditions
        if self.position == end_of_replication_position:
            self.delete()

    def change_histones(self, chromatine):
        # Simulate the histone change process by polymerase
        if 0 <= self.position < len(chromatine.histones) and chromatine.histones[self.position] == 'M':
            chromatine.histones[self.position] = 'A'
            transition_key = f"Polymerase_M_to_A"
            transitions_dict[transition_key] = transitions_dict.get(transition_key, 0) + 1

# Set seed for reproducibility
np.random.seed(42)

# Parameters for simulation
chromatine_size = 50
polymerase_count = 0
simulation_steps = 200
adding_position = 10
end_of_replication_position = chromatine_size - 10

# Initialize chromatine and polymerases with a specified temperature
chromatine = Chromatine(chromatine_size)
polymerases = [Polymerase(chromatine, temperature=1.0) for _ in range(polymerase_count)]

# Track existing polymerase positions using a list to avoid duplicates
existing_polymerase_positions = [polymerase.position for polymerase in polymerases]

# Dictionary to store transitions and their counts
transitions_dict = {}

# Simulate the process
for frame in range(simulation_steps):
    deleted_positions = []  # Keep track of deleted positions for regeneration

    for polymerase in polymerases:
        polymerase.move(chromatine)
        polymerase.change_histones(chromatine)
        deleted_positions.append(polymerase.position)

    # Change the next histones based on the influence of first neighbors
    for position in range(1, chromatine_size):
        # Use p_recruitment and p_change probabilities with decreasing probability with vicinity
        chromatine.change_next_histones(position, p_recruitment=0.2, p_change=0.3, nth_neighbor=np.random.randint(1, chromatine_size), vicinity_size=5)

    # Regenerate histones at unmodified positions
    # if np.random.random() < 0.4:
    #    chromatine.regenerate_histones(deleted_positions)

    # Randomly add new polymerase at the beginning of the chromatine with a certain probability
    if np.random.random() < chromatine.adding_poly_proba(adding_position):  # Adjust the probability as needed
        # Add new polymerases with non-overlapping random positions
        chromatine.add_polymerases(1, existing_polymerase_positions, adding_position)
        new_polymerase_positions = existing_polymerase_positions[-1:]
        new_polymerases = [Polymerase(chromatine, position=pos, temperature=1.0) for pos in new_polymerase_positions]
        polymerases.extend(new_polymerases)

# Create a DataFrame from the transitions dictionary
transitions_df = pd.DataFrame(list(transitions_dict.items()), columns=['Transition', 'Count'])

# Save the DataFrame to a CSV file
transitions_df.to_csv('transitions_data.csv', index=False)
