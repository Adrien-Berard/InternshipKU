import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def potential_gradient(positions, bonds, l0, b):
    repulsion_gradient = -4 * np.sum(np.sign(positions[:, np.newaxis] - positions) *
                                     np.exp(-4 * np.abs(positions[:, np.newaxis] - positions) / l0) / l0, axis=1, dtype=np.float128)

    bonding_gradient = 4 * b * np.sum(np.sign(positions[:, np.newaxis] - positions[bonds[:, 0]]) *
                                          np.exp(-4 * b * np.abs(positions[:, np.newaxis] - positions[bonds[:, 0]]) / l0) / l0, axis=1, dtype=np.float128)

    return repulsion_gradient - bonding_gradient

class Particle:
    def __init__(self, mass=1.0):
        self.mass = mass
        self.position = np.zeros(3)
        self.velocity = np.zeros(3)
    
    def calculate_force(self, positions, bonds, l0, b):
        return potential_gradient(positions, bonds, l0, b)

class PolymerSimulation:
    def __init__(self, num_beads, temperature, damping_coefficient):
        self.num_beads = num_beads
        self.temperature = temperature
        self.damping_coefficient = damping_coefficient

        self.particles = [Particle() for _ in range(num_beads)]
        self.trajectory = []  # List to store particle positions at each timestep

    def initialize_positions(self):
        for particle in self.particles:
            particle.position = np.random.rand(3)

    def euler_maruyama_3d(self, dt, bonds, l0, b):
        random_forces = np.random.normal(0, np.sqrt(2 * self.temperature * self.damping_coefficient / dt), (self.num_beads, 3))
        
        positions_array = np.vstack([particle.position for particle in self.particles])
        for i in range(1, self.num_beads):
            damping_forces = -self.damping_coefficient * self.particles[i - 1].velocity
            
            potential_gradient_forces = self.particles[i].calculate_force(
                positions_array,
                bonds,
                l0,
                b
            )
            
            equilibrium_distance = 1.0
            distance_correction_forces = (positions_array[i - 1, :] - positions_array[i, :]) * (np.linalg.norm(positions_array[i - 1, :] - positions_array[i, :]) - equilibrium_distance)
            
            stochastic_forces = random_forces[i, :]

            self.particles[i].position = self.particles[i - 1].position + self.particles[i - 1].velocity * dt
            self.particles[i].velocity = self.particles[i - 1].velocity + (damping_forces + potential_gradient_forces + distance_correction_forces + stochastic_forces) * dt
            
            self.trajectory.append([particle.position.copy() for particle in self.particles])

    def run_simulation(self, timesteps, dt, bonds, l0, b):
        self.initialize_positions()

        for _ in range(timesteps):
            self.euler_maruyama_3d(dt, bonds, l0, b)
            print(_)

    def plot_trajectory(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for bead in range(self.num_beads):
            ax.plot(
                [particle.position[0] for particle in self.particles],
                [particle.position[1] for particle in self.particles],
                [particle.position[2] for particle in self.particles],
                label=f'Bead {bead + 1}'
            )

        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.set_zlabel('Z-axis')
        ax.legend()
        plt.show()

    def save_trajectory_to_csv(self, filename='molecularPolymerTrajectory.xyz'):
        with open(filename, 'w') as file:
            for timestep, positions in enumerate(self.trajectory):
                file.write(f"{self.num_beads}\n")
                file.write(f"Timestep {timestep}\n")
                for position in positions:
                    file.write(f"C {position[0]} {position[1]} {position[2]}\n")

# Parameters
num_beads = 10
temperature = 300.0  # Kelvin
damping_coefficient = 2
bonds = np.array([[i, i + 1] for i in range(num_beads - 1)])
l0 = 1.0
b = 1.0

# Create and run simulation
polymer_sim = PolymerSimulation(num_beads, temperature, damping_coefficient)
polymer_sim.run_simulation(timesteps=1000, dt=0.01, bonds=bonds, l0=l0, b=b)

# Save the trajectory to a CSV file
polymer_sim.save_trajectory_to_csv()

# Plot 3D trajectory
# polymer_sim.plot_trajectory()