from genetic.genetic import TSPGeneticAlgorithm
import numpy as np
from data_loader import load_solution_file, load_problem_file

distances = [
    # 0  1  2  3
    [00, 10,  1, 10],  # 0
    [10, 00, 10,  1],  # 1
    [1, 10, 00,  1],  # 2
    [10,  1,  1, 00],  # 3
]
optimal_path = [0, 2, 3, 1]


distances, coordinates = load_problem_file('res/berlin52.tsp')
optimal_path = load_solution_file('res/berlin52.opt.tour')

objective_distance = 0.0
for i in range(len(optimal_path) - 1):
    c_from = optimal_path[i]
    c_to = optimal_path[i + 1]
    objective_distance += distances[c_to][c_from]

print('Objective distance is:', objective_distance)


TSPGeneticAlgorithm(
<<<<<<< HEAD
    population_size=100,
    mutation_rate=0.1,
    max_epochs=2000
=======
    population_size=500,
    mutation_rate=0.5,
    elitism_rate=0.2,
    max_epochs=1000
>>>>>>> 5bdc66eeb607ad5c01d9d997cce5c596897683a3
).call(distances, coordinates, optimal_path)

print('Objective distance is:', objective_distance)
