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


distances = load_problem_file('res/a280.tsp')
optimal_path = load_solution_file('res/a280.opt.tour')

objective_distance = 0.0
for i in range(len(optimal_path) - 1):
    c_from = optimal_path[i]
    c_to = optimal_path[i + 1]
    objective_distance += distances[c_to][c_from]

print('Objective distance is:', objective_distance)


TSPGeneticAlgorithm(
    population_size=1000,
    mutation_rate=0.1,
    max_epochs=50
).call(distances)

print('Objective distance is:', objective_distance)
