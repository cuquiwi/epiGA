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


coordinates = load_problem_file('res/berlin52.tsp')
optimal_path, distance = load_solution_file('res/berlin52.opt.tour', coordinates)

print('Objective distance is:', distance)

TSPGeneticAlgorithm(
    population_size=500,
    mutation_rate=0.5,
    elitism_rate=0.2,
    max_epochs=1000
).call(coordinates, optimal_path)

print('Objective distance is:', objective_distance)
