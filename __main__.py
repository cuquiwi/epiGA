from genetic.epi_genetic import EpigeneticAlgorithm
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
optimal_path, distance = load_solution_file(
    'res/berlin52.opt.tour', coordinates)

print('Objective distance is:', distance)

EpigeneticAlgorithm(
    individuals_number = 100,
    cells_number = 10,
    epi_probs = [0.8],
    nucleo_prob = 0.4,
    position_prob=0.7,
    nucleo_rad = 1,
    mechanisms = ['position'],
    max_epochs=500
).call(coordinates, optimal_path)

print('Objective distance is:', distance)
