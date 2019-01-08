from genetic.genetic import TSPGeneticAlgorithm
import numpy as np
from data_loader import load_solution_file, load_problem_file

distances = load_problem_file('res/a280.tsp')
optimal_path = load_solution_file('res/a280.opt.tour')

objective_distance = 0.0
for i in range(len(optimal_path) - 1):
    c_from = optimal_path[i]
    c_to = optimal_path[i + 1]
    objective_distance += distances[c_to][c_from]

print(objective_distance)


TSPGeneticAlgorithm(
    population_size=100,
    max_epochs=1000
).call(distances)


print(objective_distance)