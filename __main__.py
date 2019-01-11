from genetic.epi_genetic import EpigeneticAlgorithm
from util.metrics_util import GeneticMetricPrinter
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

alg = EpigeneticAlgorithm(
    individuals_number=50,
    cells_number=25,

    nucleo_prob=0.02,

    nucleo_rad=3,
    mechanisms=['position', 'imprinting', 'reprogramming'],
    epi_probs=[1, 1, 0.01],
    position_prob=.5,

    max_epochs=500
)

alg.subscribe(GeneticMetricPrinter())

alg.call(coordinates, optimal_path)

print('Objective distance is:', distance)
