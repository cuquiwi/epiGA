from genetic.epi_genetic import EpigeneticAlgorithm
from genetic.genetic import TSPGeneticAlgorithm
from util.metrics_util import GeneticMetricPrinter
import numpy as np
from data_loader import load_solution_file, load_problem_file

coordinates = load_problem_file('res/berlin52.tsp')
optimal_path, distance = load_solution_file('res/berlin52.opt.tour', coordinates)

print('Objective distance is:', distance)

alg = EpigeneticAlgorithm(
    individuals_number=100,
    cells_number=10,

    nucleo_prob=0.02,

    nucleo_rad=3,
    mechanisms=['position', 'imprinting', 'reprogramming'],
    epi_probs=[1, 1, 0.1],
    position_prob=.2,
    imprinting_prob=0.2,

    max_epochs=500
)

# alg = TSPGeneticAlgorithm(
#     population_size=100,
#     mutation_rate=0.3,
#     max_epochs=500
# )

alg.subscribe(GeneticMetricPrinter())

result = alg.call(coordinates, optimal_path)

print('Obtained distance is:', result)
