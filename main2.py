from genetic.epi_genetic import EpigeneticAlgorithm
from genetic.genetic import TSPGeneticAlgorithm
from util.metrics_util import GeneticMetricPrinter
import numpy as np
from data_loader import load_solution_file, load_problem_file

coordinates = load_problem_file('res/a280.tsp')
optimal_path, distance = load_solution_file('res/a280.opt.tour', coordinates)

print('Objective distance is:', distance)

# alg = EpigeneticAlgorithm(
#     individuals_number=50,
#     cells_number=25,

#     nucleo_prob=0.02,

#     nucleo_rad=3,
#     mechanisms=['position', 'imprinting', 'reprogramming'],
#     epi_probs=[1, 1, 0.3],
#     position_prob=.5,

#     max_epochs=500
# )

alg = TSPGeneticAlgorithm(
    population_size=200,
    mutation_rate=0.04,
    max_epochs=1000
)

alg.subscribe(GeneticMetricPrinter())

alg.call(coordinates, optimal_path)

print('Objective distance is:', distance)
