from genetic.epi_genetic import EpigeneticAlgorithm
from genetic.genetic import TSPGeneticAlgorithm
from util.metrics_util import GeneticMetricPrinter
import numpy as np
from data_loader import load_solution_file, load_problem_file
import time

coordinates = load_problem_file('res/circle50.tsp')
optimal_path, distance = load_solution_file('res/circle50.opt.tour', coordinates)

print('Objective distance is:', distance)

# alg = EpigeneticAlgorithm(
#     individuals_number=100,
#     cells_number=10,

#     nucleo_prob=0.02,

#     nucleo_rad=3,
#     mechanisms=['position', 'imprinting', 'reprogramming'],
#     epi_probs=[.2, .2, 0.1],

#     max_epochs=500
# )

alg = TSPGeneticAlgorithm(
    population_size=1000,
    mutation_rate=0.3,
    max_epochs=200
)

alg.subscribe(GeneticMetricPrinter())

perfs = []
results = []
for _ in range(20):
    start = time.time()
    result = alg.call(coordinates, optimal_path)
    perfs.append(time.time() - start)
    print('-Obtained distance is:', result)
    results.append(result)

print('Results')
print(results)
print('Perfs')
print(perfs)