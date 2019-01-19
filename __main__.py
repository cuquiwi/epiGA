from genetic.epi_genetic import EpigeneticAlgorithm
from genetic.genetic import TSPGeneticAlgorithm
from util.metrics_util import GeneticMetricPrinter
import numpy as np
from data_loader import load_solution_file, load_problem_file
import time

from multiprocessing.pool import Pool
from threading import active_count
from functools import partial

coordinates = load_problem_file('res/circle25.tsp')
optimal_path, distance = load_solution_file('res/circle25.opt.tour', coordinates)

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

def _run(i, alg, coordinates, optimal_path):
    start = time.time()
    result = alg.call(coordinates, optimal_path)
    elapsed = time.time() - start
    print(f'[{i}] Obtained distance is:', result)
    results.append(result)
    return result, elapsed

with Pool(active_count()) as pool:
    results = pool.map(
        partial(_run, alg=alg, coordinates=coordinates, optimal_path=optimal_path),
        list(range(20))
    )

print('Results')
print(results)
