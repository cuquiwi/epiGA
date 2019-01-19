from genetic.epi_genetic import EpigeneticAlgorithm
from genetic.genetic import TSPGeneticAlgorithm
from util.metrics_util import GeneticMetricPrinter
from data_loader import load_solution_file, load_problem_file
import time

from multiprocessing.pool import Pool
from threading import active_count
from functools import partial
import json

def _run(i, alg, coordinates, optimal_path):
    start = time.time()
    result = alg.call(coordinates, optimal_path)
    elapsed = time.time() - start
    print(f'[{i}] Obtained distance is: {result} in {elapsed}')
    return result, elapsed

for problem in ['circle25', 'circle50', 'berlin52', 'a280']:
    coordinates = load_problem_file(f'res/{problem}.tsp')
    optimal_path, distance = load_solution_file(f'res/{problem}.opt.tour', coordinates)
    
    for alg_name in ['epiAlg', 'gaAlg']:
        
        if alg_name == 'epiAlg':
            alg = EpigeneticAlgorithm(
                     individuals_number=100,
                     cells_number=10,
                     nucleo_prob=0.02,
                     nucleo_rad=2,
                     mechanisms=['position', 'imprinting', 'reprogramming'],
                     epi_probs=[.2, .2, 0.1],
                     max_epochs=500
                )
        else:
            alg = TSPGeneticAlgorithm(
                    population_size=1000,
                    mutation_rate=0.1,
                    max_epochs=500
                )

        with Pool(active_count()) as pool:
            results = pool.map(
                partial(_run, alg=alg, coordinates=coordinates, optimal_path=optimal_path),
                list(range(20))
            )

        with open(f'results_{problem}_{alg_name}.json', 'w') as outfile:
            json.dump(results, outfile)