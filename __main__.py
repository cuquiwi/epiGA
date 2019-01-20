from genetic.epi_genetic import EpigeneticAlgorithm
from genetic.genetic import TSPGeneticAlgorithm
from util.metrics_util import GeneticMetricPrinter
from data_loader import load_solution_file, load_problem_file

if __name__ == '__main__':

    coordinates = load_problem_file(f'res/circle25.tsp')
    optimal_path, _ = load_solution_file(f'res/circle25.opt.tour', coordinates)

    alg = EpigeneticAlgorithm(
                individuals_number=100,
                cells_number=10,
                nucleo_prob=0.02,
                nucleo_rad=int(len(optimal_path)*.05),
                mechanisms=['position', 'imprinting', 'reprogramming'],
                epi_probs=[.2, .2, 0.1],
                max_epochs=500
        )
    # alg = TSPGeneticAlgorithm(
    #         population_size=1000,
    #         mutation_rate=0.1,
    #         max_epochs=500
    #     )

    alg.subscribe(GeneticMetricPrinter())
    alg.call(coordinates, optimal_path)