from genetic.genetic import TSPGeneticAlgorithm
import numpy as np

N = 7
b = np.random.random_integers(0,200,size=(N,N))
distances = (b + b.T)/2
np.fill_diagonal(distances, 0)

TSPGeneticAlgorithm(
    population_size=100,
    max_epochs=1000
).call(distances)
