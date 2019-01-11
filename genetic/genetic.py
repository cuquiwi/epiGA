from random import sample, randint, random, shuffle
import numpy as np
import matplotlib.pyplot as plt


class TSPIndividual(object):

    def __init__(self, path):
        self.path = path
        self.fitness = None
        self.distance = None

    def __repr__(self):
        return f'<TSPIndividual({self.fitness}, {self.distance}) {self.path}>'


class ITSPGeneticAlgorithm(object):

    def __init__(self, population_size=100, mutation_rate=.01, elitism_rate=0.2,
                 max_epochs=10000):
        self.population_size = population_size
        self.mutation_rate = mutation_rate
        self.max_epochs = max_epochs
        self.elitism_rate = elitism_rate
        self.subscriptions = []

    def call(self, coordinates, optimum_path):
        distance_matrix = self.calculate_distances(coordinates)
        population = self.do_initialize_population(len(distance_matrix))

        self.on_launch(coordinates, optimum_path)

        i = 0
        for _ in range(self.max_epochs):
            self.do_calculate_fitness(population, distance_matrix)

            self.on_epoch(population, i)

            new_population = []
            for _ in range(self.population_size):
                parents = self.do_pick_parents(population)
                child = self.do_crossover(*parents)
                child = self.do_mutation(child)
                new_population.append(child)

            population = self.do_next_generation(population, new_population)

            i = i+1

        return population

    def calculate_distances(self, coordinates):
        distance_matrix = np.zeros((len(coordinates), len(coordinates)))
        for i in range(len(coordinates)):
            for j in range(len(coordinates)):
                distance_matrix[i][j] = np.linalg.norm(
                    coordinates[i] - coordinates[j]
                )
        return distance_matrix

    def do_initialize_population(self, n_cities):
        raise NotImplementedError()

    def do_calculate_fitness(self, population, distances):
        raise NotImplementedError()

    def do_pick_parents(self, population):
        raise NotImplementedError()

    def do_crossover(self, *parents):
        raise NotImplementedError()

    def do_mutation(self, individual):
        raise NotImplementedError()

    def do_next_generation(self, olders, newers):
        raise NotImplementedError()

    def on_epoch(self, population, i):
        raise NotImplementedError()

    def on_launch(self, coordinates, optimum_path):
        raise NotImplementedError()


class TSPGeneticAlgorithm(ITSPGeneticAlgorithm):

    def do_initialize_population(self, n_cities):
        return [
            TSPIndividual(sample(range(n_cities), n_cities))
            for _ in range(self.population_size)
        ]

    def do_calculate_fitness(self, population, distances):
        """Calculates the fitness of the individual based on the distance between
        the cities

        Arguments:
            population {List of TSPIndividual} -- Population
            distances {NxN float matrix} -- Distances between cities
        """

        for individual in population:
            total_distance = 0.0
            for i in range(len(individual.path)-1):
                origin = individual.path[i]
                destination = individual.path[i+1]
                total_distance += distances[origin][destination]
            individual.fitness = 1 / total_distance
            individual.distance = total_distance
        self.normalize_fitness(population)

    def normalize_fitness(self, population):
        total_fitness = sum([ind.fitness for ind in population])
        for individual in population:
            individual.fitness = individual.fitness / (total_fitness / 100)

    def do_pick_parents(self, population):
        return [self.pick_parent(population), self.pick_parent(population)]

    def pick_parent(self, population):
        """Performs the selection over the population by choosing at random based
        on the fitness of the individuals. That is an individual with high fitness
        has more probabilities of being picked.

        Arguments:
            population {List of TSPIndividual} -- Population with fitness calculated

        Returns:
            TSPIndividual -- Selected
        """

        counter = randint(0, 100)
        for individual in population:
            counter -= individual.fitness
            if counter <= 0:
                return individual
        return individual

    def do_crossover(self, *parents):
        """Performs a crossover among two parents, the crossover is designed to 
        pick a random range from the first parent, use it as it is and fill the 
        the remaining space with the cities left in the order of the second path

        [4,2,3,1]
                (pick from 1 to 2)[1,2,3,4]
        [3,2,1,4]

        Returns:
            List of TSPIndividual -- Parents to crossover
        """

        path_1 = parents[0].path
        path_2 = parents[1].path

        rand1 = randint(0, len(path_1))
        rand2 = randint(0, len(path_1))

        ini = min(rand1, rand2)
        end = max(rand1, rand2)

        new_path = [-1 for _ in range(len(path_2))]
        new_path[ini:end] = path_1[ini:end]

        for city in [city for city in path_2 if city not in new_path]:
            for i in range(len(new_path)):
                if new_path[i] == -1:
                    new_path[i] = city
                    break

        return TSPIndividual(new_path)

    def do_mutation(self, individual):
        """Performs a swap mutation in the path at random

        [ 1, 2, 3] -> [1, 3, 2]

        Arguments:
            individual {TSPIndividual} -- Individual to mutate

        Returns:
            TSPIndividual -- Mutated Individual
        """

        while random() <= self.mutation_rate:
            origin = randint(0, len(individual.path)-1)
            end = randint(0, len(individual.path)-1)
            individual.path[origin], individual.path[end] = individual.path[end], individual.path[origin]
        return individual

    def do_next_generation(self, olders, newers):
        """
        Performs the strategy to select the future generation.

        Inputs:
            - olders: old generation.
            - newers: new generation.
        Output:
            Final population of the current epoch.
        """
        shuffle(newers)
        number_of_elite_ind = int(len(olders) * self.elitism_rate)
        best_old = sorted(olders, key=lambda x: x.fitness,
                          reverse=True)[:number_of_elite_ind]
        newers[:number_of_elite_ind] = best_old
        return newers

    def on_launch(self, coordinates, optimum_path):
        [
            sub.on_launch(coordinates, optimum_path, epigenetic=False)
            for sub in self.subscriptions
        ]

    def on_epoch(self, population, i):
        sorted_pop = sorted(population, key=lambda i: i.fitness, reverse=True)
        [
            sub.on_epoch(
                [i.path for i in sorted_pop],
                [i.distance for i in sorted_pop],
                i,
                epigenetic=False
            )
            for sub in self.subscriptions
        ]

    def subscribe(self, subscriber):
        self.subscriptions.append(subscriber)
