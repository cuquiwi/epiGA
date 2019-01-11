from random import shuffle, random, uniform
from .cell import Cell
import numpy as np
import matplotlib.pyplot as plt


class EpigeneticAlgorithm(object):

    def __init__(self, individuals_number, cells_number, epi_probs,
                 nucleo_prob, nucleo_rad, mechanisms,
                 position_prob=0.4, imprinting_prob=0.2,
                 max_epochs=500):
        """
        EpiGA based on the work by D.H. Stolfi and E. Alba, (2017).
        Inputs:
            individuals_number: Number of individuals.
            cells_number: Number of cells per each individual.
            epi_probs: Epigenetic probabilities list.
            nucleo_prob: Nucleosome probability.
            nucleo_rad: Nucleosome radius.
            mechanisms: List of string containing the epigenetic algorithms
                        in the same order as the epiProb list.
            environment: Environment rules
            max_epochs: Maximun number of epochs
            position_prob: Probability for the position epigenetic mechanism
        """
        self.individuals_number = individuals_number
        self.cells_number = cells_number
        self.epi_probs = epi_probs
        self.nucleo_prob = nucleo_prob
        self.nucleo_rad = nucleo_rad
        self.mechanisms = mechanisms
        self.max_epochs = max_epochs
        self.position_prob = position_prob
        self.imprinting_prob = imprinting_prob
        self.xdata_fitness = []
        self.ydata_fitness = []
        self.mean_fitness = []
        self.min_fitness = []
        self.ydata_iter = []

    def call(self, coordinates, optimum_path):
        """Performs the actual algorithm taking into account the configuration 
        provided in the initialization and the distances matrix provided here

        Arguments:
            distances_matrix {Matrix NxN} -- Matrix containing the distances 
                                        between cities
        """

        self.distances_matrix = self.calculate_distances(coordinates)

        plt.ion()
        self.on_launch()

        population = self.init_population()
        i = 0
        fitnesses = []
        while not self.termination(i, fitnesses):

            self.print_metric(population, coordinates, optimum_path, i)

            newpop = self.selection(population)
            newpop = self.nucleosome_generation(newpop)
            newpop = self.nucleosome_reproduction(newpop)
            newpop = self.epigen_mechanism(newpop)
            

            population = self.replacement(population, newpop)

            fitnesses.append(
                np.average([self.evaluate_individual(i) for i in newpop])
            )

    def calculate_distances(self, coordinates):
        """Calculates the initial matrix of distances
        
        Arguments:
            coordinates {Array of Coordinates} -- Coordinates of the city
        
        Returns:
            Matrix NxN -- Matrix of distances
        """

        distance_matrix = np.zeros((len(coordinates), len(coordinates)))
        for i in range(len(coordinates)):
            for j in range(len(coordinates)):
                distance_matrix[i][j] = np.linalg.norm(
                    coordinates[i] - coordinates[j]
                )
        return distance_matrix

    def termination(self, i, fitnesses):
        """
        Termination condition for the EpiGA. It will stop whenever it gets to a 
        maximun epochs or the last fitnesses are equal.

        Inputs:
            - i: Current iteration of the algorithm.
            - fitnesses: fitness in all the iterations
        Output:
            True if the termination condition is accomplished, False
            otherwise.
        """
        if i >= self.max_epochs:
            return True
        if i>= 15 and np.average(fitnesses[:15]) == fitnesses[-1]:
            return True
        return False

    def evaluate_cell(self, cell):
        """
        Function that evaluates the fitness of our problem for a given cell
        and sets the cell fitness.
        Inputs:
            - cell: a cell which contains the solution.
        Outputs:
            The fitness value of a given cell.
        """
        solution = cell.solution
        fitness = 0
        for i in range(1, len(solution)):
            fitness += self.distances_matrix[cell.solution[i-1]
                                             ][cell.solution[i]]
        cell.setfitness(fitness)

        return fitness

    def evaluate_individual(self, individual):
        """
        Get the total fitness of the individual, the sum of the fitness of its cells.
        Input:
            - individual: Individual, represented by a list of cells
        Return:
            - fitness of the individual
        """
        cell_fitness = list(map(lambda cell: cell.fitness, individual))

        return np.min(cell_fitness)

    def init_population(self):
        """
        The initial population of the EpiGA.
        Return:
            A list of cell elements representing the initial population.
        """
        population = []
        for _ in range(self.individuals_number):
            individual = []
            for j in range(self.cells_number):
                solution = [k for k in range(len(self.distances_matrix))]
                shuffle(solution)
                cell = Cell(solution)
                self.evaluate_cell(cell)
                individual.append(cell)

            population.append(individual)

        return population

    def selection(self, population):
        """
        Performs the population selection in the EpiGA.
        Uses a binary tournament selection with no repetitions.
        Returns half of the population
        Inputs:
            - population: the total population.
        Return:
            The selected subset of the total population.
        """
        winners = []
        shuffle(population)
        i = 0
        while i+1 < len(population):
            if self.evaluate_individual(population[i]) > self.evaluate_individual(population[i+1]):
                winners.append(population[i])
            else:
                winners.append(population[i+1])
            i += 2
        return winners

    def k_tournament_selection(self, pop, k=2):
        """
        Performs a tournament selection, by default a binary selection
        Inputs:
            - pop: The total population
            - k: the number of individuals to compete, default is 2
        Return:
            The winner of the k-tournament
        """
        best = None
        bestfit = 0
        for i in range(1, k):
            ind = pop[random(1, len(pop))]
            if (best == None) or self.evaluate_individual(ind) > bestfit:
                best = ind
                bestfit = self.evaluate_individual(best)
        return best

    def nucleosome_generation(self, population):
        """
        Generates a new nucleosome vector as a mask for each cell in the 
        individuals of the population.
        Inputs:
            - population: the total population.
        Output:
            The population with the new nucleosomes generated.
        """
        for i in range(len(population)):
            individual = population[i]
            for j in range(len(individual)):
                cells = individual[j]
                n = np.zeros(len(cells.nucleosome[:]), dtype=bool)
                for k in range(len(n)):
                    if random() < self.nucleo_prob:
                        n = self.collapse(n, k)
                cells.nucleosome = n
                individual[j] = cells
            population[i] = individual
        return population

    def collapse(self, nucleosome, k):
        """
        Function that creates a new nucleosome.
        Inputs:
            - nucleosome: The binary list of the cell that represents the nucleosome.
            - k: the center position of the new nucleosome.
        Output:
            The new modified nucleosome.
        """
        for i in range(-self.nucleo_rad, self.nucleo_rad+1, 1):
            if (k+i in range(len(nucleosome))):
                nucleosome[k+i] = True
        return nucleosome

    def selectBestCell(self, individual):
        """Selects the best cell inside the individual

        Arguments:
            individual {List of Cells} -- Individual to look in

        Returns:
            Cell -- The best cell of the individual
        """

        fitness = list(map(lambda cell: cell.fitness, individual))
        return individual[np.argmin(fitness)]

    def crossover(self, baseSolution, secondSolution, mask):
        """
        This function performs PMX (Partial-Mapping Crossover)
        Inputs:
            - basesolution: solution that is going to be used as base
            - secondSolution: solution that is going to be used to replace
            the areas that are not bent in the chromosome of the baseSolution
            - mask: Defines the areas in which the chromosome is bent (i.e. the 
            places that are going to be passed on to the next generation)
        """
        # Define mask for the substitution of values so that they are not repeated
        mapping = {}
        for i in range(len(mask)):
            # The mask uses values from the bent area of the base solution
            if (mask[i]):
                mapping[baseSolution[i]] = secondSolution[i]

        # Start the new solution as the baseSolution
        newsolution = []
        for elem in baseSolution:
            newsolution.append(elem)

        for j in range(len(newsolution)):
            # If the chromosome is not bent, we must replace with second solution value
            if (not mask[j]):
                city = secondSolution[j]
                # However, if it is going to be a repeated value, we use the generated map
                while (city in mapping):
                    city = mapping[city]
                newsolution[j] = city

        # https://www.researchgate.net/publication/226665831_Genetic_Algorithms_for_the_Travelling_Salesman_Problem_A_Review_of_Representations_and_Operators

        return newsolution

    def removeWorstCell(self, individual, newCell):
        """
        Function that generates a new individual changing its worst cell
        Inputs: 
            - individual: List of Cells that form the base individual from 
                which the new individual is formed
            - newCell: The cell that is going to replace the worst cell of 
                the received individual
        Output:
            New individual based in the received one with the newCell
        """
        newInd = []
        for cell in individual:
            newInd.append(cell)
        fitness = list(map(lambda cell: cell.fitness, individual))
        newInd.remove(individual[np.argmax(fitness)])
        newInd.append(newCell)
        return newInd

    def nucleosome_reproduction(self, population):
        """
        Function that generates the new population based on nucleosome reproduction
        Inputs: 
            - population: list of individuals of the previous generation
        Output: 
            The list of children of the provious population
        """
        newPop = []
        for _ in range(2 * self.individuals_number):
            i1 = self.roulette_selection(population)
            i2 = self.roulette_selection(population)
            bestCell1 = self.selectBestCell(i1)
            bestCell2 = self.selectBestCell(i2)
            newNucleosome = np.logical_or(
                bestCell1.nucleosome, bestCell2.nucleosome)
            fatherBasedSolution = self.crossover(
                bestCell1.solution, bestCell2.solution, newNucleosome
            )
            motherBasedSolution = self.crossover(
                bestCell2.solution, bestCell1.solution, newNucleosome
            )
            newCellI1 = Cell(
                fatherBasedSolution, bestCell1.solution, bestCell2.solution, newNucleosome
            )
            newCellI2 = Cell(
                motherBasedSolution, bestCell2.solution, bestCell1.solution, newNucleosome
            )
            self.evaluate_cell(newCellI1)
            self.evaluate_cell(newCellI2)
            i1_child = self.removeWorstCell(i1, newCellI1)
            i2_child = self.removeWorstCell(i2, newCellI2)
            newPop.append(i1_child)
            newPop.append(i2_child)
        return newPop

    def epigen_mechanism(self, population):
        """
        This function applies to each cell of the population the epigenetic
        mechanisms with its corresponding probabilities.
        Inputs:
            - population: the total population.
            - mechanisms: A list of mechanisms that we will apply.
            - epiProb: A list of probabilities corresponding to the previous
                    specified mechanisms.
        Output:
            The new modified population.
        """
        for individual in population:
            for cell in individual:
                self.apply_mechanisms(cell)
        return population

    def apply_mechanisms(self, cell):
        """
        This function applies the epigenetic mechanisms to a given cell
        with some probability.
        Already implemented mechanisms:
            - ...
        Possible future implemented mechanisms:
            - "imprinting"
            - "reprogramming"
            - "paramutation"
            - "position"
            - "inactivation"
            - "bookmarking"
            - "silencing"
        Inputs:
            - mechanisms: List of mechanisms to be applied.
            - cell: The cell in which we will apply a mechanism.
            - epiProb: List of probabilities for every listed mechanism.
        Output:
            The new modified cell.
        """
        modified = False
        for i in range(len(self.mechanisms)):
            if random() < self.epi_probs[i]:
                modified = True
                if self.mechanisms[i] == "imprinting":
                    self.imprinting_mechanism(cell, self.imprinting_prob)
                elif self.mechanisms[i] == "reprogramming":
                    # TODO: Hacer reprogramming
                    pass
                elif self.mechanisms[i] == "paramutation":
                    # TODO: Hacer paramutation
                    pass
                elif self.mechanisms[i] == "position":
                    self.position_mechanism(cell, self.position_prob)
                elif self.mechanisms[i] == "inactivation":
                    # TODO: Hacer x-inactivation
                    pass
                elif self.mechanisms[i] == "bookmarking":
                    # TODO: Hacer bookmarking
                    pass
                elif self.mechanisms[i] == "silencing":
                    # TODO: Hacer gene silencing
                    pass
        if modified:
            self.evaluate_cell(cell)
        return cell

    def imprinting_mechanism(self, cell, probability):
        """
        Change the provenience of a gene with a certain probability, 
        then swap-it for the value where it already existed.

        Inputs:
            - Cell: to apply mechanisn
            - probability: probability of applying the mechanism
        Output:
            The new modified cell
        """
        # If there is no father or mother can't do the mechanism
        if cell.mother == None or cell.father == None:
            return cell

        for i in range(len(cell.nucleosome)):

            # If there is no change pass
            if cell.nucleosome[i] == False or random() > probability:
            # if random() > probability: # Use this line to ignore mask in the mechanism application
                continue

            # From wich parent the gene is from?
            parent_equal = cell.father[i]
            parent_change = cell.mother[i]
            if cell.solution[i] != parent_equal:
                parent_equal = parent_change
                parent_change = cell.father[i]

            # Wich is the value to swap to avoid replication
            value_replace = cell.solution[i]
            index_replace = cell.solution.index(parent_change)

            # Change the value from the other parent
            cell.solution[i] = parent_change
            # Swap-it of place where it was from
            cell.solution[index_replace] = value_replace

        return cell

    def position_mechanism(self, cell, probability):
        """
        For the genes that are collapsed change their position within a certain probability.
        The position change does not take into account if the genes are in the same collapsed zone.

        Inputs:
            - Cell: to apply mechanisn
            - probability: probability of applying the mechanism
        Output:
            The new modified cell
        """

        # Get the indexes for genes that will be affected
        affected_indexes = []
        for i in range(len(cell.nucleosome)):
            if cell.nucleosome[i] == True and random() <= probability:
                affected_indexes.append(i)

        relocation = []
        for elem in affected_indexes:
            relocation.append(elem)

        shuffle(relocation)
        newsolution = [a for a in cell.solution]

        # Change positions
        for i in range(len(relocation)):
            newsolution[affected_indexes[i]] = cell.solution[relocation[i]]

        cell.solution = newsolution

        return cell

    def replacement(self, oldpop,  newpop):
        """
        Get the best of the two populations. Pure elitism used.
        TODO: Change method
        """
        newpop = sorted(
            [*oldpop, *newpop], key=lambda x: self.evaluate_individual(x)
        )
        return newpop[:self.individuals_number]

    def on_launch(self):
        # Set up plot
        self.figure, (self.ax, self.ax2) = plt.subplots(1, 2)
        self.lines, = self.ax.plot([], [], 'go-', label="Path found")
        self.lines_optimum, = self.ax.plot(
            [], [], 'ro-.', alpha=0.5, label="Optimal path")
        self.lines2, = self.ax2.plot(
            [], [], 'b.', alpha=0.2, label="Individual fitness")
        self.mean_line, = self.ax2.plot([], [], 'y-', label="Mean fitness")
        self.min_line, = self.ax2.plot([], [], 'g-', label="Best fitness")
        # Autoscale on unknown axis and known lims on the other
        self.ax.set_autoscaley_on(True)
        self.ax2.set_autoscaley_on(True)
        self.ax2.set_yscale("log")
        self.ax2.get_yaxis().get_major_formatter().labelOnlyBase = False
        # Other stuff
        self.ax.legend()
        self.ax.grid()
        self.ax2.legend()
        self.ax2.grid()

    def on_running(self, coordinates, currentPath, optimum_path, title_string):
        # Update data (with the new _and_ the old points)
        xdata = []
        ydata = []
        xdata_opt = []
        ydata_opt = []

        for i in range(len(currentPath)):
            xdata.append(coordinates[currentPath[i]][0])
            ydata.append(coordinates[currentPath[i]][1])
        xdata.append(coordinates[currentPath[0]][0])
        ydata.append(coordinates[currentPath[0]][1])
        self.lines.set_xdata(xdata)
        self.lines.set_ydata(ydata)

        for i in range(len(optimum_path)):
            xdata_opt.append(coordinates[optimum_path[i]][0])
            ydata_opt.append(coordinates[optimum_path[i]][1])
        xdata_opt.append(coordinates[optimum_path[0]][0])
        ydata_opt.append(coordinates[optimum_path[0]][1])
        self.lines_optimum.set_xdata(xdata_opt)
        self.lines_optimum.set_ydata(ydata_opt)
        # Need both of these in order to rescale
        self.ax.set_title(title_string)
        self.ax.relim()
        self.ax.autoscale_view()
        # We need to draw *and* flush
        # plt.title(title_string)
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def on_running_fitness(self, population, iteration, min_fit, title_string):
        # Update data (with the new _and_ the old points)
        total_fitness = 0
        j = 0

        # All fitness points:
        for i in range(len(population)):
            total_fitness += self.evaluate_individual(
                population[i])/self.cells_number
            j += 1
            self.xdata_fitness.append(iteration)
            self.ydata_fitness.append(self.evaluate_individual(
                population[i])/self.cells_number)
        self.lines2.set_xdata(self.xdata_fitness)
        self.lines2.set_ydata(self.ydata_fitness)

        # Best fitness point:
        self.min_fitness.append(min_fit/self.cells_number)
        self.ydata_iter.append(iteration)
        self.min_line.set_ydata(self.min_fitness)
        self.min_line.set_xdata(self.ydata_iter)

        # Mean fitness point:
        self.mean_fitness.append(total_fitness/j)
        self.mean_line.set_ydata(self.mean_fitness)
        self.mean_line.set_xdata(self.ydata_iter)

        # Need both of these in order to rescale
        self.ax2.set_title(title_string)
        self.ax2.relim()
        self.ax2.autoscale_view()
        # We need to draw *and* flush
        # plt.title(title_string)
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def print_metric(self, population, coordinates, optimum_path, iteration):
        fitness = 99999999
        min_cell = None
        #min_ind = max(population, key=lambda x:x.fitness)
        for individual in population:
            pivot_fitness = self.evaluate_individual(individual)
            if pivot_fitness <= fitness:
                fitness = pivot_fitness
                cell_fitness = list(map(lambda cell: cell.fitness, individual))
                #print(f"Estas son las fitness de las celulas:{cell_fitness}")
                min_cell = individual[np.argmin(cell_fitness)]

        # TODO: Hacer la funcion print del mejor individual
        #print(f'Best Individual is: {min_ind}')
        self.on_running(coordinates, min_cell.solution, optimum_path,
                        "Iteration: "+str(iteration) + " Best Path: " + str(int(fitness)))
        self.on_running_fitness(population, iteration, fitness,
                                "Distances of the population")

    def roulette_selection(self, population):
        individual_fitness = [self.evaluate_individual(
            individual) for individual in population]
        minim_fitness = np.min(individual_fitness)
        maxim_fitness = np.max(individual_fitness)
        media = np.floor((minim_fitness+maxim_fitness)/2)
        translated = list(map(lambda x: x-media, individual_fitness))
        inverted = list(map(lambda x: -x, translated))
        reverted = list(map(lambda x: x+media, inverted))
        max = sum(reverted)
        pick = uniform(0, max)
        current = 0
        for i in range(len(population)):
            current += reverted[i]
            if current > pick:
                return population[i]
