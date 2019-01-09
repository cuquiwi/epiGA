from random import shuffle, random
from cell import Cell
import numpy as np

class EpigeneticAlgorithm(object):

    def __init__(self, individuals_number, cells_number, epi_probs,
                 nucleo_prob, nucleo_rad, mechanisms, environment, max_epochs=500):
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
        """
        self.individuals_number = individuals_number
        self.cells_number = cells_number
        self.epi_probs = epi_probs
        self.nucleo_prob = nucleo_prob
        self.nucleo_rad = nucleo_rad
        self.mechanisms = mechanisms
        self.enviroment = environment
        self.max_epochs = max_epochs

    def call(self, distances_matrix):
        """Performs the actual algorithm taking into account the configuration 
        provided in the initialization and the distances matrix provided here
        
        Arguments:
            distances_matrix {Matrix NxN} -- Matrix containing the distances 
                                        between cities
        """

        self.distances_matrix = distances_matrix

        population = self.init_population()
        aux_population = []
        i = 0
        termination_condition = False
        while not termination_condition:
            newpop = self.selection(population[i])
            newpop = self.nucleosome_generation(newpop)
            newpop = self.nucleosome_reproduction(newpop)
            newpop = self.epigen_mechanism(newpop)

            aux_population.append(newpop)
            population.append(self.replacement(population[i], newpop))
            termination_condition = self.termination(i)
            i = i+1

    def termination(self, i):
        """
        Termination condition for the EpiGA.
        Inputs:
            - i: Current iteration of the algorithm.
        Output:
            True if the termination condition is accomplished, False
            otherwise.
        """
        # TODO: Hacer una super funcion de terminacion. Discuss why.
        if i >= self.max_epochs:
            return True
        else:
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
            fitness += self.distances_matrix[i-1][i]
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
        fitness = 0
        for cell in individual:
            fitness += cell.fitness
        return fitness

    def init_population(self):
        """
        The initial population of the EpiGA.
        Return:
            A list of cell elements representing the initial population.
        """
        population = []
        for i in range(self.individuals_number):
            individual = []
            for j in range(self.cells_number):
                solution = [k+1 for k in range(len(self.distances_matrix))]
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
                n = cells.mask
                for k in range(len(n)):
                    if random() < self.nucleo_prob:
                        n = self.collapse(n, k)
                cells.mask = n
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
            nucleosome[k+i] = 1
        return nucleosome

    def selectBestCell(self, individual):
        # TODO:  Add description
        fitness = list(map(lambda cell: cell.fitness, individual))
        return individual(np.argmax(fitness))

    def crossover(self, baseSolution, secondSolution, mask):
        # TODO:  Add description
        # TODO:  Implement Partially-mapped Crossover (PMX)
        # https://www.researchgate.net/publication/226665831_Genetic_Algorithms_for_the_Travelling_Salesman_Problem_A_Review_of_Representations_and_Operators
        return baseSolution

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
        newInd.remove(individual[np.argmin(fitness)])
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
        for i1 in population:
            for i2 in population:
                if (not (i1 == i2)):
                    bestCell1 = self.selectBestCell(i1)
                    bestCell2 = self.selectBestCell(i2)
                    newNucleosome = np.logical_or(
                        bestCell1.nucleosome, bestCell2.nucleosome)
                    fatherBasedSolution = self.crossover(
                        bestCell1.solution, bestCell2.solution, newNucleosome)
                    motherBasedSolution = self.crossover(
                        bestCell2.solution, bestCell1.solution, newNucleosome)
                    newCellI1 = Cell(
                        fatherBasedSolution, bestCell1.solution, bestCell2.solution, newNucleosome)
                    newCellI2 = Cell(
                        motherBasedSolution, bestCell2.solution, bestCell1.solution, newNucleosome)
                    i1_child = self.removeWorstCell(i1, newCellI1)
                    i2_child = self.removeWorstCell(i2, newCellI2)
                    newPop.append(i1_child, i2_child)
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
        for i in range(len(population)):
            individual = population[i]
            for j in range(len(individual)):
                cell = individual[j]
                cell = self.apply_mechanisms(cell)
                self.evaluate_cell(cell)
                #TODO: Creo que esta línea está demas porque todos los objetos en python se pasan por referencia. 
                individual[j] = cell
                #####################
            population[i] = individual
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
        for i in range(len(self.mechanisms)):
            if random() < self.epi_probs[i]:
                if self.mechanisms[i] == "imprinting":
                    #TODO: Hacer gen imprimting
                    pass
                elif self.mechanisms[i] == "reprogramming":
                    #TODO: Hacer reprogramming
                    pass
                elif self.mechanisms[i] == "paramutation":
                    #TODO: Hacer paramutation
                    pass
                elif self.mechanisms[i] == "position":
                    self.position_mechanism(cell, 0.4) # TODO: parametize the prob?
                elif self.mechanisms[i] == "inactivation":
                    #TODO: Hacer x-inactivation
                    pass
                elif self.mechanisms[i] == "bookmarking":
                    #TODO: Hacer bookmarking
                    pass
                elif self.mechanisms[i] == "silencing":
                    #TODO: Hacer gene silencing
                    pass
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
        
        #Get the indexes for genes that will be affected
        affected_indexes = []
        for i in range(cell.nucleosome):
            if cell.nucleosome[i] == 1 & random() <= probability:
                affected_indexes.append(i)

        relocation = affected_indexes # TODO: Does it copy reference or value?
        shuffle(relocation)
        newsolution = cell.solution

        #Change positions
        for i in range(len(relocation)):
            newsolution[affected_indexes[i]] = cell.solution[relocation[i]]
        cell.solution = newsolution

        return cell

    def replacement(self, oldpop,  newpop):
        #TODO
        return newpop
