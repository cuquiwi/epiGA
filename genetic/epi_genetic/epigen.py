from random import shuffle, random
from cell import Cell
import numpy as np


distMatrix = [[1]] #TODO: Hacer global la matriz de distancias y llamarla asi.

def epigen_alg(problemMatrix, individualsNb, cellsNb, epiProb, nucleoProb, nucleoRad, mechanisms, environment, max_epoch = 500):
    """
    EpiGA based on the work by D.H. Stolfi and E. Alba, (2017).
    Inputs:
        problemMatrix: Matrix with the distances for each city
        individualsNb: Number of individuals in the problem.
        cellsNb: Number of cells per each individual.
        epiProb: Epigenetic probabilities list.
        nucleoProb: Nucleosome probability.
        nucleoRad: Nucleosome radius.
        mechanisms: List of string containing the epigenetic algorithms
                    in the same order as the epiProb list.
        environment: Environment rules
    """
    population = init_population(individualsNb, cellsNb, len(problemMatrix))
    aux_population = []
    i = 0
    termination_condition = False
    while not termination_condition:
        newpop = selection(population[i])
        newpop = nucleosome_generation(newpop, nucleoProb, nucleoRad)
        newpop = nucleosome_reproduction(newpop)
        newpop = epigen_mechanism(newpop, mechanisms, epiProb)

        aux_population.append(newpop)
        population.append(replacement(population[i], newpop))
        termination_condition = termination(i, max_epoch)
        i = i+1

def termination(i, max_epoch):
    """
    Termination condition for the EpiGA.
    Inputs:
        - i: Current iteration of the algorithm.
        - max_epoch: Maximum number of epochs.
    Output:
        True if the termination condition is accomplished, False
        otherwise.
    """
    #TODO: Hacer una super funcion de terminacion. Discuss why.
    if i>=max_epoch:
        return True
    else:
        return False

def evaluate_cell(cell, distMatrix):
    """
    Function that evaluates the fitness of our problem for a given cell
    and sets the cell fitness.
    Inputs:
        - cell: a cell which contains the solution.
        - distMatrix: Matrix with the distance
    Outputs:
        The fitness value of a given cell.
    """
    solution = cell.solution
    fitness = 0
    for i in range(1,len(solution)):
        fitness += distMatrix[i-1][i]
    cell.setfitness(fitness)
    
    return fitness

def init_population(individualsNb, cellsNb, distMatrix):
    """
    The initial population of the EpiGA.
    Inputs:
        - individualsNb: Number of initial individuals in the problem.
        - cellsNb: Number of cells per each individual.
        - distMatrix: Matrix with the distance
    Return:
        A list of cell elements representing the initial population.
    """
    population = []
    for i in range(individualsNb):
        individual = []
        for j in range(cellsNb):
            solution = [k+1 for k in range(len(distMatrix))]
            shuffle(solution)
            cell = Cell(solution)
            evaluate_cell(cell, distMatrix)
            individual.append(cell)

        population.append(individual)

    return population

def selection(population):
    """
    Performs the population selection in the EpiGA.
    Inputs:
        - population: the total population.
    Return:
        THe selected subset of the total population.
    """
    #TODO: Binary tournment or roulette. Discuss why we have chosen a particular implementation.
    return population

def nucleosome_generation(population, prob, radius):
    """
    Generates a new nucleosome vector as a mask for each cell in the 
    individuals of the population.
    Inputs:
        - population: the total population.
        - prob: Nucleosome generation probability.
        - radius: Radius of the nucleosome.
    Output:
        The population with the new nucleosomes generated.
    """
    for i in range(len(population)):
        individual = population[i]
        for j in range(len(individual)):
            cells = individual[j]
            n = cells.mask
            for k in range(len(n)):
                if random() < prob:
                    n = collapse(n,radius,k)
            cells.mask = n
            individual[j] = cells
        population[i] = individual
    return population

def collapse(nucleosome, radius, k):
    """
    Function that creates a new nucleosome.
    Inputs:
        - nucleosome: The binary list of the cell that represents the nucleosome.
        - radius: trivial.
        - k: the center position of the new nucleosome.
    Output:
        The new modified nucleosome.
    """
    for i in range(-radius, radius+1,1):
        nucleosome[k+i] = 1
    return nucleosome

def selectBestCell(individual):
    #TODO:  Add description
    fitness = list(map(lambda cell:cell.fitness, individual))
    return individual(np.argmax(fitness))

def crossover(baseSolution, secondSolution, mask):
    #TODO:  Add description
    # TODO:  Implement Partially-mapped Crossover (PMX) 
    # https://www.researchgate.net/publication/226665831_Genetic_Algorithms_for_the_Travelling_Salesman_Problem_A_Review_of_Representations_and_Operators
    return baseSolution

def removeWorstCell(individual, newCell):
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
    fitness = list(map(lambda cell:cell.fitness, individual))
    newInd.remove(individual[np.argmin(fitness)])
    newInd.append(newCell)
    return newInd
    

def nucleosome_reproduction(population):
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
                bestCell1 = selectBestCell(i1)
                bestCell2 = selectBestCell(i2)
                newNucleosome = np.logical_or(bestCell1.nucleosome, bestCell2.nucleosome)
                fatherBasedSolution = crossover(bestCell1.solution, bestCell2.solution, newNucleosome)
                motherBasedSolution = crossover(bestCell2.solution, bestCell1.solution, newNucleosome)
                newCellI1 = Cell(fatherBasedSolution, bestCell1.solution,bestCell2.solution,newNucleosome)
                newCellI2 = Cell(motherBasedSolution, bestCell2.solution,bestCell1.solution,newNucleosome)
                i1_child = removeWorstCell(i1, newCellI1)
                i2_child = removeWorstCell(i2, newCellI2)
                newPop.append(i1_child, i2_child)
    return newPop

def epigen_mechanism(population, mechanisms, epiProb):
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
            cells = individual[j]
            cells = apply_mechanisms(mechanisms, cells, epiProb)
            evaluate_cell(cells, distMatrix)
            individual[j] = cells
        population[i] = individual
    return population

def apply_mechanisms(mechanisms, cell, epiProb):
    """
    This function applies the epigenetic mechanisms to a given cell
    with some probability.
    Already implemented mechanisms:
        - ...
    Possible future implemented mechanisms:
        - "imprinting"
        - "reprogramming"
        - "permutation"
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
    for i in range(len(mechanisms)):
        if random() < epiProb[i]:
            if mechanisms[i] == "imprinting":
                #TODO: Hacer gen imprimting
                pass
            elif mechanisms[i] == "reprogramming":
                #TODO: Hacer reprogramming
                pass
            elif mechanisms[i] == "permutation":
                #TODO: Hacer gen imprimting
                pass
            elif mechanisms[i] == "position":
                #TODO: Hacer position effect
                pass
            elif mechanisms[i] == "inactivation":
                #TODO: Hacer x-inactivation
                pass
            elif mechanisms[i] == "bookmarking":
                #TODO: Hacer bookmarking
                pass
            elif mechanisms[i] == "silencing":
                #TODO: Hacer gene silencing
                pass
    return cell

def replacement(oldpop,  newpop):
    #TODO
    return newpop

