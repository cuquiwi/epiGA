from random import shuffle, random
from cell import Cell

def epigen_alg(individualsNb, cellsNb, epiProb, nucleoProb, nucleoRad, mechanisms, environment, max_epoch = 500):
    """
    EpiGA based on the work by D.H. Stolfi and E. Alba, (2017).
    Inputs:
        individualsNb: Number of individuals in the problem.
        cellsNb: Number of cells per each individual.
        epiProb: Epigenetic probabilities list.
        nucleoProb: Nucleosome probability.
        nucleoRad: Nucleosome radius.
        mechanisms: List of string containing the epigenetic algorithms
                    in the same order as the epiProb list.
        environment: Environment rules
    """
    population = init_population(individualsNb, cellsNb)
    aux_population = []
    i = 0
    termination_condition = False
    while not termination_condition:
        newpop = selection(population[i])
        newpop = nucleosome_generation(newpop, nucleoProb, nucleoRad)
        newpop = nucleosome_reproduction(newpop)
        newpop = epigen_mechanism(newpop, epiProb)

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

def evaluate_cell(solution):
    """
    Function that evaluates the fitness of our problem for a given cell.
    Inputs:
        - solution: a cell which contains the solution.
    Outputs:
        The fitness value of a given cell.
    """
    #TODO: Hacer la funcion que evalua el problema
    return 3.1415926535

def init_population(individualsNb, cellsNb):
    """
    The initial population of the EpiGA.
    Inputs:
        - individualsNb: Number of initial individuals in the problem.
        - cellsNb: Number of cells per each individual.
    Return:
        A list of cell elements representing the initial population.
    """
    population = []
    for i in range(individualsNb):
        individual = []
        for j in range(cellsNb):
            solution = shuffle([k+1 for k in range(128)]) #TODO generalize solution for other problems than a280
            cell = Cell(solution)
            evaluate_cell(cell)
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

def nucleosome_reproduction(population):
   #TODO
    return population

def epigen_mechanism(population, epiProb):
    #TODO
    return population

def replacement(oldpop,  newpop):
    #TODO
    return newpop

