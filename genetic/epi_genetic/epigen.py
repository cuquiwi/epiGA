from random import shuffle, random
from cell import Cell

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

def evaluate_individual(individual):
    """
    Get the total fitness of the individual, the sum of the fitness of its cells.
    Input:
        - individual: Individual, represented by a list of cells
    Return:
        - fitness of the individual
    """
    fitness = 0
    for cell in individual:
        fitness+=cell.fitness
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
    Uses a binary tournament selection with no repetitions.
    Returns half of the population
    Inputs:
        - population: the total population.
    Return:
        The selected subset of the total population.
    """
    winners = []
    shuffle(population)
    i=0
    while i+1<len(population):
        if evaluate_individual(population[i]) > evaluate_individual(population[i+1]):
            winners.append(population[i])
        else:
            winners.append(population[i+1])
        i+=2
    return winners

def k_tournament_selection(pop, k=2):
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
        if (best == None) or evaluate_individual(ind) > bestfit:
            best = ind
            bestfit = evaluate_individual(best)
    return best

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
            individual[j] = cells
        population[i] = individual
    return population

def apply_mechanisms(mechanisms, cells, epiProb):
    #TODO: Hacer los mecanismos
    return cells

def replacement(oldpop,  newpop):
    #TODO
    return newpop

