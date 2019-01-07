from random import shuffle
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
        - max_epoch: 
    """
    if i>=max_epoch:
        return True
    else:
        return False

def init_population(individualsNb, cellsNb):
    population = []
    for i in range(individualsNb):
        individual = []
        for j in range(cellsNb):
            solution = shuffle([k+1 for k in range(128)]) #TODO generalize solution for other problems than a280
            individual.append(Cell(solution))

        population.append(individual)

    return population

def selection(population):
    #TODO
    return population

def nucleosome_reproduction(population):
   #TODO
    return population

def generate_nucleosome(population, nucleo_prob, nucleoRad):
    #TODO
    return population

def epigen_mechanism(population, epiProb):
    #TODO
    return population

def replacement(oldpop,  newpop):
    #TODO
    return newpop
