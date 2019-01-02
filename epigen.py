def epigen_alg(individualsNb, cellsNb, epiProb, nucleoProb, nucleoRad, max_epoch = 500):
    population = init_population(individualsNb, cellsNb)
    aux_population = []
    i = 0
    while i < max_epoch: #TODO add termination condition
        newpop = selection(population[i])
        newpop = nucleosome_reproduction(newpop)
        newpop = epigen_mechanism(newpop, epiProb)

        aux_population.append(newpop)
        population.append(replacement(poulation[i], newpop))
        i = i+1

def selection(population):
    #TODO
    return population

def generate_nucleosome(population, nucleo_prob, nucleoRad):
    #TODO
    return population

def epigen_mechanism(population, epiProb):
    #TODO
    return population

def replacement(oldpop,  newpop)
    #TODO
    return newpop
