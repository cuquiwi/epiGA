import numpy as np

class Cell:

    def __init__(self, solution, father = None, mother = None, nucleosome = None, fitness = None):
        self.solution = solution
        self.father = father
        self.mother = mother
        self.fitness = fitness
        if nucleosome == None:
            self.nucleosome = np.zeros(len(solution))
        else:
            self.nucleosome = nucleosome
