import numpy as np

class Cell:

    def __init__(self, solution, father = None, mother = None, nucleosome = [], fitness = None):
        self.solution = solution
        self.father = father
        self.mother = mother
        self.fitness = fitness
        if len(nucleosome)== 0:
            self.nucleosome = np.zeros(len(solution), dtype=bool)
        else:
            self.nucleosome = nucleosome

    def setfitness(self, fitness):
        self.fitness = fitness

    #TODO: Hacer una funcion __str__