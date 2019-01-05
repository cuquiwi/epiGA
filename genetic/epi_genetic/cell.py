import numpy as np

class Cell:

    def __init__(self, solution, father = None, mother = None, mask = None):
        self.solution = solution
        self.father = father
        self.mother = mother
        if mask == None:
            self.mask = np.zeros(len(solution))
        else:
            self.mask = mask
