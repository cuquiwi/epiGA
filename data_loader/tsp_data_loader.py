import numpy as np
import re


def load_problem_file(file_name):
    """Loads a data file containg a TSP problem defined by iwr in their own
       format

    Arguments:
        file_name {str} -- Path to the file
    Returns:
        matrix NxN, matrix Nx2 -- Symetric matrix with the distances among cities, matrix containing the coordinates.
    """
    coordinates = []
    with open(file_name, 'r') as file:
        for line in file.readlines()[6:-2]:
            raw_coordinates = re.findall(r'\d+', line)[1:]
            coordinates.append(
                np.array([
                    int(raw_coordinates[0]),
                    int(raw_coordinates[1])
                ])
            )

    return coordinates


def load_solution_file(file_name):
    """Loads a solution file containing a TSP problem defined by iwr in their
    own format

    Arguments:
        file_name {str} -- Path to the file
    Returns:
        vector -- Ordered cities
    """

    path = []
    with open(file_name, 'r') as file:
        for line in file.readlines()[4:-2]:
            city_index = int(re.findall(r'\d+', line)[0]) - 1
            path.append(city_index)

    return path
