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


def load_solution_file(file_name, coordinates):
    """Loads a solution file containing a TSP problem defined by iwr in their
    own format

    Arguments:
        file_name {str} -- Path to the file
        coordinates {List} -- Coordinates
    Returns:
        vector -- Ordered cities
    """

    path = []
    with open(file_name, 'r') as file:
        for line in file.readlines()[4:-2]:
            city_index = int(re.findall(r'\d+', line)[0]) - 1
            path.append(city_index)

    objective_distance = 0.0
    for i in range(len(path) - 1):
        c_from = path[i]
        c_to = path[i + 1]
        objective_distance += np.linalg.norm(
            coordinates[c_to] - coordinates[c_from]
        )

    return path, objective_distance
