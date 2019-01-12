import matplotlib.pyplot as plt
import numpy as np


class GeneticMetricPrinter(object):

    def __init__(self):
        self.xdata_fitness = []
        self.ydata_fitness = []
        self.mean_fitness = []
        self.mean_fitness_cells = []
        self.min_fitness = []
        self.ydata_iter = []
        self.coordinates = []
        self.optimum_path = []

    def on_launch(self, coordinates, optimum_path, epigenetic=True):
        """
        Inputs:
            - epigenetic: if we are going to plot an EpiGA.
        """
        self.coordinates = coordinates
        self.optimum_path = optimum_path

        plt.ion()
        # Set up plot
        self.figure, (self.ax, self.ax2) = plt.subplots(1, 2)
        self.lines, = self.ax.plot([], [], 'go-', label="Path found")
        self.lines_optimum, = self.ax.plot(
            [], [], 'ro-.', alpha=0.5, label="Optimal path")
        self.lines2, = self.ax2.plot(
            [], [], 'b.', alpha=0.2, label="Individual fitness")
        self.mean_line, = self.ax2.plot([], [], 'm-', label="Mean fitness")
        if epigenetic:
            self.mean_line_best, = self.ax2.plot([], [], 'c-', label="Mean fitness individuals")
        self.min_line, = self.ax2.plot([], [], 'g-', label="Best fitness")
        # Autoscale on unknown axis and known lims on the other
        self.ax.set_autoscaley_on(True)
        self.ax2.set_autoscaley_on(True)
        self.ax2.set_yscale("log")
        self.ax2.get_yaxis().get_major_formatter().labelOnlyBase = False
        # Other stuff
        self.ax.legend()
        self.ax.grid()
        self.ax2.legend()
        self.ax2.grid()

    def on_epoch(self, solutions, fitnesses, iteration, epigenetic=True, minimum = None):
        """Receives solutions sorted by their fitness
        
        Arguments:
            solutions {[type]} -- [description]
            fitnesses {[type]} -- [description]
            iteration {[type]} -- [description]
        """
        best_solution = solutions[0]
        best_fitness = fitnesses[0]

        self.on_running(best_solution,
                        "Iteration: "+str(iteration) + " Best Path: " + str(int(best_fitness)))
        self.on_running_fitness(solutions, fitnesses, iteration, best_fitness,
                                "Distances of the population", epigenetic, minimum)

    def on_running(self, currentPath, title_string):
        # Update data (with the new _and_ the old points)
        xdata = []
        ydata = []
        xdata_opt = []
        ydata_opt = []

        for i in range(len(currentPath)):
            xdata.append(self.coordinates[currentPath[i]][0])
            ydata.append(self.coordinates[currentPath[i]][1])
        xdata.append(self.coordinates[currentPath[0]][0])
        ydata.append(self.coordinates[currentPath[0]][1])
        self.lines.set_xdata(xdata)
        self.lines.set_ydata(ydata)

        for i in range(len(self.optimum_path)):
            xdata_opt.append(self.coordinates[self.optimum_path[i]][0])
            ydata_opt.append(self.coordinates[self.optimum_path[i]][1])
        xdata_opt.append(self.coordinates[self.optimum_path[0]][0])
        ydata_opt.append(self.coordinates[self.optimum_path[0]][1])
        self.lines_optimum.set_xdata(xdata_opt)
        self.lines_optimum.set_ydata(ydata_opt)
        # Need both of these in order to rescale
        self.ax.set_title(title_string)
        self.ax.relim()
        self.ax.autoscale_view()
        # We need to draw *and* flush
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def on_running_fitness(self, solutions, fitnesses, iteration, min_fit, title_string, epigenetic, minimum):
        # Update data (with the new _and_ the old points)
        total_fitness = 0

        # All fitness points:

        self.xdata_fitness.extend([iteration for _ in range(len(solutions))])
        self.ydata_fitness.extend(fitnesses)
        total_fitness = sum(fitnesses)

        self.lines2.set_xdata(self.xdata_fitness)
        self.lines2.set_ydata(self.ydata_fitness)

        # Best fitness point:
        self.min_fitness.append(min_fit)
        self.ydata_iter.append(iteration)
        self.min_line.set_ydata(self.min_fitness)
        self.min_line.set_xdata(self.ydata_iter)

        # Mean fitness point:
        self.mean_fitness.append(total_fitness/len(solutions))
        self.mean_line.set_ydata(self.mean_fitness)
        self.mean_line.set_xdata(self.ydata_iter)

        # Mean fitness point for best cells in an EpiGA:
        if epigenetic:
            self.mean_fitness_cells.append(np.mean(minimum))
            self.mean_line_best.set_ydata(self.mean_fitness_cells)
            self.mean_line_best.set_xdata(self.ydata_iter)


        # Need both of these in order to rescale
        self.ax2.set_title(title_string)
        self.ax2.relim()
        self.ax2.autoscale_view()

        self.figure.canvas.draw()
        self.figure.canvas.flush_events()
