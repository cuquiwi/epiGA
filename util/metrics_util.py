import matplotlib.pyplot as plt
import numpy as np

class GeneticMetricPrinter(object):

    def __init__(self):
        self.xdata_fitness = []
        self.ydata_fitness = []
        self.mean_fitness = []
        self.min_fitness = []
        self.ydata_iter = []

    def on_launch(self):
        plt.ion()
        # Set up plot
        self.figure, (self.ax, self.ax2) = plt.subplots(1, 2)
        self.lines, = self.ax.plot([], [], 'go-', label="Path found")
        self.lines_optimum, = self.ax.plot(
            [], [], 'ro-.', alpha=0.5, label="Optimal path")
        self.lines2, = self.ax2.plot(
            [], [], 'b.', alpha=0.2, label="Individual fitness")
        self.mean_line, = self.ax2.plot([], [], 'y-', label="Mean fitness")
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

    def on_running(self, coordinates, currentPath, optimum_path, title_string):
        # Update data (with the new _and_ the old points)
        xdata = []
        ydata = []
        xdata_opt = []
        ydata_opt = []

        for i in range(len(currentPath)):
            xdata.append(coordinates[currentPath[i]][0])
            ydata.append(coordinates[currentPath[i]][1])
        xdata.append(coordinates[currentPath[0]][0])
        ydata.append(coordinates[currentPath[0]][1])
        self.lines.set_xdata(xdata)
        self.lines.set_ydata(ydata)

        for i in range(len(optimum_path)):
            xdata_opt.append(coordinates[optimum_path[i]][0])
            ydata_opt.append(coordinates[optimum_path[i]][1])
        xdata_opt.append(coordinates[optimum_path[0]][0])
        ydata_opt.append(coordinates[optimum_path[0]][1])
        self.lines_optimum.set_xdata(xdata_opt)
        self.lines_optimum.set_ydata(ydata_opt)
        # Need both of these in order to rescale
        self.ax.set_title(title_string)
        self.ax.relim()
        self.ax.autoscale_view()
        # We need to draw *and* flush
        # plt.title(title_string)
        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def on_running_fitness(self, population, fitnesses, iteration, min_fit, title_string):
        # Update data (with the new _and_ the old points)
        total_fitness = 0
        j = 0

        # All fitness points:
        for i in range(len(population)):
            total_fitness += fitnesses[i]/len(population[0])
            j += 1
            self.xdata_fitness.append(iteration)
            self.ydata_fitness.append(fitnesses[i]/len(population[0]))
        self.lines2.set_xdata(self.xdata_fitness)
        self.lines2.set_ydata(self.ydata_fitness)

        # Best fitness point:
        self.min_fitness.append(min_fit/len(population[0]))
        self.ydata_iter.append(iteration)
        self.min_line.set_ydata(self.min_fitness)
        self.min_line.set_xdata(self.ydata_iter)

        # Mean fitness point:
        self.mean_fitness.append(total_fitness/j)
        self.mean_line.set_ydata(self.mean_fitness)
        self.mean_line.set_xdata(self.ydata_iter)

        # Need both of these in order to rescale
        self.ax2.set_title(title_string)
        self.ax2.relim()
        self.ax2.autoscale_view()

        self.figure.canvas.draw()
        self.figure.canvas.flush_events()

    def on_epoch(self, population, fitnesses, coordinates, optimum_path, iteration):
        fitness = np.Infinity
        min_cell = None
        for individual, pivot_fitness in zip(population, fitnesses):
            if pivot_fitness <= fitness:
                fitness = pivot_fitness
                cell_fitness = list(map(lambda cell: cell.fitness, individual))
                min_cell = individual[np.argmin(cell_fitness)]

        self.on_running(coordinates, min_cell.solution, optimum_path,
                        "Iteration: "+str(iteration) + " Best Path: " + str(int(fitness)))
        self.on_running_fitness(population,fitnesses, iteration, fitness,
                                "Distances of the population")