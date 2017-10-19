import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
from glob import glob

import pandas as pd

sequential = "sequential"
omp = "omp"
thread = "thread"
fastflow = "fastflow"

matrix_size = "matrix_size"
algorithm = "algorithm"
workers = "workers"
init_time = "initialization time"
comp_time = "computation time"
iterations = "iterations computed"
error = "error"

directory = "./data"
filetype = "*.csv"

labels = [algorithm, matrix_size, workers, init_time, comp_time, iterations, error]


def line_graph(frame, label, xticklabel, ylabel, xlabel, figname):
    col = [x for x in labels if label in x]
    if init_time in col:
        frame['total time'] = frame[col[0]] + frame[col[1]]
        col.append('total time')
    a = frame[col].plot(kind='line', linewidth=2.0, rot=0, grid=True)
    if workers == xticklabel:
        a.set_xticks(frame[xticklabel])
        locator = ticker.AutoLocator()
        a.xaxis.set_major_locator(locator)
        ticks = a.get_xticks()
        if ticks[0] == 0:
            ticks[0] = 1
        a.set_xticks(ticks)
        plt.xlim(1, len(frame[xticklabel]))
    else:
        a.set_xticklabels(frame[xticklabel])
    a.set_yticks(frame['total time'])
    locator = ticker.AutoLocator()
    a.yaxis.set_major_locator(locator)
    plt.legend(col, loc='best')
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.title(figname)
    plt.savefig(figname)
    # plt.show()


omp_data = {}
sequential_data = {}
thread_data = {}
ff_data = {}

metrics = ['iterations', 'size', 'workers', 'xxl']


def get_data(metric, filename, label, xticklabel, ylabel, xlabel, figname):
    frame = pd.read_csv(filename)
    if omp in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join([omp, figname]))
        omp_data[metric] = frame
    if thread in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join([thread, figname]))
        thread_data[metric] = frame
    if fastflow in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join([fastflow, figname]))
        ff_data[metric] = frame
    if sequential in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join([sequential, figname]))
        sequential_data[metric] = frame


def main():
    os.chdir(directory)
    for filename in glob(filetype):
        if "workers" in filename:
            metric = 'workers'
            get_data(metric, filename, 'time', metric, metric, 'time', 'workers with small matrix')
        elif 'iterations' in filename:
            metric = 'iterations'
            get_data(metric, filename, 'time', iterations, metric, 'time', 'different iterations')
        elif 'size' in filename:
            metric = 'size'
            get_data(metric, filename, 'time', matrix_size, 'matrix size', 'time', 'matrix different sizes')
        else:
            metric = 'workers'
            get_data('xxl', filename, 'time', metric, metric, 'time', 'workers with large matrix')


if __name__ == "__main__":
    main()
