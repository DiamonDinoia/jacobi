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
total_time = "total time"
iterations = "iterations computed"
error = "error"

directory = "./data"
filetype = "*.csv"

labels = [algorithm, matrix_size, workers, init_time, comp_time, iterations, error]


def line_graph(frame, label, xticklabel, ylabel, xlabel, figname):
    col = [x for x in labels if label in x]
    if init_time in col:
        frame[total_time] = frame[col[0]] + frame[col[1]]
        col.append(total_time)
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
    a.set_yticks(frame[total_time])
    locator = ticker.AutoLocator()
    a.yaxis.set_major_locator(locator)
    plt.legend(col, loc='best')
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.title(figname)
    plt.savefig(figname)
    plt.clf()
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


def line_metric_graph(frame, label, xticklabel, ylabel, xlabel, figname):
    frame[label].plot(kind='line', linewidth=2.0, rot=0, grid=True)
    if label == 'efficiency':
        plt.ylim([0, 1])
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.title(figname if 'xxl' not in figname else 'efficiency')
    plt.savefig(figname)
    plt.clf()
    # plt.show()


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

    time = sequential_data['size'][sequential_data['size']['matrix_size'] == 1024][total_time].iloc[0]

    frames = [omp_data, thread_data, ff_data]
    for frame in frames:
        frame[workers]['ideal time'] = [time / x for x in frame[workers][workers]]
        frame[workers]['efficiency'] = frame[workers]['ideal time'] / frame[workers][total_time]
    line_metric_graph(omp_data[workers], 'efficiency', workers, 'efficiency', workers, 'omp efficiency')
    line_metric_graph(thread_data[workers], 'efficiency', workers, 'efficiency', workers, 'thread efficiency')
    line_metric_graph(ff_data[workers], 'efficiency', workers, 'efficiency', workers, 'fastflow efficiency')

    time = sequential_data['size'][sequential_data['size']['matrix_size'] == 16384][total_time].iloc[0]

    for frame in frames:
        frame['xxl']['ideal time'] = [time / x for x in frame['xxl'][workers]]
        frame['xxl']['efficiency'] = frame['xxl']['ideal time'] / frame['xxl'][total_time]
        # print frame
    line_metric_graph(omp_data['xxl'], 'efficiency', workers, 'efficiency', workers, 'omp xxl efficiency')
    line_metric_graph(thread_data['xxl'], 'efficiency', workers, 'efficiency', workers, 'thread xxl efficiency')
    line_metric_graph(ff_data['xxl'], 'efficiency', workers, 'efficiency', workers, 'fastflow xxl efficiency')


if __name__ == "__main__":
    main()
