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


def line_graph(frame, label, xticklabel, ylabel, xlabel, figname, ideal=None):
    col = [x for x in labels if label in x]
    if label == total_time:
        col = [total_time]
    if init_time in col:
        frame[total_time] = frame[col[0]] + frame[col[1]]
        col.append(total_time)
    if ideal is None:
        a = frame[col].plot(kind='line', linewidth=2.0, rot=0, grid=True)
    else:
        b = ideal[col].plot(kind='line', linewidth=2.0, rot=0, grid=True)
        a = frame[col].plot(kind='line', linewidth=2.0, rot=0, grid=True, ax=b)
        # print(ideal[col],frame[col])
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
    tmp = figname.split(' ')
    if ideal is None:
        plt.legend(col, loc='best')
    else:
        plt.legend(['Model', tmp[0]], loc='best')
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    plt.title(tmp[0])
    plt.savefig(figname)
    print(figname)
    plt.clf()
    plt.close()
    # plt.show()


omp_data = {}
sequential_data = {}
thread_data = {}
ff_data = {}

metrics = ['iterations', 'size', 'workers', 'xxl']


def get_data(metric, filename, label, xticklabel, ylabel, xlabel, figname):
    frame = pd.read_csv(filename)
    if omp in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join(['OpenMP', figname]))
        omp_data[metric] = frame
    if thread in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join(['Threads', figname]))
        thread_data[metric] = frame
    if fastflow in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join(['FastFlow', figname]))
        ff_data[metric] = frame
    if sequential in filename:
        line_graph(frame, label, xticklabel, ylabel, xlabel, ' '.join(['Sequential', figname]))
        sequential_data[metric] = frame


def line_metric_graph(frame, label, xticklabel, ylabel, xlabel, figname):
    a = frame[label].plot(kind='line', linewidth=2.0, rot=0, grid=True)
    plt.ylim(0, 1)
    a.set_xticks(frame[xticklabel])
    locator = ticker.AutoLocator()
    a.xaxis.set_major_locator(locator)
    ticks = a.get_xticks()
    if ticks[0] == 0:
        ticks[0] = 1
    a.set_xticks(ticks)
    plt.xlim(1, len(frame[xticklabel]))
    plt.ylabel(xlabel)
    plt.xlabel(ylabel)
    tmp = figname.split(' ')
    plt.title(tmp[0])
    plt.savefig(' '.join(tmp))
    print(figname)
    plt.clf()
    plt.close()
    # plt.show()


def evaluate_model(k, m, n, mult, thread_cost=0., alloc_cost=0.):
    frame = pd.DataFrame(columns=[matrix_size, workers, init_time, comp_time, total_time, iterations])

    for i in range(1, n + 1):
        init_cost = thread_cost * i + alloc_cost * m
        calc = k * mult * m * m / i
        values = [m, i, init_cost, calc, init_cost + calc, k]
        frame.loc[i] = values
    return frame


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

    multiplication_cost = sequential_data['size'][sequential_data['size'][matrix_size] == 16384][comp_time].iloc[0] / (
            50 * 16384 * 16384)
    thread_cost = thread_data['xxl'][thread_data['xxl']['workers'] == 128][init_time].iloc[0]
    thread_cost -= thread_data['xxl'][thread_data['xxl']['workers'] == 1][init_time].iloc[0]
    thread_cost /= 128
    alloc_cost = (thread_data['workers'][thread_data['workers']['workers'] == 1][init_time].iloc[
                      0] - thread_cost) / 1024

    time = sequential_data['size'][sequential_data['size']['matrix_size'] == 1024][comp_time].iloc[0]
    init = sequential_data['size'][sequential_data['size']['matrix_size'] == 1024][init_time].iloc[0]

    idealframe = evaluate_model(50, 1024, 128, multiplication_cost, 0.382702 / 1024, 0.050056 / 1638400)
    line_graph(idealframe, 'time', 'workers', 'workers', 'time', 'Model time 1024')

    frames = [omp_data, thread_data, ff_data]
    for frame in frames:
        frame[workers]['ideal time'] = [(time / x) + init for x in frame[workers][workers]]
        frame[workers]['efficiency'] = frame[workers]['ideal time'] / frame[workers][total_time]

    line_graph(omp_data[workers], total_time, 'workers', 'workers', 'time', 'OpenMP comparision', idealframe)
    line_graph(thread_data[workers], total_time, 'workers', 'workers', 'time', 'Threads comparision', idealframe)
    line_graph(ff_data[workers], total_time, 'workers', 'workers', 'time', 'FastFlow comparision', idealframe)
    line_metric_graph(omp_data[workers], 'efficiency', workers, 'efficiency', workers, 'OpenMP efficiency')
    line_metric_graph(thread_data[workers], 'efficiency', workers, 'efficiency', workers, 'Threads efficiency')
    line_metric_graph(ff_data[workers], 'efficiency', workers, 'efficiency', workers, 'FastFlow efficiency')

    time = sequential_data['size'][sequential_data['size']['matrix_size'] == 16384][comp_time].iloc[0]
    init = sequential_data['size'][sequential_data['size']['matrix_size'] == 16384][init_time].iloc[0]
    idealframe = evaluate_model(50, 16384, 128, multiplication_cost, 0.382702 / 1024, 0.050056 / 1638400)
    line_graph(idealframe, 'time', 'workers', 'workers', 'time', 'Model time 16384')

    for frame in frames:
        frame['xxl']['ideal time'] = [(time / x) + init for x in frame['xxl'][workers]]
        frame['xxl']['efficiency'] = frame['xxl']['ideal time'] / frame['xxl'][total_time]
        # print frame
    line_graph(omp_data['xxl'], total_time, 'workers', 'workers', 'time', 'OpenMP xxl comparision', idealframe)
    line_graph(thread_data['xxl'], total_time, 'workers', 'workers', 'time', 'Threads xxl comparision', idealframe)
    line_graph(ff_data['xxl'], total_time, 'workers', 'workers', 'time', 'FastFlow xxl comparision', idealframe)
    line_metric_graph(omp_data['xxl'], 'efficiency', workers, 'efficiency', workers, 'OpenMP xxl efficiency')
    line_metric_graph(thread_data['xxl'], 'efficiency', workers, 'efficiency', workers, 'Threads xxl efficiency')
    line_metric_graph(ff_data['xxl'], 'efficiency', workers, 'efficiency', workers, 'FastFlow xxl efficiency')


if __name__ == "__main__":
    main()
