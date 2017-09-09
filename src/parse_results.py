import matplotlib.pyplot as plt
import sys
from optparse import OptionParser

import numpy as np
import pandas as pd

serial_jacobi = 'serial jacobi'
parallel_for = 'parallel for'
thread_jacobi = 'thread jacobi'
map_jacobi = 'map jacobi'

algorithms = [serial_jacobi, parallel_for, thread_jacobi, map_jacobi]

iterations_computed = 'iterations computed'
error = 'error'
solution = 'solution'
total_time = 'total time'
computation_time = 'computation time'

params = [iterations_computed, solution, total_time, error, computation_time]
run = 'RUN'
workers = 'workers'

data = {}

correct_solution = {}


def get_correct_solution(name):
    with open(name, 'r') as f:
        correct_solution[name] = []
        found = False
        for line in f:
            if solution + ':' in line:
                found = True
                continue
            if found:
                correct_solution[name].append(np.float64(line))


def read_file(infile):
    with open(infile, 'r') as f:
        name = ''
        for line in f:
            line_split = line.split(' ')
            line_split[-1] = line_split[-1][:-1]
            if workers in line:
                nwokers = int(line_split[-1])
            if run in line:
                name = line_split[-1]
                get_correct_solution(name)
                data[name] = {}
                for algorithm in algorithms:
                    data[name][algorithm] = {}
            for algorithm in algorithms:
                for param in params:
                    if algorithm in line and param in line:
                        if param == error:
                            data[name][algorithm][error] = np.float64(line_split[-1])
                            data[name][algorithm][iterations_computed] = np.float64(line_split[-3])
                        elif param == solution:
                            data[name][algorithm][param] = [np.float64(x) for x in line_split[5:-1]]
                        else:
                            data[name][algorithm][param] = np.float64(line_split[-1])
    return nwokers


metrics = {}

time = 'time'
speedup = 'speedup'
ideal_speedup = 'ideal speedup'
sequential = 'sequential'
_map = 'map'
thread = 'thread'
parallel_for = 'parallel for'
ideal_time = 'ideal time'
overhead = 'overhead'
matrix_size = 'matrix size'
efficiency = 'efficiency'


def check_error(a, b):
    a = np.array(a, dtype=np.float64)
    b = np.array(b, dtype=np.float64)
    c = a - b
    c = np.abs(c)
    c = np.sum(c)
    return c


labels = []


def calculate_metrics(nworkers):
    for matrix, v in data.items():
        key = matrix
        metrics[matrix] = {}
        metrics[matrix][sequential + ' ' + time] = v[serial_jacobi][total_time]
        metrics[matrix][sequential + ' ' + error] = check_error(v[serial_jacobi][solution], correct_solution[matrix])

        metrics[matrix][matrix_size] = int(len(correct_solution[matrix]))
        metrics[matrix][ideal_speedup] = int(nworkers)
        metrics[matrix][ideal_time] = metrics[matrix][sequential + ' ' + time] / nworkers

        metrics[matrix][thread + ' ' + time] = v[thread_jacobi][total_time]
        metrics[matrix][thread + ' ' + overhead] = v[thread_jacobi][total_time] - v[thread_jacobi][computation_time]
        metrics[matrix][thread + ' ' + speedup] = metrics[matrix][sequential + ' ' + time] / v[thread_jacobi][
            total_time]
        metrics[matrix][thread + ' ' + error] = check_error(v[thread_jacobi][solution], correct_solution[matrix])
        metrics[matrix][thread + ' ' + efficiency] = metrics[matrix][ideal_time] / v[thread_jacobi][total_time]

        metrics[matrix][parallel_for + ' ' + time] = v[parallel_for][total_time]
        metrics[matrix][parallel_for + ' ' + overhead] = v[parallel_for][total_time] - v[parallel_for][computation_time]
        metrics[matrix][parallel_for + ' ' + speedup] = metrics[matrix][sequential + ' ' + time] / v[parallel_for][
            total_time]
        metrics[matrix][parallel_for + ' ' + error] = check_error(v[parallel_for][solution], correct_solution[matrix])
        metrics[matrix][parallel_for + ' ' + efficiency] = metrics[matrix][ideal_time] / v[parallel_for][total_time]

        metrics[matrix][_map + ' ' + time] = v[map_jacobi][total_time]
        metrics[matrix][_map + ' ' + overhead] = v[map_jacobi][total_time] - v[map_jacobi][computation_time]
        metrics[matrix][_map + ' ' + speedup] = metrics[matrix][sequential + ' ' + time] / v[map_jacobi][total_time]
        metrics[matrix][_map + ' ' + error] = check_error(v[map_jacobi][solution], correct_solution[matrix])
        metrics[matrix][_map + ' ' + efficiency] = metrics[matrix][ideal_time] / v[map_jacobi][total_time]

    for x in metrics[key].keys():
        labels.append(x)
    return pd.DataFrame(metrics)


def write_metrics(outfile):
    with open(outfile, 'w') as f:
        for matrix, v in data.items():
            f.write(matrix + '\n')
            for key, val in sorted(metrics[matrix].items()):
                f.write(key + ": " + str(val) + '\n')
            f.write('\n')


def write_dataframe(outfile, frame):
    with open(outfile, 'w') as f:
        frame.to_string(f)


def write_tolatex(outfile, frame):
    with open(outfile, 'w') as f:
        f.write(frame.to_latex())


sizes = ['4', '64', '512', '1024', '2048', '4096']


def bar_graph_time(frame):
    col = [x for x in labels if 'time' in x]
    frame[col].plot(kind='bar', rot=0)
    col = [x.replace(' time', '') for x in col]
    plt.legend(col, loc='best')
    plt.ylabel('Time (s)')
    plt.xlabel('Matrix size')
    plt.savefig('time')
    # plt.show()


def bar_graph_parallel_benchmark(frame):
    col = [x for x in labels if 'time' in x and sequential not in x]
    frame[col].plot(kind='bar', rot=0)
    col = [x.replace(' time', '') for x in col]
    plt.legend(col, loc='best')
    plt.ylabel('Time (s)')
    plt.xlabel('Matrix size')
    plt.savefig('benchmark')
    # plt.show()


def speedup_graph(frame):
    col = [x for x in labels if 'speedup' in x]
    x = frame[col].plot(kind='line', rot=0, grid=True, linewidth=3.0, color=['blue', 'green', 'red', 'black'])
    col = [x.replace(' speedup', '') for x in col]
    plt.legend(col, loc='best')
    plt.ylabel('Speedup')
    plt.xlabel('Matrix size')
    plt.savefig('speedup')
    # plt.show()


def efficiency_graph(frame):
    col = [x for x in labels if 'efficiency' in x]
    frame[col].plot(kind='line', rot=0, grid=True, linewidth=3.0, color=['blue', 'green', 'red', 'black'])
    col = [x.replace(' efficiency', '') for x in col]
    plt.legend(col, loc='best')
    plt.ylabel('Efficiency')
    plt.xlabel('Matrix size')
    axes = plt.gca()
    axes.set_ylim([0, 1])
    plt.savefig('efficiency')
    # plt.show()


def main(argv):
    parser = OptionParser()
    parser.add_option('-i', dest='infile')
    parser.add_option('-o', dest='outfile')
    options, _ = parser.parse_args(argv)
    nworkers = read_file(options.infile)
    frame = calculate_metrics(nworkers)
    write_dataframe(options.outfile, frame)
    write_tolatex('latex' + options.outfile, frame)
    frame.columns = sizes
    frame = frame.transpose()
    bar_graph_time(frame)
    bar_graph_parallel_benchmark(frame)
    speedup_graph(frame)
    efficiency_graph(frame)
    # print frame.transpose()
    return 0


if __name__ == "__main__":
    main(sys.argv[1:])
