import os
from glob import glob

sequential = "sequential"
omp = "omp"
thread = "thread"
fastflow = "fastflow"

matrix_size = "matrix_size:"
algorithm = "algorithm:"
nworkers = "workers:"

iterations_computed = "iterations_computed:"
error_s = "error:"
computation_time_s = "computation_time:"
init_time_s = "initialization_time:"

directory = "./data"
filetype = "*.txt"


def get_metrics(filename):
    pass


def main():
    os.chdir(directory)
    for filename in glob(filetype):
        get_metrics(filename)


if __name__ == "__main__":
    main()
