import numpy as np


def generate_matrix(size, name):
    matrix = np.random.rand(size, size)
    terms = np.array([x[0] for x in np.random.rand(size, 1)])
    for i in range(size):
        value = sum(map(abs, matrix[i])) - matrix[i][i]
        while value > abs(matrix[i][i]):
            matrix[i][i] += 1.0
    solution = np.linalg.solve(matrix, terms)
    with open(name, 'w') as f:
        f.write('size: {0} \n'.format(size))
        f.write('matrix:\n')
        np.savetxt(f, matrix)
        f.write('terms:\n')
        np.savetxt(f, terms)
        f.write('solution:\n')
        np.savetxt(f, solution)


i = 0
sizes = [4, 64, 512, 1024, 2048, 4096]

for size in sizes:
    generate_matrix(size, 'matrix' + str(i) + '.txt')
    i += 1
    print 'generated matrix number: ', i, 'size: ', size
