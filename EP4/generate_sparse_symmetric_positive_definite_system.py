import sys
import random
import numpy as np
from sklearn.datasets import make_sparse_spd_matrix

if len(sys.argv) != 2:
    print("Use the command 'python3 generate_sparse_symmetric_positive_definite_system #n' to generate a n-dimensional Ax=b system with A being a sparse symmetric positive definite coefficient matrix")
    exit()

n = int(sys.argv[1])
fd = open('testfiles/system-' + str(n), 'w')

fd.write(str(n) + '\n')

# generating coefficients matrix
prng = np.random.RandomState(13)
prec = make_sparse_spd_matrix(n, alpha=.98, random_state=prng)
np.savetxt(fd, prec, fmt="%lf")

# generating b random vector
for i in range(n):
    fd.write(str(random.uniform(0, 1)))

    if i != n:
        fd.write(' ')

fd.write('\n')

fd.close()