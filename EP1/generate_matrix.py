import sys
import random

if len(sys.argv) != 2:
    print("Utilize o comando 'python3 generate_matrix #n' para gerar uma matriz n x n")
    exit()

n = int(sys.argv[1])
fd = open('matrix-test/matrix-' + str(n), 'w')

fd.write(str(n) + '\n')

for i in range(n):
    for j in range(n):
        fd.write(str(random.uniform(-100000, 100000)))

        if j != n:
            fd.write(' ')

    fd.write('\n')

fd.close()
