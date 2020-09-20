import sys
import random

if len(sys.argv) != 2:
    print("Utilize o comando 'python3 generate_matrix #n' para gerar uma matriz n x n")
    exit()

n = int(sys.argv[1])
test = range(n)
fd = open('matrix-test/matrix-' + str(n), 'w')

fd.write(str(n) + '\n')

for i in test:
    for j in test:
        fd.write(str(random.randint(-10, 10)))

        if j != n:
            fd.write(' ')

    fd.write('\n')

fd.close()
