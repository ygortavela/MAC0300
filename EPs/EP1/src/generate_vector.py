import sys
import random

if len(sys.argv) != 2:
    print("Utilize o comando 'python3 generate_vector #n' para gerar um vetor de tamanho n")
    exit()

n = int(sys.argv[1])
fd = open('vector-test/vector-' + str(n), 'w')

fd.write(str(n) + '\n')

for i in range(n):
    fd.write(str(random.uniform(-10, 10)))

    if i != n:
        fd.write(' ')

fd.write('\n')

fd.close()
