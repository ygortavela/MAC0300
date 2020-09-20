import sys
import random

if len(sys.argv) != 2:
    print("Utilize o comando 'python3 generate_vector #n' para gerar um vetor de tamanho n")
    exit()

n = int(sys.argv[1])
test = range(n)
fd = open('vector-test/vector-' + str(n), 'w')

fd.write(str(n) + '\n')

for i in test:
    fd.write(str(random.randint(-10, 10)))

    if i != n:
        fd.write(' ')

fd.write('\n')

fd.close()
