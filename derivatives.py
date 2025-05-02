from libpyanalytical import get_partitions, partition_multiplicity, Polynomial, PolyFrac, PolyExp

N = 5
partitions = get_partitions(N)
print(f'Partitions of {N} :')
for p in partitions:
    print(f'{partition_multiplicity(p):<5} : {p}')

print('\n------------------------')
print('PolyExp and PolyFrac :\n')
p1 = Polynomial(-2, 5, [1, 2, 3, -1, -2, -3, 2.5, 3.1], 1)
p2 = Polynomial(-10, -4, [1, -2, 3, -4], 2)

pexp = PolyExp(p1, p2) # p1(x) * exp[p2(x)]
pfrac = PolyFrac(p1, p2) # p1(x) / p2(x)

x = 2.3145
dx = 1e-8
for n in range(6):
    print(f'Derivative {n} : {pexp.derivative(x, n)} / {pexp.numerical_derivative(x, n, 1e-8)}')
    print(f'Derivative {n} : {pfrac.derivative(x, n)} / {pfrac.numerical_derivative(x, n, 1e-8)}\n')

