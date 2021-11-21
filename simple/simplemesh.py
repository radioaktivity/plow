import numpy as np

def pretty_print(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))

N = 10
boxsize = 1


# Mesh
dx = boxsize / N
vol = dx**2
xlin = np.linspace(0.5*dx, boxsize-0.5*dx, N)
Y, X = np.meshgrid( xlin, xlin )

rho = np.ones(X.shape)
u = np.zeros(X.shape)
v = np.zeros(X.shape)
p = np.ones(X.shape)*1013*100


print(np.round(u,decimals=3))
print('_'*50)
print(np.round(rho,decimals=3))
