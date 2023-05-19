import numpy as np
import matplotlib.pyplot as plt

X_n = [*range(1, 5, 1)]

n = len(X_n)

L = [[0 for n in range(n)]for j in range(n)]

n_list = []

L_list = []

for n in X_n:
    n_range = np.linspace(-1, 1, 1 + n)
    n_list.append(n_range)

print(n_list)


def n_func(n):
    return (2**(n+1))/(n*np.log(n))

other_n_values = [2, 3, 4, 5, 6]

func_list = []

def func(n):
    return (2**(n+1))/(n*np.log(n))

for n in other_n_values:
    y = func(n)
    func_list.append(y)

plt.plot(other_n_values, func_list, 'ro')
plt.show()