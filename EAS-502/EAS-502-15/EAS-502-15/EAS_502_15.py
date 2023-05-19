import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def function(x):
    return abs(x)*np.exp(x)

exact_value = quad(function, 0, 1)[0]

print('%.16f ' % exact_value)

def comp_simp(a, b, n, func):
    h = (b - a)/n

    XI0 = func(a) + func(b)
    XI1 = 0
    XI2 = 0

    for i in range(n):
        X = a + i*h
        if (i % 2) == 0:
            XI2 = XI2 + func(X)
        else:
            XI1 = XI1 + func(X)

    XI = h*(XI0 + 2*XI2 + 4*XI1)/3

    return XI

############################################

def comp_trap(a, b, n, func):
    h = (b - a)/n

    XI0 = func(a) + func(b)
    XI1 = 0

    for i in range(n):
        X = a + i*h
        XI1 = XI1 + func(X)

    XI = h*(XI0 + 2*XI1)/2

    return XI

############################################

def comp_mid(a, b, n, func):
    h = (b - a)/(n + 2)

    XI1 = 0

    for i in range(n):
        X = a + i*h
        XI1 = XI1 + func(X)

    XI = 2*h*X

    return X

simp_list = []
trap_list = []
mid_list = []
n_list = [n for n in range(1, 2001) if n % 2 == 0]


for n in n_list:
    simp_result = comp_simp(0, 1, n, function)
    trap_result = comp_trap(0, 1, n, function)
    mid_result = comp_mid(0, 1, n, function)

    simp_error = abs(exact_value - simp_result)
    trap_error = abs(exact_value - trap_result)
    mid_error = abs(exact_value - mid_result)

    simp_list.append(simp_error)
    trap_list.append(trap_error)
    mid_list.append(mid_error)

plt.plot(n_list, simp_list, '-', label = 'Composite Simpson', )
plt.plot(n_list, trap_list, '--', label = 'Composite Trapazoid', )
plt.plot(n_list, mid_list, '-.', label = 'Composite Midpoint')
plt.xlabel('Number of iterations')
plt.ylabel('Absolute error')
plt.legend()
plt.yscale('log')
plt.show()

print(min(simp_list))
print(min(trap_list))
print(min(mid_list))