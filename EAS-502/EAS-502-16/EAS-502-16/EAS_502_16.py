import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def function(x):
    return np.exp(x)

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

def comp_trap(a, b, n, func):
    h = (b - a)/n

    XI0 = func(a) + func(b)
    XI1 = 0

    for i in range(n):
        X = a + i*h
        XI1 = XI1 + func(X)

    XI = h*(XI0 + 2*XI1)/2

    return XI

m_list = [m for m in range(1, 2005) if m % 2 == 0]

for m in m_list:
    simp_result = comp_simp(0, 1, m, function)
    simp_error = abs(exact_value - simp_result)

    error = 5e-4

    if abs(simp_error) <= error:
        print('The desired error was achived after ', m, ' iterations.')
        break

for m in m_list:
    trap_result = comp_trap(0, 1, m, function)
    trap_error = abs(exact_value - trap_result)

    if abs(trap_error) <= 5e-4:
        print('The desired error was achived after ', m, ' iterations.')
        break