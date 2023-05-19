import numpy as np
import matplotlib.pyplot as plt

def y_actual(t):
    return (1/(2*(t**2)))*(4 + np.cos(2) - np.cos(2*t))

def y_prime(t, y):
    return (t**(-2))*(np.sin(2*t) - 2*t*y)

omega_list = []
t_list = []

def Forward_Euler(a, b, alph, h, func):
    N = (b - a)/h
    t = a
    omega = alph
    
    for i in np.arange(N + 1):
       omega_list.append(omega)
       t_list.append(t)
       omega = omega + h*func(t, omega)
       t = t + h

    return

Forward_Euler(1, 2, 2, 0.25, y_prime)

y_short_list = []

for t in t_list:
    y = y_actual(t)
    y_short_list.append(y)

y_actual_list = []
t_actual_list = []

for t in np.linspace(1, 2, 25):
    y = y_actual(t)
    y_actual_list.append(y)
    t_actual_list.append(t)


plt.plot(t_list, omega_list, 'bp--', label = 'Forward Euler')
plt.plot(t_actual_list, y_actual_list, 'r', label = 'y(t)')
plt.legend()
plt.xlabel('t')
plt.ylabel('y')
plt.show()
plt.clf()

error_list = []
error_bound_list = []
M = 7.53
L = 2
h = 0.25
a = 1

for y, omega, t in zip(y_short_list, omega_list, t_list):
    error = abs(y - omega)
    error_list.append(error)
    error_bound = h*(M/(2*L))*(np.exp(L*(t - a) - 1))
    error_bound_list.append(error_bound)

plt.plot(t_list, error_list, 'bp--', label = 'real error')
plt.plot(t_list, error_bound_list, 'ro-', label = 'theoretical bound')
plt.xlabel('t')
plt.ylabel('Absolute error')
plt.legend()
plt.show()