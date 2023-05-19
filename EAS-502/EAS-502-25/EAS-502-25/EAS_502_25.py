import numpy as np
import matplotlib.pyplot as plt

omega_2_list = []
omega_4_list = []
t_list = []

def y_actual(t):
    return (1/(2*(t**2)))*(4 + np.cos(2) - np.cos(2*t))

def f(t, y):
    return (t**(-2))*(np.sin(2*t) - 2*t*y)

def f_prime(t, y):
    return (-2*(t**(-3)))*(np.sin(2*t) - 2*t*y) + (t**(-2))*(2*np.cos(2*t) - 2*(y + t*f(t, y)))

def f_double_prime(t, y):
    return (6/(t**4))*(np.sin(2*t) - 2*t*y) + (-2*(t**(-3)))*(2*np.cos(2*t) - 2*(y + t*f(t, y))) + (-2*(t**(-3)))*(2*np.cos(2*t) - 2*(y + t*f(t, y))) + (t**(-2))*(-4*np.sin(2*t) - 2*(f(t, y) + f(t, y) + t*f_prime(t, y)))

def f_triple_prime(t, y):
    a = (24/(t**(-5)))*(np.sin(2*t) - 2*t*y) + (6/(t**4))*(2*np.cos(2*t) - 2*(y + t*f(t, y)))
    b = (6/(t**4))*(2*np.cos(2*t) - 2*(y + t*f(t, y))) + (-2*(t**(-3)))*(-4*np.sin(2*t) - 2*(f(t, y) + f(t, y) + t*f_prime(t, y)))
    c = (6/(t**4))*(2*np.cos(2*t) - 2*(y + t*f(t, y))) + (-2*(t**(-3)))*(-4*np.sin(2*t) - 2*(f(t, y) + f(t, y) + t*f_prime(t, y)))
    d = (-2*(t**(-3)))*(-4*np.sin(2*t) - 2*(f(t, y) + f(t, y) + t*f_prime(t, y))) + (t**(-2))*(-8*np.cos(2*t) -2*(f_prime(t, y) + f_prime(t, y) + f_prime(t, y) + t*f_double_prime(t, y)))
    return a + b + c + d

def Taylor(a, b, alph, h, order):
    N = (b - a)/h
    t = a

    omega = alph

    if order == 2:
        for i in np.arange(N + 1):
            T2 = f(t, omega) + (h/2)*f_prime(t, omega)
            omega_2_list.append(omega)
            t_list.append(t)
            omega = omega + h*T2
            t = t + h
    elif order == 4:
        for i in np.arange(N + 1):
            T4 = f(t, omega) + (h/2)*f_prime(t, omega) + ((h**2)/6)*f_double_prime(t, omega) + ((h**3)/24)*f_triple_prime(t, omega)
            omega_4_list.append(omega)
            t_list.append(t)
            omega = omega + (h)*T4
            t = t + h
    else:
        print('Not an accepted order')
        return

Taylor(1, 2, 2, 0.25, 2)
t_list = []
Taylor(1, 2, 2, 0.25, 4)

y_actual_list = []

print(omega_2_list)

for t in t_list:
    y = y_actual(t)
    y_actual_list.append(y)

plt.plot(t_list, omega_2_list, 'bp--', label = 'Taylor 2')
plt.plot(t_list, omega_4_list, 'gh-.', label = 'Taylor 4')
plt.plot(t_list, y_actual_list, 'ro-', label = 'y(t)')
plt.legend()
plt.show()
plt.clf()

o_2_error_list = []
o_4_error_list = []

for o_2, o_4, y in zip(omega_2_list, omega_4_list, y_actual_list):
    o_2_error = abs(y - o_2)
    o_4_error = abs(y - o_4)
    o_2_error_list.append(o_2_error)
    o_4_error_list.append(o_4_error)

plt.plot(t_list, o_2_error_list, 'bp--', label = 'Taylor 2')
plt.plot(t_list, o_4_error_list, 'gh-.', label = 'Taylor 4')
plt.legend()
plt.show()