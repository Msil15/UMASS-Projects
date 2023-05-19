import numpy as np
import matplotlib.pyplot as plt

def S(x, n):

    term_list = []

    for k in range(1, n):
        term = (2/(k*np.pi))*(1 - np.cos(np.pi*k))*np.sin(k*x)
        term_list.append(term)

    S_n = sum(term_list)

    print(S_n)

    return S_n

input_points = np.linspace(-np.pi, np.pi, 200)

list_0 = []
list_1 = []
list_2 = []
list_3 = []

for x in input_points:
    output_0 = S(x, 0)
    list_0.append(output_0)

    output_1 = S(x, 1)
    list_1.append(output_1)

    output_2 = S(x, 2)
    list_2.append(output_2)

    output_3 = S(x, 3)
    list_3.append(output_3)

plt.hlines(y= -1, xmin=-np.pi, xmax=0, label = 'f(x)')
plt.hlines(y = 1, xmin=0, xmax=np.pi)
#plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.plot(input_points, list_0, 'y', label = '$S_{0}$')
plt.plot(input_points, list_1, 'r-.', label = '$S_{1}$')
plt.plot(input_points, list_2, 'c', label = '$S_{2}$')
plt.plot(input_points, list_3, 'm--', label = '$S_{3}$')
plt.legend()
plt.title('Increasing convergance for trig polynomials on a piecewise function')
plt.xlabel('x')
plt.ylabel('y')
plt.show()