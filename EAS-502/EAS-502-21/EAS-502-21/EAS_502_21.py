import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.exp(x)*np.cos(2*x)

x_points = np.linspace(-np.pi, np.pi, 8)
y_points = []

for x in x_points:
    y = func(x)
    y_points.append(y)

def a_k(x_list, y_list, k):

    term_list = []

    for x, y, in zip(x_list, y_list):
        term = y*np.cos(k*x)
        term_list.append(term)

    a = (1/4)*sum(term_list)
    return a

def b_k(x_list, y_list, k):

    term_list = []

    for x, y, in zip(x_list, y_list):
        term = y*np.sin(k*x)
        term_list.append(term)

    b = (1/4)*sum(term_list)
    return b

a0 = a_k(x_points, y_points, 0)
a1 = a_k(x_points, y_points, 1)
a2 = a_k(x_points, y_points, 2)
a3 = a_k(x_points, y_points, 3)

b1 = b_k(x_points, y_points, 1)
b2 = b_k(x_points, y_points, 2)

def S3(x):
    S = a0/2 + a3*np.cos(3*x) + a1*np.cos(x) + b1*np.sin(x) + a2*np.cos(2*x) + b2*np.sin(2*x)
    return S

test_points = np.linspace(-np.pi, np.pi, 1000)
function_points = []
S_points = []


for x in test_points:
    func_y = func(x)
    S_y = S3(x)

    function_points.append(func_y)
    S_points.append(S_y)
plt.plot(test_points, function_points)
plt.plot(test_points, S_points, '--')
plt.show()

term_list = []

for func_y, S_y in zip(function_points, S_points):

    term = (func_y - S_y)**2
    term_list.append(term)

E = sum(term_list)
print(E)
