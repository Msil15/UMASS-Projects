import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return x*np.log(x)

a = 1
b = 3

def trans_domain(a, b, x):
    t_x = (1/2)*((b - a)*x + a + b)
    return t_x

def T3(x):
    return x**3 - (3/4)*x

n = 3

x_k_list = []
x_k_t_list = []
y_k_t_list = []

for k in range(1, n + 1):
    x_k = np.cos(((2*k - 1)/(2*n))*np.pi)
    x_k_list.append(x_k)

    xkt = trans_domain(a, b, x_k)
    x_k_t_list.append(xkt)

    ykt = func(xkt)
    y_k_t_list.append(ykt)

print(x_k_list)
print(x_k_t_list)
print(y_k_t_list)

#LAGRANGE INTERPOLATION
def Lagrange(x):
    term01 = (x - x_k_t_list[1])/(x_k_t_list[0] - x_k_t_list[1])
    term02 = (x - x_k_t_list[2])/(x_k_t_list[0] - x_k_t_list[2])

    L0 = term01*term02

    term11 = (x - x_k_t_list[0])/(x_k_t_list[1] - x_k_t_list[0])
    term12 = (x - x_k_t_list[2])/(x_k_t_list[1] - x_k_t_list[2])

    L1 = term11*term12

    term21 = (x - x_k_t_list[0])/(x_k_t_list[2] - x_k_t_list[0])
    term22 = (x - x_k_t_list[1])/(x_k_t_list[2] - x_k_t_list[1])

    L2 = term21*term22

    LT = y_k_t_list[0]*L0 + y_k_t_list[1]*L1 + y_k_t_list[2]*L2

    return LT

z = np.linspace(a, b, 22)
L_value = []

for x_val in z:
    L = Lagrange(x_val)
    L_value.append(L)

plt.plot(z, L_value, 'b--')
plt.plot(x_k_t_list, y_k_t_list, 'ro')
plt.title('Lagrange & Chebyshev Approximation for f(x) = xln(x)')
plt.xlabel('x')
plt.ylabel('y')
plt.show()