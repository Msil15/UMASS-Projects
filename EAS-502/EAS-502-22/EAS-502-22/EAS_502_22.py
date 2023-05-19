import numpy as np
import matplotlib.pyplot as plt

a = -1
b = 1

m = 1000
p = 3

def func(x):
    return x

x_list = []

for j in range(2*m):
    x = -np.pi + j*np.pi/m
    x_list.append(x)

y_list = x_list

'''
def FFT(m, p, y_list):
    
    c_list = np.zeros_like(y_list)
    eta_list = np.zeros_like([*range(m + 1)])
    
    M = m
    q = p
    zeta = np.exp(np.pi*1j/m)

    for j in range(2*m + 1):
        c_list[j] = y_list[j]

    for j in range(M + 1):
        eta_list[j] = zeta**j
        eta_list[j + M] = -eta_list[j]

    K = 0
    eta_list[0] = 1

    for L in range(1, p + 2):
        while K < 2*m - 1:
            for j in range(M + 1):
'''

c_k_list = []

for k in range(2*m):
    c_k_term_list = []
    for j in range(2*m):
        c_k_term = y_list[j]*np.exp(1j*k*np.pi*j/4)
        c_k_term_list.append(c_k_term)

    ck_sum = sum(c_k_term_list)
    c_k_list.append(ck_sum)

print(c_k_list)

ck_mag_list = []
k_list = range(0, 2*m)

for ck in c_k_list:
    mag = np.sqrt(ck.real**2 + ck.imag**2)
    ck_mag_list.append(mag)

plt.plot(k_list, ck_mag_list, 'ro--')
plt.ylabel('ck')
plt.xlabel('k')
plt.show()