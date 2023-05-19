import random
import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return np.sin(2*np.pi*x)

x_values = np.linspace(-1, 1, 22)

y_values = []
y_polluted_values = []
noise = []

for x in x_values:
    y = func(x)
    y_values.append(y)

y_values[0] = y_values[-1] = 0

i = 0

while i <= len(y_values):
    n = 9.5e-4*random.uniform(-1, 1)
    noise.append(n)
    i += 1

for n, y in zip(noise, y_values):
    y_polluted = y + n
    y_polluted_values.append(y_polluted)

z_values = np.linspace(-1, 1, 1000)

def Newton_DD(nodes, values, z):
    n = len(nodes)
    m = len(z)

    a = np.zeros_like(nodes)
    p = np.zeros_like(z)

    F = [[0 for j in range(n)]for i in range(n)]

    for i in range(n):
        F[i][0] = values[i]

    for i in range(1, n):
        for j in range(1, i + 1):
            F[i][j] = (F[i][j-1] - F[i-1][j-1])/(nodes[i] - nodes[i - j])

    for i in range(n):
        a[i] = F[i][i]

    for i in range(m):
        p[i] = a[-1]*(z[i] - nodes[-2]) + a[-2]

    for j in range(2, n -1)[::-1]:
        for i in range(m):
            p[i] = p[i]*(z[i] - nodes[j - 1]) + a[j - 1]

    return p 

interpol_points = Newton_DD(x_values, y_values, z_values)

interpol_points_noise = Newton_DD(x_values, y_polluted_values, z_values)

rel_error_list = []

for x, y in zip(z_values, interpol_points):
    v_o = y #observed value
    v_e = np.sin(x) #expected value

    rel_error = abs((v_o - v_e)/v_e) * 100 #relative percentage error
    rel_error_list.append(rel_error)

rel_error_noise_list = []

for x, y in zip(z_values, interpol_points_noise):
    v_o = y
    v_e = np.sin(x)

    rel_error_noise = abs((v_o - v_e)/v_e) * 100
    rel_error_noise_list.append(rel_error_noise)

rel_error_avg = np.mean(rel_error_list)

rel_error_noise_avg = np.mean(rel_error_noise_list)

print('Relative error is, ', rel_error_avg)
print('Relative noisy error is, ', rel_error_noise_avg)

plt.plot(x_values, y_values, 'ro', label = 'Nodes')
plt.plot(z_values, interpol_points, '--r', label = 'Nodes Interpolated')
plt.ylim(-1.1, 0.3)
plt.xlim(0.4, 1.1)
plt.plot(x_values, y_polluted_values, 'g*', label = 'Polluted Nodes')
plt.plot(z_values, interpol_points_noise, ':b', label = 'Polluted Nodes Interpolated')
plt.xticks(fontsize = 13)
plt.yticks(fontsize = 13)
plt.title('Comparison of Interpolated Data (Insert)', fontsize = 14)
plt.xlabel('x', fontsize = 14)
plt.ylabel('y', fontsize = 14)
#plt.legend()
plt.show()
