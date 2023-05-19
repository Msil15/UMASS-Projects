import numpy as np
import matplotlib.pyplot as plt

Newton_List = []

actual_zero = 0

def g(x):
    return np.exp(x) - x + 1

def g_prime(x):
    return np.exp(x) - 1

def Newtons_Method(a, b, p0, TOL, N0, g, gprime):
    i = 1

    f_a = g(a) #defining the function at the bounds
    f_b = g(b)

    if f_a == 0: #checking if the bounds are roots
        return a
    if f_b == 0:
        return b

    while i <= N0:
        p = p0 - g(p0)/gprime(p0)

        if abs(p - p0) < TOL:
            Newton_List.append(abs(p - actual_zero))
            print('Finished in ', i, 'iterations')
            return p
        else:
            i += 1
            Newton_List.append(p)
            p0 = p

    print('Method failed after ', N0, 'iterations')

Newton = Newtons_Method(-0.5, 0.5, 1, 10e-10, 20, g, g_prime)

n = 0

lambda_list = []

while n <= 18:
    lambda_value = abs(Newton_List[n+1] - actual_zero)/abs(Newton_List[n] - actual_zero)
    lambda_list.append(lambda_value)
    n += 1

print(lambda_list)

avg = np.mean(lambda_list)
n_list = []

print(avg)
