import numpy as np
import matplotlib.pyplot as plt

Bisection_List = []
Newton_List = []
Secant_List = []
Chord_List = []
False_Position_List = []
actual_zero = -1.27

def Bisection(a, b, TOL, N0, g):
    i = 1
    
    f_a = g(a) #defining the function at the bounds
    f_b = g(b)

    if f_a == 0: #checking if the bounds are roots
        return a
    if f_b == 0:
        return b

    if np.sign(f_a) == np.sign(f_b): #checking if function crosses the x-axis by comparing bound signs.
        return 'Function does not cross x-axis on this interval.'

    while i <= N0: 
        p = a + (b-a)/2 #defining p and f(p)
        f_p = g(p) 
        if f_p == 0 or (b-a)/2 < TOL: #checking if f(p) is a root or outside tolerance.
            Bisection_List.append(abs(p - actual_zero))
            print('Finished in ', i, 'iterations')
            return p
        else:
            i += 1
            Bisection_List.append(abs(p - actual_zero))
            f_check = f_a * f_p #if root is not found at bounds or p, loop through, narrowing region
            if f_check > 0:
                a = p
                f_a = f_p
            else:
                b = p

    print('Method failed after ', N0, 'iterations')

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
            Newton_List.append(abs(p - actual_zero))
            p0 = p

    print('Method failed after ', N0, 'iterations')

def Secant_Method(a, b, p0, p1, TOL, N0, g):
    i = 2

    f_a = g(a) #defining the function at the bounds
    f_b = g(b)

    if f_a == 0: #checking if the bounds are roots
        return a
    if f_b == 0:
        return b

    while i <= N0:
       q0 = g(p0)
       q1 = g(p1)

       p = p1 - q1*(p1-p0)/(q1-q0)

       if abs(p - p1) < TOL:
           Secant_List.append(abs(p - actual_zero))
           print('Finished in ', i, 'iterations')
           return p
       else:
           i += 1
           Secant_List.append(abs(p - actual_zero))
           p0 = p1
           q0 = q1
           p1 = p
           q1 = g(p)

    print('Method failed after ', N0, 'iterations')

def Chord_Method(a, b, p0, TOL, N0, g):
    i = 1

    f_a = g(a) #defining the function at the bounds
    f_b = g(b)

    if f_a == 0: #checking if the bounds are roots
        return a
    if f_b == 0:
        return b

    while i <= N0:
        p = p0 - (g(p0)*(b-a))/(g(b) - g(a))

        if abs(p - p0) < TOL:
            Chord_List.append(abs(p - actual_zero))
            print('Finished in ', i, 'iterations')
            return p
        else:
            i += 1
            Chord_List.append(abs(p - actual_zero))
            p0 = p

    print('Method failed after ', N0, 'iterations')

def False_Position_Method(a, b, p0, p1, TOL, N0, g):
    i = 2

    f_a = g(a) #defining the function at the bounds
    f_b = g(b)

    if f_a == 0: #checking if the bounds are roots
        return a
    if f_b == 0:
        return b

    while i <= N0:
       q0 = g(p0)
       q1 = g(p1)

       p = p1 - q1*(p1-p0)/(q1-q0)

       if abs(p - p1) < TOL:
           False_Position_List.append(abs(p - actual_zero))
           print('Finished in ', i, 'iterations')
           return p
       else:
           i += 1
           False_Position_List.append(abs(p - actual_zero))
           q = g(p)

       if q*q1 < 0:
           p0 = p1
           q0 = q1
       
       p1 = p
       q1 = q

    print('Method failed after ', N0, 'iterations')

def myfzeros(a, b, TOL, p0, p1, N0, g, gprime, method):
    if method == 'b':
        root = Bisection(a, b, TOL, N0, g)
        print('This root was found with the bisection method!')
        return root
    elif method == 'n':
        root = Newtons_Method(a, b, p0, TOL, N0, g, gprime)
        print("This root was found with Newton's method!")
        return root
    elif method == 'c':
        root = Chord_Method(a, b, p0, TOL, N0, g)
        print('This root was found with the chord method!')
        return root
    elif method == 's':
        root = Secant_Method(a, b, p0, p1, TOL, N0, g)
        print('This root was found with the secant method!')
        return root
    elif method == 'f':
        root = False_Position_Method(a, b, p0, p1, TOL, N0, g)
        print('This root was found with the false position method!')
        return root

def g(x):
    return x**9 + x + 10

def g_prime(x):
    return 9*x**8 + 1

root_list = ['b', 'n', 'c', 's', 'f']

for root in root_list:
    root_check = myfzeros(-1.5, -1, 10e-10, -0.9, -1.1, 300, g, g_prime, root)
    print('the root is', root_check)

B_n_list = [*range(1, (len(Bisection_List) +1))]
N_n_list = [*range(1, (len(Newton_List) +1))]
C_n_list = [*range(1, (len(Chord_List) +1))]
S_n_list = [*range(1, (len(Secant_List) +1))]
F_n_list = [*range(1, (len(False_Position_List) +1))]
plt.plot(B_n_list, Bisection_List, '-o', label = 'Bisection')
plt.plot(N_n_list, Newton_List, '-*', label = 'Newton')
plt.plot(C_n_list, Chord_List, '-v', label = 'Chord')
plt.plot(S_n_list, Secant_List, '-s', label = 'Secant')
plt.plot(F_n_list, False_Position_List, '-d', label = 'False Position')
plt.legend()
plt.yscale('log')
plt.title('Comparison of Error Convergence', fontsize = 16)
plt.xticks(fontsize = 13)
plt.yticks(fontsize = 13)
plt.ylabel('Absolute Error', fontsize = 14)
plt.xlabel('Number of iterations', fontsize = 14)
plt.show()