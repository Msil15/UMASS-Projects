import numpy as np
import time

'''
def function(x):
    return np.sin(x)

def Bisection(a, b, TOL, N0, func):
    i = 1
    
    f_a = func(a) #defining the function at the bounds
    f_b = func(b)

    if f_a == 0: #checking if the bounds are roots
        return a
    if f_b == 0:
        return b

    if np.sign(f_a) == np.sign(f_b): #checking if function crosses the x-axis by comparing bound signs.
        return 'Function does not cross x-axis on this interval.'

    while i <= N0: 
        p = a + (b-a)/2 #defining p and f(p)
        f_p = func(p) 
        if f_p == 0 or p < TOL: #checking if f(p) is a root or outside tolerance.
            return p
        else:
            i += 1
            f_check = f_a * f_p #if root is not found at bounds or p, loop through, narrowing region
            if f_check < 0:
                a = p
            else:
                b = p

root = Bisection(-np.pi/2, np.pi/2, 10e-3, 10, function)

print(root)
'''

def function_g(x): #defining function g to be passed to FPM
    return x - x**5
    
def FixedPoint(p0, TOL, N0, g):
    i = 1

    while i <= N0:
        p = g(p0)

        if abs(p - p0) < TOL:
            print(i)
            return p
        else:
            i += 1
            p0 = p

    print('Method failed after ', N0, 'iterations')

fixed_root = FixedPoint(0.9, 10e-8, 300, function_g)

#function works for g1 = x - x^5, not g2 = x + x^5
#function converges fastest with a p0 within 1 of p

print(fixed_root)