import numpy as np

def func(x, t):
    k = 6.22e-19
    n1 = 2e3
    n2 = 2e3
    n3 = 3e3
    dxdt = k*((n1 - (x/2))**2)*((n2 - (x/2))**2)*((n3 - (3*x/4))**3)
    return dxdt

#Method based on 
#'Computational modeling and visualization with Python', Jay Wang
def RK4(func, omega, t, h):
    k1 = h*func(omega, t)                    
    k2 = h*func(omega+k1/2, t + h/2)      
    k3 = h*func(omega+k2/2, t + h/2)      
    k4 = h*func(omega+k3, t + h)             
    return omega + (1/6)*(k1 + 2*k2 + 2*k3 + k4)


omega = RK4(func, 0, 0, 0.2)
print(omega)

