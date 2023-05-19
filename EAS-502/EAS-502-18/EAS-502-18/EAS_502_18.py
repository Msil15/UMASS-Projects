import numpy as np
from scipy import special, integrate

def omega_x(x):
    return np.exp(-x)

def top_B(x, omega, L):
    return x*omega(x)*(L(x))**2

def bottom_B(x, omega, L):
    return omega(x)*(L(x))**2

def top_C(x, omega, L1, L2):
    return x*omega(x)*L1(x)*L2(x)

def bottom_C(x, omega, L2):
    return omega(x)*(L2(x))**2

def L0(x):
    return x**0

B1 = integrate.quad(top_B, 0, np.inf, args=(omega_x, L0))[0]/integrate.quad(bottom_B, 0, np.inf, args=(omega_x, L0))[0]

def L1(x):
    return x - B1

B2 = integrate.quad(top_B, 0, np.inf, args=(omega_x, L1))[0]/integrate.quad(bottom_B, 0, np.inf, args=(omega_x, L1))[0]
C2 = integrate.quad(top_C, 0, np.inf, args=(omega_x, L1, L0))[0]/integrate.quad(bottom_C, 0, np.inf, args=(omega_x, L0))[0]

def L2(x):
    return (x - B2)*L1(x) - C2

B3 = integrate.quad(top_B, 0, np.inf, args=(omega_x, L2))[0]/integrate.quad(bottom_B, 0, np.inf, args=(omega_x, L2))[0]
C3 = integrate.quad(top_C, 0, np.inf, args=(omega_x, L2, L1))[0]/integrate.quad(bottom_C, 0, np.inf, args=(omega_x, L1))[0]

def L3(x):
    return (x - B3)*L2(x) - C3*L1(x)

def d_alpha(x, omega, L):
    return omega(x)*L(x)

def d_coef(x, omega, L, func):
    return omega(x)*func(x)*L(x)

def function(x):
    return x**2

alph_0 = integrate.quad(d_alpha, 0, np.inf, args=(omega_x, L0))[0]
alph_1 = integrate.quad(d_alpha, 0, np.inf, args=(omega_x, L1))[0]
alph_2 = integrate.quad(d_alpha, 0, np.inf, args=(omega_x, L2))[0]
alph_3 = integrate.quad(d_alpha, 0, np.inf, args=(omega_x, L3))[0]

coef_0 = integrate.quad(d_coef, 0, np.inf, args=(omega_x, L0, function))[0]
coef_1 = integrate.quad(d_coef, 0, np.inf, args=(omega_x, L1, function))[0]
coef_2 = integrate.quad(d_coef, 0, np.inf, args=(omega_x, L2, function))[0]
coef_3 = integrate.quad(d_coef, 0, np.inf, args=(omega_x, L3, function))[0]


