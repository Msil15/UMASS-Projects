# NOTE: Numerov's method caused the eigenvalues to vary by about a tenth.
# It also improved the time it took to find the eigenvalues.
#

import numpy as np
import matplotlib.pyplot as plt

def RK45n(diffeq, y0, t, h):      # RK45 method, based on NR, Press
    a2, a3, a4, a5, a6 = 0.2, 0.3, 0.6, 1.0, 0.875
    b21, b31, b32, b41, b42, b43 = 0.2, 3./40., 9./40., 0.3, -0.9, 1.2
    b51, b52, b53, b54 = -11./54., 2.5,  -70./27., 35./27.
    b61, b62, b63, b64, b65 = [1631./55296., 175./512., 575./13824.,
                               44275./110592., 253./4096.]
    c1, c3, c4, c6 = 37./378., 250./621., 125./594., 512./1771.
    
    n, y1 = len(y0), [0.0]*len(y0)
    k1 = diffeq(y0, t)
    for i in range(n):
        y1[i] = y0[i] + h*b21*k1[i]
    k2 = diffeq(y1, t + a2*h)
    for i in range(n):             
        y1[i] = y0[i] + h*(b31*k1[i] + b32*k2[i])
    k3 = diffeq(y1, t + a3*h)
    for i in range(n):             
        y1[i] = y0[i] + h*(b41*k1[i] + b42*k2[i] + b43*k3[i])
    k4 = diffeq(y1, t + a4*h)
    for i in range(n):             
        y1[i] = y0[i] + h*(b51*k1[i] + b52*k2[i] + b53*k3[i] + b54*k4[i])
    k5 = diffeq(y1, t + a5*h)
    for i in range(n):             
        y1[i] = y0[i] + h*(b61*k1[i] + b62*k2[i] + b63*k3[i] + b64*k4[i]
                         + b65*k5[i])
    k6 = diffeq(y1, t + a6*h)
    for i in range(n): 
        y1[i] = y0[i] + h*(c1*k1[i] + c3*k3[i] + c4*k4[i] + c6*k6[i])
    return y1

def bisect(f, a, b, eps=1.e-6):  # user-defined f(x) and bracket [a,b]
    fa, fb, gap = f(a), f(b), abs(b-a)  # end points and initial gap
    if (fa*fb > 0.0):                   # no root in bracket
        print('Bisection error: no root bracketed')
        return None       
    elif fa == 0.0:   return a
    elif fb == 0.0:   return b
    
    while (True):
        xmid = 0.5*(a+b)
        fmid = f(xmid)
        if (fa*fmid > 0.0):         # root in [xmid, b]
            a, fa = xmid, fmid      # set a=xmid and save a function call
        else: b=xmid                # root in [a, xmid]
        if (fmid == 0.0 or abs(b-a) < eps*gap): break   # root found @ \label{line:biseceps} @

    return xmid

def V(x):                   # potential
    return -V0 if (-1.5 < x < - 0.9 or -0.7 < x < -0.1 or 0.1 < x < 0.7 or 0.9 < x < 1.5) else 0
    
def sch(x):            # Schrodinger eqn
    m = 1 #amu
    h_bar = 1 #natural units
    return (2*m/(h_bar**2))*(E-V(x))
    
def numerov(f, u, n, x, h):     # Numerov integrator for $u''+f(x)u=0$
    nodes, c = 0, h*h/12.       # given $[u_0,u_1]$, return $[u_0,u_1,...,u_{n+1}]$
    f0, f1 = 0., f(x+h)
    for i in range(n):
        x += h
        f2 = f(x+h)             # Numerov method below, 
        u.append((2*(1-5*c*f1)*u[i+1] - (1+c*f0)*u[i])/(1+c*f2))  
        f0, f1 = f1, f2
        if (u[-1]*u[-2] < 0.0): nodes += 1
    return u, nodes             # return u, nodes
    
def shoot(En):
    global E                    # E needed in f(r)
    E, c, xm = En, (h*h)/6., xL + M*h
    wfup, nup = numerov(sch, [0,.1], N, xL, h)
    wfdn, ndn = numerov(sch, [0,.1], N, xR, -h)     # $f'$ from 
    dup = ((1+c*sch(xm+h))*wfup[-1] - (1+c*sch(xm-h))*wfup[-3])/(h+h)
    ddn = ((1+c*sch(xm+h))*wfdn[-3] - (1+c*sch(xm-h))*wfdn[-1])/(h+h)
    return dup*wfdn[-2] - wfup[-2]*ddn

a, b, V0 = 4.0, 1.0, 6.             # double well widths, depth
xL, xR, N = -4*a, 4*a, 500          # limits, intervals
xa = np.linspace(xL, xR, 2*N+1)     # grid
h, M = xa[1]-xa[0], 2               # step size, M=matching point
E1, dE, j = -V0, 0.01, 1

plt.figure()
while (E1 < 0):     # find E, calc and plot wave function
    if (shoot(E1) * shoot(E1 + dE) < 0):        # bracket E
        E = bisect(shoot, E1, E1 + dE, 1.e-8)
        print ('Energy found: %.3f' %(E))
        wfup, psiup = numerov(sch, [0., .1], N+M-1, xL, h)      # compute WF
        wfdn, psidn = numerov(sch, [0., .1], N-M-1, xR, -h)
        psix = np.concatenate((wfup[:-1], wfdn[::-1]))  # combine WF
        psix[M:] *= wfup[-1]/wfdn[-1]                 # match WF
        ax = plt.subplot(2,2,j)
        ax.plot(xa, psix/max(psix))                     # plot WF
        ax.plot(xa, np.vectorize(V)(xa)/(2*V0))         # overlay V
        ax.set_xlim(-a,a)
        ax.text(2.2,0.7, '%.3f' %(E))
        if (j == 1 or j == 3): ax.set_ylabel(r'$\psi$')
        if (j == 3 or j == 4): ax.set_xlabel('$x$')
        if (j<4): j += 1            # 4 plots max
    E1 += dE
plt.show()