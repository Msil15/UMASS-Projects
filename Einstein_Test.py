import numpy as np
import matplotlib.pyplot as plt
import vpython as vp
vec = vp.vector
from gravipy.tensorial import *
from time import sleep

#define some symbolic variables
t, r, theta, phi, M = symbols('t, r, \\theta, \phi, M')

#create a coordinate 4-vector object instantiating
#the Coordinates class
coord = Coordinates('\chi', [t, r, phi])

#Schwarzchild metric (simplified for 2D motion, G=M=c=1)
Metric = diag(-(1-2/r), 1/(1-2/r), r**2)

#create a metric tensor object instantiating the MetricTensor class
g = MetricTensor('g', coord, Metric)

#calculate Geodesic
tau = Symbol('\\tau')
w = Geodesic('w', g, tau)

def pos_def(posvel, t): #equations of motion for test body
    pos, v = posvel[0], posvel[1]
    acc_r = -M/(pos[0]**2) + pos[0]*(v[1]**2) - 3*M*(v[1]**2)
    acc_phi = -(2*v[0]*v[1])/pos[0]
    acc = [acc_r, acc_phi]
    return np.array([v, acc])

def RK4(diffeq, y0, t, h):
    """ RK4 method for ODEs:
        Given y0 at t, returns y1 at t+h """
    k1 = h*diffeq(y0, t)                    # dy/dt at t
    k2 = h*diffeq(y0+0.5*k1, t + h/2.)      # dy/dt at t+h/2
    k3 = h*diffeq(y0+0.5*k2, t + h/2.)      # dy/dt at t+h/2
    k4 = h*diffeq(y0+k3, t + h)             # dy/dt at t+h
    return y0 + (k1+k4)/6.0 + (k2+k3)/3.0

def pol2cart(r, phi):
    x = r * np.cos(phi)
    y = r * np.sin(phi)
    return [x, y]

def set_scene(r, phi, Sch_rad, t):
    cart_pos = pol2cart(r, phi)
    vp.canvas(title= 'Schwarzschild Black Hole', background = vec(1, 1, 1))
    body = vp.sphere(pos = vec(*cart_pos, 0), color = vec(0,0,1), radius=0.01, make_trail=1)
    Blk_hle = vp.sphere(pos = vec(0,0,0), color = vec(0,0,0), radius = Sch_rad, emissive=True)
    ISCO = vp.ring(pos = vec(0,0,0), color = vec(1, 1, 0), thickness = 0.05, axis = vec(0,0,1), radius = r_isco, opacity=0.5)
    marginally_bound = vp.ring(pos = vec(0,0,0), color = vec(0, 1, 1), thickness = 0.05, axis = vec(0,0,1), radius = r_mb, opacity=0.5)
    photon_sphere = vp.ring(pos = vec(0,0,0), color = vec(0, 1, 0), thickness = 0.05, axis = vec(0,0,1), radius = r_ph, opacity=0.5)
    global coordinate_time
    coordinate_time = vp.label(pos = vec(4*rs,4,0), text = '0')
    global proper_time
    proper_time = vp.label(pos = vec(4*rs,3,0), text = '0')
    return body

def testbody_motion(posvel):
    k, i, h = 0.0, 0.0, 0.001 #initial time & timestep
    testbody = set_scene(*posvel[0], rs, i)
    upper_x_bound = np.sqrt(rs)
    lower_x_bound = -np.sqrt(rs)
    upper_y_bound = np.sqrt(rs)
    lower_y_bound = -np.sqrt(rs)
    sleep(0.5)
    while i < 160000:
        vp.rate(5000)
        posvel = RK4(pos_def, posvel, i, h)
        cart_pos = pol2cart(posvel[0,0], posvel[0,1])
        testbody.pos = vec(*cart_pos, 0)
        coordinate_time.text = 'Coordinate Time ' + str(round(i, 1))
        proper_time.text = 'Proper Time ' + str(round(k, 1))
        #print(testbody.pos)
        i += 1
        k += 1*np.sqrt(1-rs/posvel[0,0])
        if lower_x_bound <= testbody.pos.x <= upper_x_bound and lower_y_bound <= testbody.pos.y <= upper_y_bound:
            print('EVENT HORIZON')
            break
 

G = 1
M = 1 #solar masses
c = 1
rs = 2*M
r_isco = 3*rs #innermost stable circular orbit
r_mb = 2*rs #marginally bound orbits
r_ph = 1.5*rs #photon sphere
r, phi, v = 3*rs-0.2, np.pi/2, [0.0, 0.1] #initial pos, vel

#ISCO CASE; 
#r, phi, v = 3*rs, np.pi/2, [0.0, 0.1]
#INFALL CASE; 
#r, phi, v = 3*rs-0.2, np.pi/2, [0.0, 0.1]
#5M SPHERE CASE;
#r, phi, v = 3*rs-0.2, np.pi/2, [0.0, 0.02] (remove 0 from thickness of each object)
#PHOTON CASE
#r, phi, v = 15*rs, 2/15, [-1, 0] #initial pos, vel

testbody_motion(np.array([[r, phi], v]))

