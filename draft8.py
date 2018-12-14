# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 13:48:48 2018

@author: James Doering
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import RK45
from scipy.spatial.distance import cdist
import matplotlib.animation as animation


def get_acc(t, y):
    pos = y[:, :2]  # gets positions of all planets
    accel = np.zeros_like(pos)  # creates a dummy array for accelerations
    dists = cdist(pos, pos)  # calculates all distances between objects
    for i, m_i in enumerate(masses):
        for j, m_j in enumerate(masses):
            if i == j:
                continue
            elif m_i == 0:
                continue
            r_ij = pos[j] - pos[i]
            accel[j] += -G * m_i * (r_ij / dists[i,j]**3)
    v = y[:, 2:]
    out = np.hstack((v, accel))
    return out
    


def sim_euler(y0, t0, tf, n, masses, ndim):
    
    nbody = len(masses)
    times = np.linspace(t0, tf, n+1)
    all_dt = np.diff(times)
    
    all_y = np.zeros((n+1, nbody, ndim*2))
    all_y[0] = y0
    
    for i, t in enumerate(times[1:]):
        dt = all_dt[i]
        all_y[i+1] = all_y[i] + dt*get_acc(t, all_y[i])
    
    return times, all_y  # return trajectories


'''CONSTANTS'''
G = 0.667*10**(-10)  # Gravitational Constant G
ms = 0.1984*10**31  # Mass of the Sun
me = 0.5976*10**25  # Mass of the Earth
mj = 0.1903*10**28  # Mass of Jupiter
mb = 0

Re = 0.1495*10**12  # Radius of the Earth's Orbit
Rj = 0.7778*10**12  # Radius of Jupiter's orbit

re = 0.6368*10**7  # Radius of the Earth
rj = 0.6985*10**8  # Radius of Jupiter

T_Jup = 374669228.20514125  # secs it takes for Jupiter to orbit once
T_sat = 126822355.53262155  # secs it takes for satellite to reach Jup
T_Earth = 31572374.99486783  # secs it takes for earth to orbit

earthspeed = 2*np.pi*Re/T_Earth  # earths speed
jupspeed = 2*np.pi*Rj/T_Jup  # jupiters speed

'''STARTING CONDITIONS'''
anglejup = 1.0147877256222952  # angle of jup initial position
anglejupv = 3.6976012547623944  # angle of jup initial velocity

initpos_sun = [0, 0]  # initial position of sun
initpos_earth = [Re, 0]  # initial position of earth
initpos_jup = [Rj*np.cos(anglejup), Rj*np.sin(anglejup)]  # initial pos of jup
initsunv = [0, 0]  # sun initial velocity
initearthv = [0, earthspeed]  # earth initial velocity
initjupv = [jupspeed*np.cos(anglejupv), -jupspeed*np.sin(anglejupv)]  # jup v
'''ADDING A SATELLITE - THE BEBOP'''
initpos_bbop = [Re + re, 0]
initbbopv = initearthv

'''DEFINING Y'''
'''
0 = Sun
1 = Earth
2 = Jupiter
3 = Satellite'''
masses = np.array([ms, me, mj, mb])  # 
y0 = np.zeros([len(masses), 4])
y0[0] = initpos_sun + initsunv
y0[1] = initpos_earth + initearthv
y0[2] = initpos_jup + initjupv
y0[3] = initpos_bbop + initearthv

'''TODO - Set centre of mass velocity zero, stop slow creep of plot'''


times, trajs = sim_euler(y0, 0, 10*36500*86400, 5*36500, masses, 2)
sun_traj = trajs[:, 0, :2]
ear_traj = trajs[:, 1, :2]
jup_traj = trajs[:, 2, :2]
sat_traj = trajs[:, 3, :2]

sunx = sun_traj[:, 0]
suny = sun_traj[:, 1]
earx = ear_traj[:, 0]
eary = ear_traj[:, 1]
jupx = jup_traj[:, 0]
jupy = jup_traj[:, 1]
satx = sat_traj[:, 0]
saty = sat_traj[:, 1]

plt.figure(figsize=[10, 10])
plt.xlim(-10**12, 10**12), plt.ylim(-10**12, 10**12)
plt.plot(sunx, suny, 'ro', label='Sun')
plt.plot(earx, eary, 'b', label='Earth')
plt.plot(jupx, jupy, 'y', label='Jupiter')
plt.plot(satx, saty, 'g', label='Satellite')
plt.show()






