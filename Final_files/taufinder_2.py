#! /usr/bin/env python
'''
Uses secant method to find the half-time to steadystate


'''

import numpy as np
from scipy.integrate import odeint
from mpmath import findroot

def find_tau(odefun, x0, which_component=0, tol = 10**-5, findroot_solver = 'illinois', findroot_t0 = [0,10000], ss = False):
    '''
    odefun has to be of the form odefun(t, Y) = dY/dt, just like the input for scipy.integrate.ode,
    even though we're actually using odeint
    '''
    f = lambda Y, t: odefun(t, Y)
    F = lambda t: odeint(f, x0, [0, t])[-1]

    Y_ss = F(1e+10) # 10^10 seconds is a few thousand years. Should be enough to get to the steady state.
    half_ss = Y_ss[which_component]/2

    def diff_from_half_ss(t):
        return F(t)[which_component]-half_ss
    if ss == True:
        return half_ss * 2

    return findroot(f=diff_from_half_ss, x0 = [0, 10**4], tol = tol, solver = findroot_solver), Y_ss
