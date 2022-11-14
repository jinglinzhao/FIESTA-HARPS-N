import numpy as np
from scipy.optimize import minimize
from datetime import timedelta, date
from datetime import datetime
from scipy import stats
import random
from sklearn.model_selection import KFold

import george
from george import kernels
import celerite
from celerite import terms


def neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)

def GP(t, y, yerr):
    w0 = 2*np.pi/5.2
    Q = 10
    S0 = np.var(y) / (w0 * Q)
    bounds = dict(log_S0=(-15, 15), log_Q=(-15, 15), log_omega0=(-15, 15))
    kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0), bounds=bounds)
    # kernel.freeze_parameter("log_omega0")

    gp = celerite.GP(kernel, mean=np.mean(y))
    gp.compute(t, yerr) 

    initial_params = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()

    r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(y, gp))
    gp.set_parameter_vector(r.x)
    
    return gp


def GP_fit_Matern52Kernel(x, yerr, r):
    kernel 	= kernels.Matern52Kernel(r**2)
    gp 		= george.GP(kernel)
    gp.compute(x, yerr)
    return gp


def GP_fit_p1(t, y, yerr, p):
    w0 = 2*np.pi/p
    Q = 6
    S0 = np.var(y) / (w0 * Q)
    bounds = dict(log_S0=(-15, 15), log_Q=(-15, 15), log_omega0=(-15, 15))
    kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0), bounds=bounds)
    kernel.freeze_parameter("log_omega0")
    kernel.freeze_parameter("log_Q")

    gp = celerite.GP(kernel, mean=np.mean(y))
    gp.compute(t, yerr)  # You always need to call compute once.

    initial_params = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()

    r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(y, gp))
    gp.set_parameter_vector(r.x)
    
    return gp


def GP_fit_Q1(t, y, yerr, Q):
    S0 = np.var(y) / (w0 * Q)
    bounds = dict(log_S0=(-15, 15), log_Q=(-15, 15), log_omega0=(np.log(2*np.pi/7), np.log(2*np.pi/3)))
    kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0), bounds=bounds)
    kernel.freeze_parameter("log_Q")

    gp = celerite.GP(kernel, mean=np.mean(y))
    gp.compute(t, yerr)  # You always need to call compute once.

    initial_params = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()

    r = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(y, gp))
    gp.set_parameter_vector(r.x)
    
    return gp


