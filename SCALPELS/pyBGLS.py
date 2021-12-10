# -*- coding: utf-8 -*-
#
#   pyBGLS - Tools for computing Bayesian Generalised Lomb-Scargle periodograms.
#
#   Copyright (C) 2021  Prof Andrew Collier Cameron, University of St Andrews

#   Uses algorithms published by Mortier et al, A&A 573, A101 (2015)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

"""
    import sys
    sys.path.insert(0, '/Users/acc4/OneDrive - University of St Andrews/CodeDevelopment/PythonTools/pyBGLS/')
    from pyBGLS import bglsconst, bglsgrid, bglsfreq
    
    """

import numpy as np

# bglsconst(t,d,e) performs summations on a set of data with errors e
# observed at times t, returning (W,Y,YYh) which are needed by bglsfreq.
def bglsconst(t, d, e):
    #  local vars w, u, W, Y, dd, YYh
    w = 1/e/e
    u = t*0 + 1
    W = np.dot(w,u)
    Y = np.dot(w,d)
    dd = d * d
    YYh = np.dot(w,dd)
    return W, Y, YYh

# bglsgrid(t, c) is an optional function to set up frequency grid.
def bglsgrid(t, c):
    # local vars n, nc, tlen, phigh, df, tinc, plow, fmax, nf, f, j
    n = len(t)
    nc = c[n-1] - c[0]
    tlen = t[n-1] - t[0]
    phigh = 2 * tlen
    df = 1/phigh
    tinc = tlen/nc
    plow = 2 * tinc
    fmax = 1/plow
    nf = fmax/df
    j = np.linspace(1,nf)
    f = (j+1)*df
    return f

# bglsfreq(t,d,e,f,W,Y) returns a vector of log relative likelihoods
# corresponding to the input grid of frequencies f, for a set of data
# with errors e observed at times t. The sums W and Y should have been
# precalculated using the bglsconst function.
def bglsfreq(t, d, e, f, W, Y):
    # local vars w, y, x, theta, c, s, YCh, YSh, C, S, CCh, SSh, K, L, M
    w = 1/e/e 
    y = np.dot(w,np.sin(4*np.pi*f*t))
    x = np.dot(w,np.cos(4*np.pi*f*t))
    theta = np.arctan2(y,x)/2
    c = np.cos(2*np.pi*f*t - theta)
    s = np.sin(2*np.pi*f*t - theta)
    YCh = np.dot(w,(d*c))
    YSh = np.dot(w,(d*s))
    C = np.dot(w,c)
    S = np.dot(w,s)
    CCh = np.dot(w,(c*c))
    SSh = np.dot(w,(s*s))
    #CSh = np.dot(w,(c*s)) # should be zero
    #print('DBG> CSh =',CSh)
    K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2*CCh*SSh)
    L = (Y*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
    M = (YCh*YCh*SSh + YSh*YSh*CCh)/(2*CCh*SSh)
    logP = M - (L*L)/(4*K) - np.log(np.abs(K)*CCh*SSh)/2
    return logP

# ggls(t,d,e,f,W,Y) returns a vector of logl, power and sinusoid amplitudes
# corresponding to the input grid of frequencies f, for a set of data
# with errors e observed at times t. The sums W and Y should have been
# precalculated using the bglsconst function.
def gls(t, d, e, f, W, Y):
    # local vars w, y, x, theta, c, s, YCh, YSh, C, S, CCh, SSh, K, L, M
    w = 1/e/e * 1/W
    c = np.cos(2*np.pi*f*t)
    s = np.sin(2*np.pi*f*t)
    C = np.dot(w,c)
    S = np.dot(w,s)
    CSh = np.dot(w,(c*s)) # ZK2009
    CS = CSh - C*S # ZK2009
    CCh = np.dot(w,(c*c))
    SSh = np.dot(w,(s*s))
    CC = CCh - C*C
    SS = SSh - S*S
    theta = np.arctan2(2*CS,CC-SS)/2
    # Use orthogonalised quantities from now on
    c = np.cos(2*np.pi*f*t - theta)
    s = np.sin(2*np.pi*f*t - theta)
    C = np.dot(w,c)
    S = np.dot(w,s)
    CSh = np.dot(w,(c*s)) # ZK2009
    CS = CSh - C*S # ZK2009
    CCh = np.dot(w,(c*c))
    SSh = np.dot(w,(s*s))
    CC = CCh - C*C
    SS = SSh - S*S
    YCh = np.dot(w,(d*c))
    YSh = np.dot(w,(d*s))
    YC = YCh - Y*C
    YS = YSh - Y*S
    YYh = np.dot(w,d*d)
    YY = YYh - Y*Y
    D = CC*SS - CS*CS
    power = (SS*YC*YC + CC*YS*YS - 2*CS*YC*YS)/(YY*D)
    Kc = (YC*SS-YS*CS)/D
    VarKc = SS/D/W # Don't forget that the inverse-variance weights were scaled by 1/W!
    VarKs = CC/D/W # Don't forget that the inverse-variance weights were scaled by 1/W!
    Ks = (YS*CC-YC*CS)/D
    Kamp = np.sqrt(Kc*Kc+Ks*Ks)
    Kerr = np.sqrt(Kc*Kc*VarKc + Ks*Ks*VarKs)/Kamp
    gamma = Y - Kc*C - Ks*S
    K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2*CCh*SSh)
    L = (Y*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
    M = (YCh*YCh*SSh + YSh*YSh*CCh)/(2*CCh*SSh)
    logP = M - (L*L)/(4*K) - np.log(np.abs(K)*CCh*SSh)/2
    
    return Kamp,power,Kerr



