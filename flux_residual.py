import numpy as np
import os
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy import stats
import copy

import sys
sys.path.append("../")
from FIESTA_functions import *
#--------------------------------------------------------------------
# functions 
#--------------------------------------------------------------------

def FT(signal, spacing):
	n 			= signal.size
	fourier 	= np.fft.rfft(signal, n)
	freq 		= np.fft.rfftfreq(n, d=spacing)
	A 			= np.abs(fourier)
	phi 		= np.angle(fourier)
	return [A, phi, freq]

#--------------------------------------------------------------------
# read data
#--------------------------------------------------------------------
valid_ccf = sorted(glob.glob('../../AstroData/FIESTA-HARPS-N/valid_ccfs/*.fits'))
N_file 	= len(valid_ccf)
V_grid 	= (np.arange(49) - 24) * 0.82
idx 	= (-16<V_grid) & (V_grid<16)
V_grid 	= V_grid[idx]
spacing = np.diff(V_grid)[0]

hdulist     = fits.open(valid_ccf[1234])
data        = hdulist[1].data	
CCF 		= 1 - data[69,:] / np.mean(data[69,~idx])
CCF 		= CCF[idx]

#--------------------------------------------------------------------
# DFT 
#--------------------------------------------------------------------
A, phi, freq= FT(CCF, spacing)
res_rms 	= np.zeros(freq.size)
res 		= np.zeros((freq.size, CCF.size))
ft 			= np.fft.rfft(CCF, CCF.size)
for i in range(ft.size):
	pseudo_ft = copy.copy(ft)
	if i < (ft.size-1):
		pseudo_ft[(i+1):] = 0
	pseudo_ift 	= np.fft.irfft(pseudo_ft, len(CCF))
	res[i,:] 	= pseudo_ift - CCF 
	res_rms[i] 	= np.std(res[i,:])


#--------------------------------------------------------------------
# plots 
#--------------------------------------------------------------------

plt.rcParams.update({'font.size': 12})
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']


for i in range(6):
	plt.plot(V_grid, res[i,:], '-', label=str(i))
	plt.xlabel(r'$v$ [km/s]')
	plt.ylabel('CCF residual')	
	plt.legend(loc=4)
plt.savefig('Flux_residual.png')
plt.show()


k_max = 6

plt.rcParams.update({'font.size': 14})
plt.subplots(figsize=(10, k_max))
plt.gcf().subplots_adjust(wspace=0.35)

for i in range(k_max):

    ax = plt.subplot(k_max, 2, 2*i + 1)
    if i == 0:
    	plt.plot(V_grid, (np.mean(CCF))*np.ones(len(V_grid)), 'k-')	
    else:
    	plt.plot(V_grid, res[i,:] - res[i-1,:], 'k-')
    plt.ylabel(r'$k$=%d' %(i))
    if i == 0:
        plt.title('Basis functions')

    if i!=k_max-1:
        ax.set_xticks([])
    else:
        plt.xlabel('V grid [km/s]')

    ax = plt.subplot(k_max, 2, 2*i+2)
    plt.plot(V_grid, -res[i,:], 'k-')
    if i == 0:
        plt.title('Residual')

    if i!=k_max-1:
        ax.set_xticks([])
    else:
        plt.xlabel('V grid [km/s]')

plt.savefig('FIESTA_demo.pdf')
plt.show()



# Plot with two axes 
# https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.secondary_yaxis.html#matplotlib.axes.Axes.secondary_yaxis
# https://matplotlib.org/stable/gallery/subplots_axes_and_figures/secondary_axis.html

plt.rcParams.update({'font.size': 14})
fig, ax = plt.subplots(figsize=(10, k_max))

def SNR(x):
	return np.median(1-CCF)/x

inverse = SNR

for i in range(freq.size-1):
	ax.plot(i, res_rms[i], 'ko', alpha = 1)
	ax.set_xticks(range(freq.size-1))
	ax.set_xlabel('Number of Fourier modes used')
	ax.set_yscale('log')
	ax.set_ylabel('Modelling error')
	ax.grid(b=True, which='both')
	secax = ax.secondary_yaxis('right', functions=(SNR, inverse))	
	secax.set_ylabel('SNR')
plt.savefig('Flux_residual_rms.pdf')
plt.show()




#--------------------------------------------------------------------
# archived  
#--------------------------------------------------------------------
for i in range(freq.size-1):
	plt.plot(i, res_rms[i], 'ko', alpha = 1)
	# plt.plot(i, res_rms[i], '-', alpha = 0.5)
	plt.text(i+0.2, res_rms[i]*1.1, '{:.0f}'.format(np.median(1-CCF)/res_rms[i]), fontsize=8)
	plt.xticks(range(freq.size-1))
	plt.xlabel('Number of Fourier modes used')
	plt.yscale('log')
	plt.ylabel('Residual rms')
plt.savefig('Flux_residual_rms.png')
plt.show()

#--------------------------------------------------------------------
# testing 
#--------------------------------------------------------------------

fig, ax1 	= plt.subplots()
ax2 		= ax1.twinx()
for i in range(freq.size-1):
	ax1.plot(i, res_rms[i], 'ko-', alpha = 0.5)
	# plt.plot(i, 1/res_rms[i], 'ko-', alpha = 0)
ax1.set_xlabel('mode')
ax1.set_ylabel('Residual rms')
ax2.set_ylim((res_rms[freq.size-2], res_rms[0]))
ax2.set_ylabel('SNR')
ax2.set_ylim((1/res_rms[freq.size-2], 1/res_rms[0]))
plt.show()




x = np.arange(0, 10, 0.1)
y1 = 0.05 * x**2
y2 = -1 *y1

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(x, y1, 'g-')
ax2.plot(x, y2, 'b-')

ax1.set_xlabel('X data')
ax1.set_ylabel('Y1 data', color='g')
ax2.set_ylabel('Y2 data', color='b')

plt.show()

