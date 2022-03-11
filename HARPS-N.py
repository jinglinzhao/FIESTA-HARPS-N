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
from HARPS_N_functions import *
# import FIESTA_functions from the directory above 
# https://gist.github.com/MRobertEvers/55a989b4883ea8d7715d2e156627f034


# Read the CCFs

import pandas as pd
import shutil

df = pd.read_csv('./harpn_sun_release_timeseries_2015-2018.csv')

if 0: 
	file_ccf = sorted(glob.glob('./ccfs/*.fits'))
# min(df['obs_quality']) = 0.99 
# => all the files in the csv table is eligible for "no clouds"
# min(df['drs_quality']) = True
# => all files have good quality of the data reduction software

# simply check if files are matching 
	count = 0
	for n in range(len(df)):
		print(n)
		for j in n+np.arange(80)-31:
			if df['filename'][n][:-5] == file_ccf[j][7:38]:
				shutil.move(file_ccf[j], "./valid_ccfs")
				count += 1

	for n in np.arange(34534, len(df)):
		print(n)
		for j in range(len(file_ccf)):
			if df['filename'][n][:-5] == file_ccf[j][7:38]:
				shutil.move(file_ccf[j], "./valid_ccfs")
				count += 1

# All the matching files are moved to the valid_ccfs folder. 
valid_ccf = sorted(glob.glob('../../AstroData/FIESTA-HARPS-N/valid_ccfs/*.fits'))

if 0:
	for i in range(len(df)):
		if (valid_ccf[i][13:44] == df['filename'][i][:-5]) == False: 
			print(i)

N_file 	= len(valid_ccf)
V_grid 	= (np.arange(49) - 24) * 0.82
idx 	= (-10<V_grid) & (V_grid<10)
# V_grid 	= V_grid[idx]
CCF 	= np.zeros((len(V_grid), N_file))
eCCF 	= np.zeros((len(V_grid), N_file))

# Normalised CCF for the first 100 files 
if 0: 
	for n in range(100):
	 	hdulist     = fits.open(valid_ccf[n])
	 	data         = hdulist[1].data
	 	plt.plot(V_grid[idx], data[0,idx] / np.mean(data[0,~idx]))
	plt.show()	


if 0: 
	#--------------------------------------------------------------------
	# Part 1: individual analysis 
	#--------------------------------------------------------------------
	# Go to part 2 for daily analysis

	# FWHM = 2.355 sigma for a Gaussian function?

	fwhm_raw 	= np.array(df['fwhm_raw'])
	sigma 		= np.zeros(N_file)
	V_centre 	= np.zeros(N_file)
	CCF 		= np.zeros((len(V_grid), N_file))

	# Read and normalise the CCFs
	for n in range(len(valid_ccf)):
	# for n in range(100):	
		print(n)
		hdulist     = fits.open(valid_ccf[n])
		data        = hdulist[1].data	
		CCF[:,n] 	= 1 - data[69,:] / np.mean(data[69,~idx])
		data2       = hdulist[2].data
		eCCF[:,n] 	= data2[69,:] / np.mean(data[69,~idx])
		popt, pcov 	= curve_fit(gaussian, V_grid, CCF[:,n], p0=[0.5, (max(V_grid)+min(V_grid))/2, 3, 0])
		sigma[n] 	= popt[2]
		V_centre[n] = popt[1]


	# plot the sigma of all orders // orders range from 0 to 68
	sigma_by_order = np.zeros((69, len(valid_ccf)))
	for n in range(len(valid_ccf)):
	# for n in range(100):	
		print(n)
		hdulist     = fits.open(valid_ccf[n])
		data        = hdulist[1].data	
		for order in range(69):
			if not data[order,:].any():
				sigma_by_order[order, n] = 0
			else:		
				ccf 	= 1 - data[order,:] / np.mean(data[order,~idx])
				popt, pcov 	= curve_fit(gaussian, V_grid, ccf, p0=[0.5, (max(V_grid)+min(V_grid))/2, 1, 0])
				sigma_by_order[order, n] = popt[2]


	for order in range(69):
		fig, axes = plt.subplots(figsize=(12, 4))
		plt.plot(bjd, sigma_by_order[order, :], '.', alpha=0.1)
		plt.xlabel('date_bjd')
		plt.ylabel('sigma')
		plt.title('order' + str(order))
		plt.savefig('sigma_by_order' + str(order) + '.png')
		plt.close()



	# =============================================================================
	# for order in range(69):
	# 	fig, axes = plt.subplots(figsize=(12, 4))
	# 	plt.plot(V_grid, data[order, :], '-')
	# 	plt.xlabel('V_grid')
	# 	plt.ylabel('last_CCF_by_order')
	# 	plt.title('order' + str(order))
	# 	plt.savefig('last_CCF_by_order' + str(order) + '.png')
	# 	plt.close()	
	# =============================================================================

	# max(bjd) - min(bjd) = 1093.819484889973
	sigma_by_order_itp = np.zeros((69, 10000))
	bjd_itp = np.linspace(min(bjd), max(bjd), num=10000)
	for order in range(69):
		f = interp1d(bjd, sigma_by_order[order, :])
		sigma_by_order_itp[order, :] = f(bjd_itp)



	X, Y = np.mgrid[0:9999:complex(0, 10000), 0:68:complex(0, 69)]
	Z = np.transpose(sigma_by_order_itp)
	for i in range(69):
		Z[:, i] = Z[:, i] - np.mean(Z[:, i])

	fig, ax0, = plt.subplots()

	# c = ax0.pcolor(X, Y, Z,
	#                norm=LogNorm(vmin=Z.min(), vmax=Z.max()), cmap='PuBu_r')
	# fig.colorbar(c, ax=ax0)

	# viridis = plt.cm.get_cmap('viridis', 12)
	c = ax0.pcolor(X, Y, Z, cmap='PuBu_r', alpha=.9)
	fig.colorbar(c, ax=ax0)
	plt.savefig('sigma_by_order_2D.png')	
	plt.close()


	clean_orders = [2, 4, 6, 8, 9, 11, 18, 19, 20, 23, 24, 27, 28, 29, 30, 34,37, 42, 43, 44, 45, 46,47,48,49,50,51,52,53, 55,56,57, 59,61,62, 64,65,66]

	CCF_clean 		= np.zeros((len(V_grid), N_file))
	sigma_clean 	= np.zeros(N_file)
	# Read and normalise the CCFs
	for n in range(len(valid_ccf)):
	# for n in range(100):	
		print(n)
		hdulist     = fits.open(valid_ccf[n])
		data        = hdulist[1].data
		ccf_temp 	= sum(data[clean_orders,:])
		CCF_clean[:,n] 	= 1 - ccf_temp / np.mean(ccf_temp[~idx])
		# eCCF[:,n] = (data2[n+1,idx]**0.5 / max(data[n+1,idx])) 
		popt, pcov 	= curve_fit(gaussian, V_grid, CCF_clean[:,n], p0=[0.5, (max(V_grid)+min(V_grid))/2, 1, 0])
		sigma_clean[n] 	= popt[2]



	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, sigma, '.', alpha=0.2, label='all_orders')
	plt.plot(bjd, sigma_clean, '.', alpha=0.1, label='clean_orders')
	plt.xlabel('date_bjd')
	plt.ylabel('sigma')
	plt.legend()
	plt.savefig('sigma_clean.png')
	# plt.show()
	plt.close()

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, sigma - np.mean(sigma), '.', alpha=0.1, label='all_orders')
	plt.plot(bjd, sigma_clean - np.mean(sigma_clean), '.', alpha=0.05, label='clean_orders')
	plt.xlabel('date_bjd')
	plt.ylabel('sigma')
	plt.legend()
	plt.savefig('sigma_clean2.png')
	# plt.show()
	plt.close()

	def moving_average(x, w):
	    return np.convolve(x, np.ones(w), 'valid') / w

	fig, axes = plt.subplots(figsize=(12, 4))
	diff = sigma - np.mean(sigma) - (sigma_clean - np.mean(sigma_clean))
	plt.plot(bjd, diff, '.', alpha=0.1)
	plt.plot(bjd[45:34506], moving_average(diff, 90), '-', alpha=1)
	plt.xlabel('date_bjd')
	plt.ylabel('sigma_diff')
	plt.savefig('sigma_clean3.png')
	# plt.show()
	plt.close()


	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, df['rv_raw'], '.', alpha=0.2)
	plt.xlabel('date_bjd')
	plt.ylabel('rv_raw')
	plt.savefig('rv_raw.png')
	# plt.show()
	plt.close()

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, df['rv'], '.', alpha=0.2)
	plt.xlabel('date_bjd')
	plt.ylabel('rv (heliocentric)')
	plt.savefig('rv.png')
	# plt.show()
	plt.close()


	fig, axes = plt.subplots(figsize=(12, 4))
	# plt.plot(bjd, df['rv'], '.', alpha=0.2, label="raw, i.e. df.['rv']")
	# plt.plot(bjd_daily[rv_daily!=0], rv_daily[rv_daily!=0], 'x', alpha=0.8, label='daily averages')
	# plt.plot(bjd, df['rv_raw'] - df['berv_bary_to_helio'] - df['rv_diff_extinction'] -df['rv'], '.', alpha=0.1, label="__")
	plt.plot(bjd, df['rv_diff_extinction'], '-', alpha=0.1, label="__")
	plt.legend()
	plt.xlabel('date_bjd')
	plt.ylabel('rv (heliocentric) [m/s]')
	plt.savefig('rv2.png')
	# plt.show()
	# plt.close()


	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, fwhm_raw, '.', alpha=0.2)
	plt.xlabel('bjd')
	plt.ylabel('fwhm_raw')
	plt.savefig('fwhm_raw.png')
	# plt.show()
	plt.close()

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, df['fwhm'], '.', alpha=0.2)
	plt.xlabel('bjd')
	plt.ylabel('fwhm')
	plt.savefig('fwhm.png')
	# plt.show()
	plt.close()

	# compare the fwhm and sigma
	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, sigma, '.', alpha=0.2)
	plt.xlabel('bjd')
	plt.ylabel('sigma')
	plt.savefig('sigma_correct.png')
	# plt.show()
	plt.close()

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, fwhm_raw/sigma/1000, '.', alpha=0.2)
	plt.xlabel('bjd')
	plt.ylabel('ratio')
	plt.savefig('ratio.png')
	# plt.show()
	plt.close()


#----------------------------------
# Part 2: Daily average analysis 
#----------------------------------

# first, generate a df file only with abs(df['rv_diff_extinction']) < 0.1
array_index = np.arange(N_file)
df 			= df.drop(array_index[abs(df['rv_diff_extinction']) >= 0.1])

rv_raw 		= np.array(df['rv_raw'])
rv 			= np.array(df['rv'])
rv_err 		= np.array(df['rv_err'])
bjd 		= np.array(df['date_bjd'])
fwhm_raw 	= np.array(df['fwhm_raw'])

# if everything goes alright
# then correct the measured CCFs

bjd_daily 			= np.zeros(632)
CCF_daily 			= np.zeros((len(V_grid), 632))
eCCF_daily			= np.zeros((len(V_grid), 632))
rv_daily			= np.zeros(632)
rv_raw_daily		= np.zeros(632)
erv_daily			= np.zeros(632)
widening_var 		= (df['fwhm']**2 - fwhm_raw**2) / 1e6
widening_var_daily 	= np.zeros(632)

bis_daily 			= np.zeros(632)
fwhm_daily 			= np.zeros(632)
ebis_daily			= np.zeros(632)
efwhm_daily 		= np.zeros(632)


date0 		= df.iloc[0][0][8:18] # first elegible file 
k 			= 0
i0			= 0 
CCF_daily_sum = 0
eCCF_daily_sum = 0

for i in range(len(df)):
	hdulist     = fits.open('../../AstroData/FIESTA-HARPS-N/valid_ccfs/' + df.iloc[i][0][:-5] + '_CCF_A.fits')
	data        = hdulist[1].data		
	data2       = hdulist[2].data	
	date 		= df.iloc[i][0][8:18]
	if date == date0:
		CCF_daily_sum 	+= data[69,:]
		eCCF_daily_sum 	+= data2[69,:]**2
	else:
		# eCCF_daily[:, k]= eCCF_daily_sum[idx]**0.5 / np.mean(CCF_daily_sum[~idx])
		# CCF_daily[:, k] = 1 - CCF_daily_sum[idx] / np.mean(CCF_daily_sum[~idx])
		eCCF_daily[:, k]= eCCF_daily_sum**0.5 / np.mean(CCF_daily_sum[~idx])
		CCF_daily[:, k] = 1 - CCF_daily_sum / np.mean(CCF_daily_sum[~idx])
		rv_daily[k] 	= np.average(rv[i0:i], weights=1/rv_err[i0:i]**2)
		rv_raw_daily[k] = np.average(rv_raw[i0:i], weights=1 / rv_err[i0:i] ** 2)
		erv_daily[k] 	= (1/sum(1/rv_err[i0:i]**2))**0.5
		bjd_daily[k] 	= np.average(bjd[i0:i], weights=1/rv_err[i0:i]**2)
		widening_var_daily[k] = np.average(widening_var[i0:i], weights=1/rv_err[i0:i]**2)

		bis_daily[k] 	= np.average(df["bis_span"][i0:i], weights=1/df["bis_span_err"][i0:i]**2)
		fwhm_daily[k]	= np.average(df["fwhm"][i0:i], weights=1/df["fwhm_err"][i0:i]**2)
		ebis_daily[k] 	= (1 / sum(1 / df["bis_span_err"][i0:i] ** 2)) ** 0.5
		efwhm_daily[k] 	= (1 / sum(1 / df["fwhm_err"][i0:i] ** 2)) ** 0.5

		print(k, i0, i, date)
		k += 1
		i0 = i
		date0 = date
		CCF_daily_sum 	= data[69,:]
		eCCF_daily_sum 	= data2[69,:]**2
		
eCCF_daily[:, k] 		= eCCF_daily_sum**0.5 / np.mean(CCF_daily_sum[~idx])
CCF_daily[:, k] 		= 1 - CCF_daily_sum / np.mean(CCF_daily_sum[~idx])
# eCCF_daily[:, k] 		= eCCF_daily_sum[idx]**0.5 / np.mean(CCF_daily_sum[~idx])
# CCF_daily[:, k] 		= 1 - CCF_daily_sum[idx] / np.mean(CCF_daily_sum[~idx])
rv_daily[k] 			= np.average(rv[i0:i+1], weights=1/rv_err[i0:i+1]**2)
rv_raw_daily[k] 		= np.average(rv_raw[i0:i+1], weights=1/rv_err[i0:i+1]**2)
erv_daily[k] 			= (1/sum(1/rv_err[i0:i+1]**2))**0.5
bjd_daily[k] 			= np.average(bjd[i0:i+1], weights=1/rv_err[i0:i+1]**2)
widening_var_daily[k] 	= np.average(widening_var[i0:i+1], weights=1/rv_err[i0:i+1]**2)

bis_daily[k] 			= np.average(df["bis_span"][i0:i+1], weights=1 / df["bis_span_err"][i0:i+1] ** 2)
fwhm_daily[k] 			= np.average(df["fwhm"][i0:i+1], weights=1 / df["fwhm_err"][i0:i+1] ** 2)
ebis_daily[k] 			= (1/sum(1/df["bis_span_err"][i0:i+1]**2))**0.5
efwhm_daily[k] 			= (1/sum(1/df["fwhm_err"][i0:i+1]**2))**0.5

idx_0 		= (rv_daily==0)
rv_daily 	= rv_daily[~idx_0]
rv_raw_daily= rv_raw_daily[~idx_0]
erv_daily 	= erv_daily[~idx_0]
bjd_daily 	= bjd_daily[~idx_0]
CCF_daily 	= CCF_daily[:, ~idx_0]
eCCF_daily 	= eCCF_daily[:, ~idx_0]

bis_daily 	= bis_daily[~idx_0]
fwhm_daily  = fwhm_daily[~idx_0]
ebis_daily 	= ebis_daily[~idx_0]
efwhm_daily = efwhm_daily[~idx_0]

# Because files with high diff_extinction (>0.1) are excluded, 
# the number of daily averages becomes 567!
# In [64]: rv_daily.shape
# Out[64]: (567,)

#----------------------------------
# then go to FIESTA
#----------------------------------

# The following are tests regarding part 2
 
idx_ccf = (V_grid<=15.58) & (V_grid>-16.40)
plt.plot(V_grid[idx_ccf], 1-CCF_daily[idx_ccf,:])
plt.title('CCF (overplotted)')
plt.xlabel('V grid [km/s]')
plt.ylabel('Normalised flux')
plt.savefig('CCF.png')
plt.show()

if 0:
	for i in range(len(bjd_daily)):
		plt.plot(V_grid, CCF_daily[:,i]-CCF_daily[:,0], alpha=0.1)
	plt.xlabel('V_grid')
	plt.ylabel('delta_CCF')
	plt.savefig('delta_ccf.png')
	plt.show()

	sigma_daily 		= np.zeros(632)
	V_centre_daily 	= np.zeros(632)

	for n in range(632):
		# print(n)
		popt, pcov 	= curve_fit(gaussian, V_grid, CCF_daily[:,n], p0=[0.5, (max(V_grid)+min(V_grid))/2, 1, 0])
		sigma_daily[n] 	= popt[2]
		V_centre_daily[n] = popt[1]


	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd_daily, sigma_daily, '.', alpha=0.5)
	plt.plot(bjd_daily, sigma_daily, '-', alpha=0.2)
	plt.xlabel('bjd')
	plt.ylabel('sigma_daily')
	plt.savefig('sigma_daily.png')
	plt.show()
	plt.close()

	#overplot
	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, fwhm_raw/2.355/1000, '.', alpha=0.2, label='fwhm_raw_scaled')
	plt.plot(bjd_daily, sigma_daily, '.', alpha=0.5, label='sigma_daily')
	plt.plot(bjd_daily, sigma_daily, '-')
	plt.legend()
	plt.xlabel('bjd')
	plt.ylabel('sigma_daily')
	plt.savefig('sigma_daily_vs_fwhm_raw.png')
	# plt.show()
	plt.close()


	# How much the line profile should be widened? 
	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, (df['fwhm']**2 - fwhm_raw**2)**0.5 / 1000, '.', alpha=0.5)
	plt.xlabel('bjd')
	plt.ylabel('widened [km/s]')
	plt.savefig('widened.png')
	plt.close()


	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, widening_var, '.', alpha=0.5, label='widening_var')
	plt.plot(bjd_daily, widening_var_daily, '.', alpha=0.5, label='widening_var_daily')
	plt.legend()
	plt.xlabel('bjd')
	plt.ylabel('widening_fwhm_var [km/s]^2')
	plt.savefig('widened_comparison.png')
	plt.close()


	widening_sigma2_daily = widening_var_daily / 2.355**2


	CCF_daily_centre = np.zeros(632)
	c_correction = np.zeros(632)
	for n in range(632):
		popt, pcov 	= curve_fit(gaussian, V_grid, CCF_daily[:,n], p0=[0.5, (max(V_grid)+min(V_grid))/2, 1, 0])
		CCF_daily_centre[n] = popt[1]
		c_correction[n] = popt[3]


	if 0: # testing the broading of a profile 
		# def gaussian(x, amp, mu, sig, c):
		#     return amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) + c

		broadening_profile0 = gaussian(V_grid, 1, 0, widening_sigma2_daily[0]**0.5, 0)

		plt.plot(V_grid, CCF_daily[:,0], label='CCF_daily0')
		plt.plot(V_grid, broadening_profile, label='win')
		plt.plot(V_grid, np.convolve(CCF_daily[:,0], broadening_profile0, 'same'), label='broadened_profile')
		plt.legend()
		plt.xlabel('V_grid')
		plt.savefig('CCF_broadening.png')
		plt.close()


		from scipy import signal

		V_grid_inter 	= (np.arange(4801) - 2400) * (0.82/100)
		f = interp1d(V_grid, CCF_daily[:,0], kind='cubic')
		CCF_inter = f(V_grid_inter)
		broadening_profile = gaussian(V_grid_inter, 1, 0, widening_sigma2_daily[0]**0.5, 0)

		plt.plot(V_grid_inter, CCF_inter, label='CCF_daily0')
		plt.plot(V_grid_inter, broadening_profile, label='win')
		broadened_profile = signal.convolve(CCF_inter, broadening_profile,mode='same') / sum(broadening_profile)
		plt.plot(V_grid_inter, broadened_profile, label='broadened_profile')
		plt.legend()
		plt.xlabel('V_grid')
		plt.savefig('CCF_broadening_inter.png')
		plt.close()


		if 0:
			V_grid_inter 	= (np.arange(4801) - 2400) * (0.82/100)
			f = interp1d(V_grid, CCF_daily[:,0], kind='cubic')
			CCF_inter = f(V_grid_inter)
			broadening_profile = gaussian(V_grid_inter, 1, 0, widening_sigma2_daily[0]**0.5, 0)

			plt.plot(V_grid_inter, CCF_inter, label='CCF_daily0')
			plt.plot(V_grid_inter, broadening_profile, label='broadening_profile')
			plt.plot(V_grid_inter, np.convolve(CCF_inter, broadening_profile, 'same') / sum(broadening_profile), label='broadening_profile')
			plt.legend()
			plt.xlabel('V_grid')
			plt.savefig('CCF_broadening_inter.png')
			plt.close()



	CCF_daily_corrected = np.zeros(CCF_daily.shape)
	from scipy import signal
	V_grid_inter = (np.arange(4801) - 2400) * (0.82/100)

	for n in range(632):
		f 					= interp1d(V_grid, CCF_daily[:,n], kind='cubic')
		CCF_inter 			= f(V_grid_inter)
		win 				= gaussian(V_grid_inter, 1, 0, widening_sigma2_daily[n]**0.5, 0)
		broadened_profile 	= signal.convolve(CCF_inter, win, mode='same') / sum(win)
		CCF_daily_corrected[:, n] = broadened_profile[np.arange(49)*100]


	# check 
	sigma_daily_corrected = np.zeros(sigma_daily.shape)
	V_centre_daily_corrected = np.zeros(sigma_daily.shape)
	for n in range(632):
		popt, pcov 	= curve_fit(gaussian, V_grid, CCF_daily_corrected[:,n], p0=[0.5, (max(V_grid)+min(V_grid))/2, 1, 0])
		sigma_daily_corrected[n] 	= popt[2]
		V_centre_daily_corrected[n] = popt[1]

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd_daily, sigma_daily_corrected / sigma_daily, '.', alpha=0.5, label='ratio_sigma')
	plt.plot(bjd, df['fwhm'] / fwhm_raw, '.', alpha=0.5, label='ratio_fwhm')
	plt.legend()
	plt.xlabel('bjd_daily')
	plt.savefig('ratio_check.png')
	plt.close()




	fig, axes = plt.subplots(figsize=(12, 4))
	# plt.plot(bjd, rv_raw, '.', alpha=0.2, label='rv raw')
	plt.plot(bjd_daily, CCF_daily_centre*1000, '.', alpha=0.2, label='CCF_daily_rv')
	plt.plot(bjd_daily, V_centre_daily_corrected*1000, '.', alpha=0.2, label='CCF_daily_corrected')
	plt.legend()
	plt.xlabel('date_bjd')
	plt.ylabel('rv [m/s]')
	plt.savefig('rv_comparison.png')
	plt.show()
	plt.close()

	fig, axes = plt.subplots(figsize=(12, 4))
	# plt.plot(bjd, rv_raw, '.', alpha=0.2, label='rv raw')
	plt.plot(bjd_daily, (CCF_daily_centre-V_centre_daily_corrected)*1000, '.', alpha=0.2, label='CCF_daily_rv_diff')
	plt.plot(bjd_daily, widening_sigma2_daily**0.5, '.', alpha=0.2, label='widening_sigma2_daily**0.5')
	plt.legend()
	plt.xlabel('date_bjd')
	plt.ylabel('rv_diff [m/s]')
	plt.savefig('rv_comparison_difference.png')
	plt.show()
	plt.close()


	# =============================================================================
	# # analyse a subset of the time series 
	# =============================================================================
	# let's take segment #2
	if 0:
		idx2 = (bjd>57475) & (sigma < 3.47) & (bjd<57500)
		idx1 = (bjd>57315) & (bjd<bjd[idx2][0]) 

		fig, axes = plt.subplots(figsize=(12, 4))
		# plt.plot(bjd[idx1], sigma[idx1], '.', alpha=0.5)
		# plt.plot(bjd[idx2], sigma[idx2], '.', alpha=0.5)
		plt.plot(bjd[idx1], sigma[idx1], '.', alpha=0.5)
		plt.plot(bjd[idx1], reg.coef_ * bjd[idx1] + reg.intercept_, '-')
		plt.xlabel('bjd')
		plt.savefig('simga_s2.png')
		plt.close()


		from sklearn.linear_model import LinearRegression
		reg = LinearRegression().fit(bjd[idx1].reshape(-1, 1), sigma[idx1])


		sigma_new = sigma[idx1] - (reg.coef_ * bjd[idx1] + reg.intercept_) + np.max(sigma[idx1])
		fig, axes = plt.subplots(figsize=(12, 4))
		plt.plot(bjd[idx1], sigma_new, '.', alpha=0.5)
		# plt.plot(bjd[idx1], reg.coef_ * bjd[idx1] + reg.intercept_, '-')
		plt.xlabel('bjd')
		plt.savefig('simga_s2_corrected.png')
		plt.close()


		# update the width of CCF 

		widening_sigma = (sigma_new**2 - sigma[idx1]**2 - min(sigma_new**2 - sigma[idx1]**2) +0.01)**0.5
		CCF_corrected = np.zeros(CCF[:,idx1].shape)

		for n in range(len(bjd[idx1])):
			f 					= interp1d(V_grid, CCF[:, idx1][:,n], kind='cubic')
			CCF_inter 			= f(V_grid_inter)
			win 				= gaussian(V_grid_inter, 1, 0, widening_sigma[n], 0)
			broadened_profile 	= signal.convolve(CCF_inter, win, mode='same') / sum(win)
			CCF_corrected[:, n] = broadened_profile[np.arange(49)*100]

	idx_t =  (-15<V_grid) & (V_grid<15)
	CCF_daily_corrected_trunked = CCF_daily_corrected[idx_t, :]


#==============================================================================
# Feed CCFs into FIESTA
#==============================================================================
# eCCF = np.zeros(CCF_daily_corrected.shape)
# eCCF = np.zeros(CCF_daily.shape)

if 0:
	CCF_daily = CCF_daily[:, 0::5]
	eCCF_daily = eCCF_daily[:, 0::5]
	bjd_daily = bjd_daily[0::5]
	rv_daily = rv_daily[0::5]
	rv_raw_daily = rv_raw_daily[0::5]
	erv_daily = erv_daily[0::5]

if 0:
	np.savetxt('V_grid.txt', V_grid)
	np.savetxt('CCF_daily.txt', CCF_daily)
	np.savetxt('eCCF_daily.txt', eCCF_daily)
	np.savetxt('bjd_daily.txt', bjd_daily)
	np.savetxt('rv_daily.txt', rv_daily)
	np.savetxt('rv_raw_daily.txt', rv_raw_daily)
	np.savetxt('erv_daily.txt', erv_daily)

	np.savetxt('bis_daily.txt', bis_daily)
	np.savetxt('fwhm_daily.txt', fwhm_daily)
	np.savetxt('ebis_daily.txt', ebis_daily)
	np.savetxt('efwhm_daily.txt', efwhm_daily)

V_grid = np.loadtxt('V_grid.txt')
CCF_daily = np.loadtxt('CCF_daily.txt')
eCCF_daily = np.loadtxt('eCCF_daily.txt')
bjd_daily = np.loadtxt('bjd_daily.txt')
rv_daily = np.loadtxt('rv_daily.txt')
rv_raw_daily = np.loadtxt('rv_raw_daily.txt')
erv_daily = np.loadtxt('erv_daily.txt')

bis_daily = np.loadtxt('bis_daily.txt')
fwhm_daily = np.loadtxt('fwhm_daily.txt')
ebis_daily = np.loadtxt('ebis_daily.txt')
efwhm_daily = np.loadtxt('efwhm_daily.txt')

# read from the 5 years solar data
import pickle
with open('./SCALPELS/CCFdata_daily_15min.pkl','rb') as ccfile:
    ccfpkl,errpkl,badpkl,bjdpkl,rvbpkl,rvhpkl,rvepkl,bervpkl,sn60pkl,qualpkl,vcorpkl,btohpkl,vel = pickle.load(ccfile)

# Plot the RVs 
plt.rcParams.update({'font.size': 14})
alpha=0.5
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

fig, axes = plt.subplots(figsize=(15, 3))
plt.gcf().subplots_adjust(bottom=0.2)
plt.errorbar(bjd_daily, rv_daily-np.mean(rv_daily), erv_daily, c='purple', marker='.', ls='none', alpha= 0.3, label='3 years RV')
plt.errorbar(bjdpkl[rowmask]-2400000 , (rvhpkl[rowmask]-np.mean(rvhpkl[rowmask]))*1000+10, rvepkl[rowmask]*1000, c='black', marker='.', ls='none', alpha= 0.3, label='5 years RV')
plt.legend()
plt.xlabel('BJD - 2400000 [d]')
plt.ylabel('RV [m/s]')
plt.savefig('rv_daily.pdf')
plt.show()

# some statistics 
if 0:
	bjdpkl = bjdpkl-2400000
	len(bjdpkl[rowmask][bjdpkl[rowmask]<58100]) #431
	len(bjd_daily[bjd_daily<58100]) #469

# change the parameter names 
rowmask		= np.load('./SCALPELS/rowmask.npy')
CCF_daily 	= 1-ccfpkl[rowmask,:].T
eCCF_daily 	= errpkl[rowmask,:].T # note that noise is overestaimated in Cameron's data
eCCF_daily 	= eCCF_daily / (np.mean(eCCF_daily) / np.mean(np.loadtxt('eCCF_daily.txt')))
bjd_daily 	= bjdpkl[rowmask]-2400000
rv_daily 	= rvhpkl[rowmask]*1000
rv_raw_daily= rvbpkl[rowmask]*1000
erv_daily 	= rvepkl[rowmask]*1000


if 0:
	# plt.plot(V_grid, CCF_daily - np.mean(CCF_daily, axis=1), '.')
	# plt.show()

	plt.plot(V_grid, 1-CCF_daily, 'k')
	# plt.plot(V_grid, eCCF_daily, 'k')
	plt.show()

	plt.plot(V_grid, np.loadtxt('eCCF_daily.txt'), 'b')
	plt.show()

	plt.plot(V_grid, ccfpkl[rowmask,:].T, 'b')
	plt.show()

	plt.plot(V_grid, errpkl[rowmask,:].T, 'b')
	plt.show()


# normality tests
if 0:
	shift_spectrum, err_shift_spectrum, power_spectrum, err_power_spectrum, RV_gauss = FIESTA(V_grid, CCF_daily, eCCF_daily)

	noise_power_spectrum = noise_power_spectrum * 1000
	for k in range(noise_power_spectrum.shape[1]):
		plt.hist(noise_power_spectrum[:,k], bins = 30)
		plt.savefig('histogram_power_spectrum_' + str(k+1) + '.png')
		plt.close()

	noise_shift_spectrum = noise_shift_spectrum * 1000
	for k in range(noise_shift_spectrum.shape[1]):
		plt.hist(noise_shift_spectrum[:,k], bins = 50)
		plt.savefig('histogram_shift_spectrum_' + str(k+1) + '.png')
		plt.close()

	#----------------------------------
	# Normality Tests
	#----------------------------------
	# refer to https://machinelearningmastery.com/a-gentle-introduction-to-normality-tests-in-python/
	from scipy.stats import shapiro # Shapiro-Wilk Test
	from scipy.stats import normaltest # D’Agostino’s K^2 Test
	from scipy.stats import anderson # Anderson-Darling Test

	for k in range(noise_shift_spectrum.shape[1]):
		data = noise_shift_spectrum[:,k]
		alpha = 0.05

		stat, p = shapiro(data)
		if p < alpha:
			print('Shapiro-Wilk Test: not Gaussian')
		# if p > alpha:
		# 	print('Sample looks Gaussian (fail to reject H0)')
		# else:
		# 	print('Sample does not look Gaussian (reject H0)')

		stat, p = normaltest(data)
		if p < alpha:
			print('D’Agostino’s K^2 Test: not Gaussian')
		# if p > alpha:
		# 	print('Sample looks Gaussian (fail to reject H0)')
		# else:
		# 	print('Sample does not look Gaussian (reject H0)')

		result = anderson(data)
		# print('Statistic: %.3f' % result.statistic)
		p = 0
		for i in range(len(result.critical_values)):
			sl, cv = result.significance_level[i], result.critical_values[i]
			if result.statistic > result.critical_values[i]:
				print(str(k) + ' Anderson-Darling Test: not Gaussian')
			# if result.statistic < result.critical_values[i]:
			# 	print('%.3f: %.3f, data looks normal (fail to reject H0)' % (sl, cv))
			# else:
			# 	print('%.3f: %.3f, data does not look normal (reject H0)' % (sl, cv))



#----------------------------------
# testing GP for filtering 
#----------------------------------


# example 
gp, st, lt 		= long_short_divide(x=bjd_daily, y=my_pca_score[:,0], yerr=my_err_pca_score[:,0], r=100)
x_pred 			= np.linspace(min(bjd_daily), max(bjd_daily), 1000)
pred, pred_var 	= gp.predict(my_pca_score[:,0], x_pred, return_var=True)

plt.fill_between(x_pred, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var), color="k", alpha=0.2)
plt.plot(x_pred, pred, "k", lw=1.5, alpha=0.5)
plt.errorbar(bjd_daily, my_pca_score[:,0], yerr=my_err_pca_score[:,0], fmt=".k", capsize=0)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# residual
y_pred, y_pred_var = gp.predict(y, x, return_var=True)
plt.errorbar(x, y_pred-y, yerr=yerr, fmt=".k", capsize=0)
plt.show()


for i in range(3):

plt.rcParams.update({'font.size': 12})

def time_series(x, y, dy, N=None,
				ylabel='k=',
				title='Time series',
				file_name='Time_series.png'):
	if N==None:
		N = y.shape[1]
	plt.subplots(figsize=(12, N))

	for i in range(N):
		ax = plt.subplot(N, 1, i+1)
		if i == 0:
			plt.title(title)
		plt.errorbar(x, y[:, i], dy[:, i], marker='.', ls='none', alpha=0.5, ms=5)
		plt.ylabel(ylabel+str(i+1))
		if i != N-1:
			ax.set_xticks([])
		else:
			plt.xlabel('BJD - 2400000 [d]')
	plt.savefig(file_name)
	plt.show()

gp, st, lt 		= long_short_divide(x=bjd_daily, y=my_pca_score[:,0], yerr=my_err_pca_score[:,0], r=100)
x_pred 			= np.linspace(min(bjd_daily), max(bjd_daily), 1000)
pred, pred_var 	= gp.predict(my_pca_score[:,0], x_pred, return_var=True)
bjd_daily_pred, bjd_daily_pred_var 	= gp.predict(my_pca_score[:,0], bjd_daily, return_var=True)

plt.subplots(figsize=(12, 3))
ax = plt.subplot(3, 1, 1)
plt.errorbar(bjd_daily, my_pca_score[:,0]-bjd_daily_pred, yerr=my_err_pca_score[:,0], fmt=".k", capsize=0, alpha=0.5)
ax.set_xticks([])
ax = plt.subplot(3, 1, 2)
plt.fill_between(x_pred, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var), color="k", alpha=0.2)
plt.plot(x_pred, pred, "k", lw=1.5, alpha=0.5)
ax.set_xticks([])
ax = plt.subplot(3, 1, 3)
gp, st, lt 		= long_short_divide(x=bjd_daily, y=my_pca_score[:,1], yerr=my_err_pca_score[:,1], r=100)
pred, pred_var 	= gp.predict(my_pca_score[:,1], x_pred, return_var=True)
plt.fill_between(x_pred, pred - np.sqrt(pred_var), pred + np.sqrt(pred_var), color="k", alpha=0.2)
plt.plot(x_pred, pred, "k", lw=1.5, alpha=0.5)
plt.show()



std_matrix = np.zeros(15)

k_max = 7
for k_max in (np.arange(16)+5):
	df, shift_spectrum, err_shift_spectrum, power_spectrum, err_power_spectrum, RV_gauss = FIESTA(V_grid, CCF_daily, eCCF_daily, k_max=k_max)
	# Convertion from km/s to m/s
	shift_spectrum 		*= 1000
	err_shift_spectrum 	*= 1000
	RV_gauss 			*= 1000

	shift_function 	= np.zeros(shift_spectrum.shape)
	# short_variation = np.zeros(shift_spectrum.shape)
	# long_variation 	= np.zeros(shift_spectrum.shape)

	for i in range(shift_spectrum.shape[0]):
		shift_function[i,:] = shift_spectrum[i,:] - rv_raw_daily # look back; change rv_raw_daily

		# _, short_variation[i,:], long_variation[i,:] = long_short_divide(
		# 	x=bjd_daily, y=shift_function[i,:], yerr=err_shift_spectrum[i,:], r=100)

	if 0: # Compute ICA		
		ica = FastICA(n_components=3)
		S_ = ica.fit_transform(shift_function.T)  # Reconstruct signals
		time_series(x=bjd_daily, y=S_, dy=err_shift_spectrum[0:3,:].T*0, ylabel='ICA',
						title='ICA',
						file_name='ICA.png')	
		periodogram(x=bjd_daily, y=S_, dy=err_shift_spectrum[0:3,:].T,
					plot_min_t=2, study_min_t=5, max_f=1, spp=100,
					ylabel='ICA',
					title='ICA periodogram',
					file_name='ICA_periodogram.png')

		# ICA on PCA
		ica = FastICA(n_components=3)
		S_ = ica.fit_transform(my_pca_score[:,0:4])  # Reconstruct signals
		time_series(x=bjd_daily, y=S_, dy=my_pca_score[:,0:3]*0, ylabel='ICA',
						title='ICA_pca',
						file_name='ICA_pca.png')	
		periodogram(x=bjd_daily, y=S_, dy=my_err_pca_score[:,0:3],
					plot_min_t=2, study_min_t=5, max_f=1, spp=100,
					ylabel='ICA',
					title='ICA_pca periodogram',
					file_name='ICA_pca_periodogram.png')

	if 0:
		plt.rcParams.update({'font.size': 14})
		plot_all(k_mode=20, t=bjd_daily, rv=rv_daily, erv=erv_daily, 
			ind=power_spectrum, eind=err_power_spectrum, 
			ts_xlabel='BJD - 2400000 [d]', 
			rv_xlabel='$RV_{HARPS}$', 
			pe_xlabel='Period [days]',
			ind_yalbel=r'$A$',
			file_name='Amplitude_time-series_correlation_periodogram_SCALPELS.pdf')
			# file_name='Amplitude_time-series_correlation_periodogram.pdf')
		plt.show()

		plot_all(k_mode=20, t=bjd_daily, rv=rv_daily, erv=erv_daily, 
			ind=shift_function, eind=err_shift_spectrum, 
			ts_xlabel='BJD - 2400000 [d]', 
			rv_xlabel='$RV_{HARPS}$', 
			pe_xlabel='Period [days]',
			ind_yalbel=r'$\Delta RV$',
			# file_name='shift_time-series_correlation_periodogram_SCALPELS.pdf')
			file_name='shift_time-series_correlation_periodogram.pdf')
		plt.show()

		# --------------------------------- #
		# divide the time-series in the gap #
		# --------------------------------- #
		idx_bjd = bjd_daily<58100
		plot_all(k_mode=20, t=bjd_daily[idx_bjd], rv=rv_daily[idx_bjd], erv=erv_daily[idx_bjd], 
			ind=shift_function[:,idx_bjd], eind=err_shift_spectrum[:,idx_bjd], 
			ts_xlabel='BJD - 2400000 [d]', 
			rv_xlabel='$RV_{HARPS}$', 
			pe_xlabel='Period [days]',
			ind_yalbel=r'$\Delta RV$',
			file_name='shift_time-series_correlation_periodogram_SCALPELS1.pdf')
		plt.show()

		plot_all(k_mode=20, t=bjd_daily[~idx_bjd], rv=rv_daily[~idx_bjd], erv=erv_daily[~idx_bjd], 
			ind=shift_function[:,~idx_bjd], eind=err_shift_spectrum[:,~idx_bjd], 
			ts_xlabel='BJD - 2400000 [d]', 
			rv_xlabel='$RV_{HARPS}$', 
			pe_xlabel='Period [days]',
			ind_yalbel=r'$\Delta RV$',
			file_name='shift_time-series_correlation_periodogram_SCALPELS2.pdf')
		plt.show()	

		# hack!
		# shift_function = short_variation
		# shift_function = long_variation

	my_P, my_pca_score, my_err_pca_score, n_pca = my_pca(X=shift_function.T, X_err=err_shift_spectrum.T, nor=True)

	if 0: # temporarily skip
		plot_all(k_mode=3, t=bjd_daily, rv=rv_daily, erv=erv_daily, 
			ind=my_pca_score.T, eind=my_err_pca_score.T, 
			ts_xlabel='BJD - 2400000 [d]', 
			rv_xlabel='$RV_{HARPS}$', 
			pe_xlabel='Period [days]',
			ind_yalbel='PC',
			file_name='PCA_delta_RV_k_max={:d}.pdf'.format(k_max))
		plt.close()

	if 0:
		my_P, my_pca_score, my_err_pca_score, n_pca = my_pca(X=power_spectrum.T, X_err=err_power_spectrum.T, nor=True)

		plot_all(k_mode=3, t=bjd_daily, rv=rv_daily, erv=erv_daily, 
			ind=my_pca_score.T, eind=my_err_pca_score.T, 
			ts_xlabel='BJD - 2400000 [d]', 
			rv_xlabel='$RV_{HARPS}$', 
			pe_xlabel='Period [days]',
			ind_yalbel='PC',
			file_name='PCA_A_k_max={:d}.pdf'.format(k_max))
		plt.close()

	k_feature = 3
	short_variation = np.zeros((k_feature, shift_spectrum.shape[1]))
	long_variation 	= np.zeros((k_feature, shift_spectrum.shape[1]))

	for i in range(k_feature):
		_, short_variation[i,:], long_variation[i,:] = long_short_divide(
			x=bjd_daily, y=my_pca_score[:,i], yerr=my_err_pca_score[:,i], r=100)

	if 0:

		plot_all_but_corr(k_mode=3, t=bjd_daily,
			ind=short_variation, eind=my_err_pca_score[:,0:3].T, 
			ts_xlabel='BJD - 2400000 [d]', 
			height_ratio = 0.8,
			pe_xlabel='Period [days]',
			ind_yalbel=['S-PC'],
			file_name='Rotation_modulated_variation.pdf')
		plt.close()

		bgls_periodogram(x=bjd_daily, y=power_spectrum.T, dy=err_power_spectrum.T, N=None,
					plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=None,
					ylabel='A',
					title = 'Periodogram',
					file_name='Periodogram.png')

		plot_all_but_corr(k_mode=3, t=bjd_daily,
			ind=long_variation, eind=my_err_pca_score[:,0:3].T, 
			ts_xlabel='BJD - 2400000 [d]', 
			height_ratio = 0.8,
			pe_xlabel='Period [days]',
			ind_yalbel=['L-PC'],
			file_name='Long-term_variation.pdf')
		plt.close()

		# further combined 
		plot_all_but_corr(k_mode=6, t=bjd_daily,
			ind=np.vstack([short_variation, long_variation]), 
			eind=np.vstack([my_err_pca_score[:,0:3].T, my_err_pca_score[:,0:3].T]),
			ts_xlabel='BJD - 2400000 [d]', 
			height_ratio = 0.8,
			pe_xlabel='Period [days]',
			ind_yalbel=['S-PC$_1$','S-PC$_2$','S-PC$_3$','L-PC$_1$','L-PC$_2$','L-PC$_3$'],
			file_name='LS-term_variation.pdf')
		plt.close()


	if 0:
		rot_P, rot_pca_score, rot_err_pca_score, n_pca = my_pca(X=short_variation[0:3,:].T, X_err=my_err_pca_score[:,0:3], nor=True)

		time_series(x=bjd_daily, y=rot_pca_score, dy=rot_err_pca_score, N=None, ylabel='PCA',
						title='PCA on rotation modulated variation time series',
						file_name='PCA_rotation_modulated_variation_time_series.png')

		periodogram(x=bjd_daily, y=rot_pca_score, dy=rot_err_pca_score, N=None,
					plot_min_t=2, study_min_t=5, max_f=1, spp=100,
					ylabel='PCA',
					title='Rotation modulated variation time series PCA periodogram',
					file_name='Rotation_modulated_variation_time_series_PCA_periodogram.png')		

		plt.plot(short_variation[0,:], rot_pca_score[:,0].T, '.')
		plt.xlabel('PCA1')
		plt.ylabel('new PCA1')
		plt.show()

		# hack #
		my_pca_score = np.vstack((long_variation[0:3,:], rot_pca_score[:,0].reshape(1,567))).T
		my_err_pca_score = np.vstack((my_err_pca_score[:,0:3].T, rot_err_pca_score[:,0].reshape(1,567))).T

		# my_P, my_pca_score, my_err_pca_score, n_pca = my_pca(X=power_spectrum.T, X_err=err_power_spectrum.T, nor=False)

		# plot_all(k_mode=6, t=bjd_daily, rv=rv_daily, erv=erv_daily, 
		# 	ind=my_pca_score.T, eind=my_err_pca_score.T, 
		# 	ts_xlabel='BJD - 2400000 [d]', 
		# 	rv_xlabel='$RV_{HARPS}$', 
		# 	pe_xlabel='Period [days]',
		# 	ind_yalbel='$\Delta RV_k$',
		# 	file_name='PCA_delta_RV_k_max={:d}.png'.format(k_max))

		bjd_daily 		= np.loadtxt('bjd_daily.txt')
		
		n_pca = 4
		plot_all(k_mode=4, t=bjd_daily, rv=rv_daily, erv=erv_daily, 
			ind=my_pca_score.T, eind=my_err_pca_score.T, 
			ts_xlabel='BJD - 2400000 [d]', 
			rv_xlabel='$RV_{HARPS}$', 
			pe_xlabel='Period [days]',
			ind_yalbel='PCA',
			file_name='PCA_long_short_trend.png')
		plt.show()


	#---------------------------#
	# Multiple Regression Model 
	#---------------------------#
	if 0: # temporarily skip
		# Model 1 #
		feature_matrix = my_pca_score[:,0:3] #(567,3)
		y_hat1, w_std_all, w_rms1, score1, df1 = mlr(feature_matrix, target_vector=rv_daily, etarget_vector=erv_daily, lag='False')

	# Model 2 #
	feature_matrix = np.vstack([short_variation, long_variation]).T #(567,6)
	y_hat2, w_std_all, w_rms2, score2, df2 = mlr(feature_matrix, target_vector=rv_daily, etarget_vector=erv_daily, lag='False')

	if 0: # temporarily skip
		# Model 4 #
		fwhm_bis 	= np.vstack((fwhm_daily, bis_daily)).T
		from sklearn.preprocessing import StandardScaler
		scaler 		= StandardScaler()
		fwhm_bis 	= scaler.fit_transform(fwhm_bis)
		efwhm_bis 	= np.vstack((efwhm_daily/np.std(fwhm_daily), ebis_daily/np.std(bis_daily))).T
		y_hat4, w_std_all, w_rms4, score4, df4 = mlr(fwhm_bis, target_vector=rv_daily, etarget_vector=erv_daily, lag='False')

		# Model 5 #
		k_feature 		= fwhm_bis.shape[1]
		short_variation = np.zeros((k_feature, fwhm_bis.shape[0]))
		long_variation 	= np.zeros((k_feature, fwhm_bis.shape[0]))

		for i in range(k_feature):
			_, short_variation[i,:], long_variation[i,:] = long_short_divide(
				x=bjd_daily, y=fwhm_bis[:,i], yerr=efwhm_bis[:,i], r=100)

		feature_matrix = np.vstack([short_variation, long_variation]).T #(567,4)
		y_hat5, w_std_all, w_rms5, score5, df5 = mlr(feature_matrix, target_vector=rv_daily, etarget_vector=erv_daily, lag='False')

	#----------------
	if 0: # visually check how the daily binned observation hours are distributed 
		time = bjd_daily - [int(bjd_daily[i]+0.5) for i in np.arange(len(bjd_daily))] + 0.5
		len(time[abs((time-np.median(time))*24)<=1]) / len(time) # 92%
		len(time[abs((time-np.median(time))*24)<=2]) / len(time) # 99%
		# no need to interpolate the daily binned RV 
		_ = plt.hist((time)*24, bins=20, color='black', alpha=0.5)
		plt.savefig('observation_hour.png')
		plt.title('Daily binned observation hour')
		plt.show()
	#----------------

	#---------------------------------------------------------------------------------#
	# Multiple Linear Regression (only interpolate the indicators)
	#---------------------------------------------------------------------------------#
	from scipy.interpolate import CubicSpline

	day = 5
	# bjd_daily 		= np.loadtxt('bjd_daily.txt')
	idx_bjd 		= bjd_daily<58100
	bjd_daily_part1 = bjd_daily[idx_bjd]
	bjd_daily_part2 = bjd_daily[~idx_bjd]

	t1_min = min(bjd_daily_part1) 	#= 57233.05440362564
	t1_max = max(bjd_daily_part1)	#= 58068.0739790099
	t2_min = min(bjd_daily_part2) 	#= 58163.022739323365
	t2_max = max(bjd_daily_part2) 	#= 58315.979416495415

	bjd_daily_lag1 	= bjd_daily_part1[(bjd_daily_part1>=t1_min+day-0.3) & (bjd_daily_part1<=t1_max-day+0.3)]
	bjd_daily_lag2 	= bjd_daily_part2[(bjd_daily_part2>=t2_min+day-0.3) & (bjd_daily_part2<=t2_max-day+0.3)]

	bjd_daily_lag 	= np.hstack((bjd_daily_lag1, bjd_daily_lag2))

	k_feature2 	= feature_matrix.shape[1]
	k_feature 	= int(k_feature2/2)
	feature_matrix_int_lag 	= np.zeros((len(bjd_daily_lag), k_feature2*(day+1)))
	
	for i in range((2*day+1)):	
		for j in range(k_feature):
			f = interp1d(bjd_daily, feature_matrix[:,j], fill_value='extrapolate')
			feature_matrix_int_lag[:, k_feature*i+j] = f(bjd_daily_lag-day+i)

			# cs = CubicSpline(bjd_daily, feature_matrix[:,j], extrapolate=True)
			# feature_matrix_int_lag[:, k_feature*i+j] = cs(bjd_daily_lag-day+i)

	index = ((bjd_daily>=t1_min+day-0.3) & (bjd_daily<=t1_max-day+0.3)) \
	| ((bjd_daily>=t2_min+day-0.3) & (bjd_daily<=t2_max-day+0.3))
	
	feature_matrix_int_lag[:,k_feature*(2*day+1):] = feature_matrix[index, k_feature:]

	# y_hat3, w_std_all3, res_wrms, score, variance_matrix3 = mlr(feature_matrix_int_lag, target_vector=rv_daily[index], etarget_vector=erv_daily[index], feature_matrix2=feature_matrix[index])
	# imshow_matrix(variance_matrix3, file_name='fiesta_multi_coef') # working :) 

	if 0: # determine np.median(alpna_n) = 0.019 (lambda)

		n_folds = 5
		N = 100
		alpna_n = np.zeros(N)
		for n in range(N):

			size 	= feature_matrix_int_lag.shape[0]
			idx_random = shuffle(np.arange(size))
			X 		= feature_matrix_int_lag[idx_random,:]
			y 		= rv_daily[index][idx_random]
			sw 		= 1/erv_daily[index][idx_random]**2

			lasso = Lasso()
			alphas = np.logspace(-3, -1, 100)

			tuned_parameters = [{"alpha": alphas}]

			clf = GridSearchCV(lasso, tuned_parameters, cv=n_folds, scoring='neg_root_mean_squared_error', refit=True, return_train_score=True)
			clf.fit(X, y, sample_weight=sw)
			alpna_n[n] = clf.best_estimator_.alpha
			# clf.best_score_

		import scipy
		_, bins, _ = plt.hist(alpna_n, 30, density=1, alpha=0.5)
		mu, sigma = scipy.stats.norm.fit(alpna_n)
		best_fit_line = scipy.stats.norm.pdf(bins, mu, sigma)
		plt.plot(bins, best_fit_line)
		plt.axvline(np.median(alpna_n))
		plt.text(np.median(alpna_n), 10, '{:.3f}'.format(np.median(alpna_n)))
		plt.title(r'$\lambda$ histogram')
		plt.xlabel(r'$\lambda$')
		plt.savefig('histogram_alpha.png')
		plt.show()
		#----------------------------------------



	from sklearn.linear_model import Lasso
	from sklearn.model_selection import train_test_split
	from sklearn.utils import shuffle
	from sklearn.model_selection import GridSearchCV

	size 		= feature_matrix_int_lag.shape[0]
	N_folds 	= 5
	N 			= 1000
	coef 		= np.zeros((36,N))
	alphas 		= np.logspace(-4, -1, 100)

	clf_wrms 	= np.zeros(N_folds)
	CLF_WRMS 	= np.zeros(N)
	δCLF_WRMS 	= np.zeros(N)
	CLF_aWRMS 	= np.zeros(len(alphas))
	δCLF_aWRMS 	= np.zeros(len(alphas))
	
	clf_wrms_train 	= np.zeros(N_folds)
	CLF_WRMS_train 	= np.zeros(N)
	δCLF_WRMS_train 	= np.zeros(N)
	CLF_aWRMS_train 	= np.zeros(len(alphas))
	δCLF_aWRMS_train 	= np.zeros(len(alphas))

	coef_median = np.zeros((36,len(alphas)))
	count 		= np.zeros(len(alphas))

	for n_alpha in range(len(alphas)):
		alpha = alphas[n_alpha]

		for n in range(N):

			for nf in range(N_folds): 
				idx_random = shuffle(np.arange(size))
				X 		= feature_matrix_int_lag
				y 		= rv_daily[index]
				sw 		= 1/erv_daily[index]

				X_train, X_test, y_train, y_test, sw_train, sw_test = train_test_split(X, y, sw, test_size=1/N_folds, shuffle=True)

				clf = Lasso(alpha=alpha)
				clf.fit(X_train, y_train, sample_weight=sw_train)
				y_hat = clf.predict(X_test)
				y_hat_train = clf.predict(X_train)

				_, clf_wrms[nf] = weighted_avg_and_std(y_hat-y_test, weights=sw_test)
				_, clf_wrms_train[nf] = weighted_avg_and_std(y_hat_train-y_train, weights=sw_train)
				coef[:,n] = clf.coef_ 
		
			CLF_WRMS[n] 	= np.mean(clf_wrms)
			δCLF_WRMS[n] 	= np.std(clf_wrms) / np.sqrt(N_folds)
			CLF_WRMS_train[n] 	= np.mean(clf_wrms_train)
			δCLF_WRMS_train[n] 	= np.std(clf_wrms_train) / np.sqrt(N_folds)

		CLF_aWRMS[n_alpha] 	= np.mean(CLF_WRMS)
		δCLF_aWRMS[n_alpha] = np.std(CLF_WRMS)
		CLF_aWRMS_train[n_alpha] 	= np.mean(CLF_WRMS_train)
		δCLF_aWRMS_train[n_alpha] = np.std(CLF_WRMS_train)

		coef_median[:,n_alpha] = np.median(coef,axis=1)
		count[n_alpha] = sum(coef_median[:,n_alpha] != 0)

	plt.figure().set_size_inches(8, 6)
	plt.axvline(np.median(alpna_n), color='k')
	plt.rcParams['font.size'] = '16'
	# plt.title(r'$n_{fold} = $' + str(N_folds))
	plt.xlabel(r"$\lambda$")
	plt.ylabel("Residual WRMS [m/s]")
	plt.xscale('log')
	plt.plot(alphas, CLF_aWRMS, "r--", label='testing')
	plt.fill_between(alphas, CLF_aWRMS+δCLF_aWRMS, CLF_aWRMS-δCLF_aWRMS, alpha=0.2, color='r')
	plt.plot(alphas, CLF_aWRMS_train, "b--", label='training')
	plt.fill_between(alphas, CLF_aWRMS_train+δCLF_aWRMS_train, CLF_aWRMS_train-δCLF_aWRMS_train, alpha=0.2, color='b')
	# plt.savefig('CLF_CLF_aWRMS_test_train.png')
	plt.legend()
	plt.show()

	plt.figure().set_size_inches(8, 6)
	plt.semilogx(alphas, count/max(count), lw=5, alpha=0.3)
	for i in range(len(count)):
		plt.text(alphas[i], count[i]/max(count), str(int(count[i])))
	plt.axvline(np.median(alpna_n))
	plt.xscale('log')
	plt.rcParams['font.size'] = '16'
	# plt.title(r'$n_{fold} = $' + str(N_folds))
	plt.xlabel(r"$\lambda$")
	plt.semilogx(alphas, coef_median.T)
	# plt.savefig('coeff_path_.png')
	plt.show()







	from sklearn.linear_model import LassoCV
	from sklearn.linear_model import Lasso
	from sklearn.model_selection import KFold
	from sklearn.model_selection import GridSearchCV
	

	# for n_folds in np.arange(1)+5:
	# for n_folds in np.arange(10)+5:
		# print(n_folds)
		n_folds = 5
		size 	= feature_matrix_int_lag.shape[0]
		idx_random = shuffle(np.arange(size))
		X 		= feature_matrix_int_lag[idx_random,:]
		y 		= rv_daily[index][idx_random]
		sw 		= 1/erv_daily[index][idx_random]**2

		lasso = Lasso(random_state=0, max_iter=10000)
		alphas = np.logspace(-3, -1, 100)

		tuned_parameters = [{"alpha": alphas}]

		clf = GridSearchCV(lasso, tuned_parameters, cv=n_folds, scoring='neg_root_mean_squared_error', refit=True, return_train_score=True)
		clf.fit(X, y, sample_weight=sw)
		scores = -clf.cv_results_["mean_test_score"]
		scores2 = -clf.cv_results_["mean_train_score"]
		scores_std = clf.cv_results_["std_test_score"]
		scores_std2 = clf.cv_results_["std_train_score"]
		plt.figure().set_size_inches(8, 6)
		plt.rcParams['font.size'] = '16'
		plt.title(r'$n_{fold} = $' + str(n_folds))
		plt.semilogx(alphas, scores, label='testing')
		plt.semilogx(alphas, scores2, label='training')

		# plot error lines showing +/- std. errors of the scores
		std_error = scores_std / np.sqrt(n_folds)
		std_error2 = scores_std2 / np.sqrt(n_folds)

		plt.semilogx(alphas, scores + std_error, "b--")
		plt.semilogx(alphas, scores - std_error, "b--")
		plt.semilogx(alphas, scores2 + std_error2, "b--")
		plt.semilogx(alphas, scores2 - std_error2, "b--")	

		# alpha=0.2 controls the translucency of the fill color
		plt.fill_between(alphas, scores + std_error, scores - std_error, alpha=0.2)
		plt.fill_between(alphas, scores2 + std_error2, scores2 - std_error2, alpha=0.2)
		xx=clf.best_estimator_.alpha
		plt.axvline(xx)
		plt.text(xx, scores[alphas==xx], '{:.3f}'.format(xx))

		plt.ylabel("CV score +/- std error")
		plt.xlabel(r"$\lambda$")
		plt.axhline(np.min(scores), linestyle="--", color=".5")
		plt.xlim([alphas[0], alphas[-1]])
		plt.legend()

		if 0:
			lasso_cv = LassoCV(alphas=alphas, random_state=0, max_iter=10000, fit_intercept=True)
			k_fold = KFold(3)

			for k, (train, test) in enumerate(k_fold.split(X, y)):
			    lasso_cv.fit(X[train], y[train])
			    # print(X[train])
			    print(
			        "[fold {0}] alpha: {1:.5f}, score: {2:.5f}".format(
			            k, lasso_cv.alpha_, lasso_cv.score(X[test], y[test])
			        )
			    )

		plt.savefig('n_fold=' + str(n_folds) + '.png')
		# plt.close()
		plt.show()

# linspace 
		n_folds = 5
		size 	= feature_matrix_int_lag.shape[0]
		idx_random = shuffle(np.arange(size))
		X 		= feature_matrix_int_lag[idx_random,:]
		y 		= rv_daily[index][idx_random]
		sw 		= 1/erv_daily[index][idx_random]**2

		lasso = Lasso(random_state=0, max_iter=10000)
		alphas = np.linspace(0.01, 0.05, 50)
		

		tuned_parameters = [{"alpha": alphas}]

		clf = GridSearchCV(lasso, tuned_parameters, cv=n_folds, refit=True, return_train_score=True)
		clf.fit(X, y, sample_weight=sw)

		scores = clf.cv_results_["mean_test_score"]
		scores2 = clf.cv_results_["mean_train_score"]
		scores_std = clf.cv_results_["std_test_score"]
		scores_std2 = clf.cv_results_["std_train_score"]
		
		plt.figure().set_size_inches(8, 6)
		plt.rcParams['font.size'] = '16'
		plt.title(r'$n_{fold} = $' + str(n_folds))
		plt.plot(alphas, scores, label='testing')
		plt.plot(alphas, scores2, label='training')

		# plot error lines showing +/- std. errors of the scores
		std_error = scores_std / np.sqrt(n_folds)
		std_error2 = scores_std2 / np.sqrt(n_folds)

		plt.plot(alphas, scores + std_error, "b--")
		plt.plot(alphas, scores - std_error, "b--")
		plt.plot(alphas, scores2 + std_error2, "b--")
		plt.plot(alphas, scores2 - std_error2, "b--")	

		# alpha=0.2 controls the translucency of the fill color
		plt.fill_between(alphas, scores + std_error, scores - std_error, alpha=0.2)
		plt.fill_between(alphas, scores2 + std_error2, scores2 - std_error2, alpha=0.2)
		xx=clf.best_estimator_.alpha
		plt.axvline(xx)
		plt.text(xx, scores[alphas==xx], '{:.3f}'.format(xx))

		plt.ylabel("CV score +/- std error")
		plt.xlabel(r"$\lambda$")
		plt.axhline(np.max(scores), linestyle="--", color=".5")
		plt.xlim([alphas[0], alphas[-1]])
		plt.legend()
		plt.savefig('n_fold=' + str(n_folds) + '.png')
		plt.show()
		# plt.close()
		




		
		
		



	#----------------------------------------		
	N = 1
	coef = np.zeros((36,N))
	alphass = np.logspace(-3, 0, 50)
	clf_wrms = np.zeros(N)
	clf_mean = np.zeros(N)
	CLF_WRMS = np.zeros(len(alphass))
	CLF_MEAN = np.zeros(len(alphass))
	coef_median = np.zeros((36,len(alphass)))
	count = np.zeros(len(alphass))

	for n_alpha in range(len(alphass)):
		alphas = np.array([alphass[n_alpha]])
		for n in range(N):
			n_folds = 5
			size 	= feature_matrix_int_lag.shape[0]
			idx_random = shuffle(np.arange(size))
			X 		= feature_matrix_int_lag[idx_random,:]
			y 		= rv_daily[index][idx_random]
			sw 		= 1/erv_daily[index][idx_random]**2

			lasso = Lasso(random_state=0, max_iter=10000)
			
			tuned_parameters = [{"alpha": alphas}]

			clf = GridSearchCV(lasso, tuned_parameters, cv=n_folds, refit=True, return_train_score=True)
			clf.fit(X, y, sample_weight=sw)

			y_hat = clf.predict(X)
			_, clf_wrms[n] = weighted_avg_and_std(y_hat-y, weights=sw)
			coef[:,n] = clf.best_estimator_.coef_ 
		
		CLF_MEAN[n_alpha] = np.mean(clf_wrms)
		CLF_WRMS[n_alpha] = np.mean(clf_wrms) / np.sqrt(n_folds)
		coef_median[:,n_alpha] = np.median(coef,axis=1)
		count[n_alpha] = sum(coef_median[:,n_alpha] != 0)

	plt.figure().set_size_inches(8, 6)
	plt.semilogx(alphass, count/max(count), lw=5, alpha=0.3)
	for i in range(len(count)):
		plt.text(alphass[i], count[i]/max(count), str(int(count[i])))
	plt.axvline(np.median(alpna_n))
	plt.xscale('log')
	plt.rcParams['font.size'] = '16'
	plt.title(r'$n_{fold} = $' + str(n_folds))
	plt.xlabel(r"$\lambda$")
	plt.semilogx(alphass, coef_median.T)
	plt.savefig('coeff_path.png')
	plt.show()
	plt.close()

	plt.figure().set_size_inches(8, 6)
	plt.axvline(np.median(alpna_n))
	plt.xscale('log')
	plt.rcParams['font.size'] = '16'
	plt.title(r'$n_{fold} = $' + str(n_folds))
	plt.xlabel(r"$\lambda$")

	plt.plot(alphass, CLF_MEAN, "b--")
	plt.fill_between(alphass, CLF_MEAN+CLF_WRMS, CLF_MEAN-CLF_WRMS, alpha=0.2)

	plt.savefig('CLF_MEAN_WRMS.png')
	plt.show()
	plt.close()








	n_folds = 5
	size 	= feature_matrix_int_lag.shape[0]
	idx_random = shuffle(np.arange(size))
	X 		= feature_matrix_int_lag[idx_random,:]
	y 		= rv_daily[index][idx_random]
	sw 		= 1/erv_daily[index][idx_random]**2

	lasso = Lasso(random_state=0, max_iter=10000)
	alphas = np.logspace(-5, 0, 100)

	tuned_parameters = [{"alpha": alphas}]

	clf = GridSearchCV(lasso, tuned_parameters, cv=n_folds, refit=True, return_train_score=True)
	clf.fit(X, y, sample_weight=sw)
	clf.best_estimator_.coef_



# from sklearn.linear_model import LassoCV
# from sklearn.datasets import make_regression
# X, y = make_regression(noise=4, random_state=0)
# reg = LassoCV(cv=5, random_state=0).fit(X, y, sample_weight=y)
# reg.score(X, y)

# reg.predict(X[:1,])


	# --- explore the res_wrms vs lambda
	N = 20
	res_wrms_array = np.zeros(N)
	for i in range(N):
		alpha = i/N/20
		y_hat3, w_std_all3, res_wrms, score, variance_matrix3 = mlr(feature_matrix_int_lag, target_vector=rv_daily[index], etarget_vector=erv_daily[index], alpha=alpha, feature_matrix2=feature_matrix[index])
		res_wrms_array[i] = res_wrms
	plt.plot(np.arange(20)/20/20, res_wrms_array, label=str(k_max))	
plt.legend(loc='best')
plt.savefig('rms_vs_lambda.pdf')
plt.show()
	# ---


	y_hat6, w_std_all6, res_wrms, score, variance_matrix6 = mlr(feature_matrix_int_lag, target_vector=rv_daily[index], etarget_vector=erv_daily[index], feature_matrix2=feature_matrix[index])
	imshow_matrix(variance_matrix6, file_name='fwhm_bis_multi_coef') 

	if 0: # test 
		for i in range(11):
			plt.plot(bjd_daily_lag, feature_matrix_int_lag[:,i*3], '.', alpha=0.5)
			plt.plot(bjd_daily_lag, feature_matrix_int_lag[:,i*3], '-', alpha=0.2)	
		plt.savefig('test.png')
		plt.show()

		plt.plot(bjd_daily_lag, y_hat3 - rv_daily[index], '.')
		plt.title('residual')
		plt.show()



	#---------------------------------------------------------------------------------#
	# Compare the performance of two models
	#---------------------------------------------------------------------------------#
	# start with a square Figure
	plt.rcParams.update({'font.size': 14})
	fig = plt.figure(figsize=(8, 8))

	# Add a gridspec with two rows and two columns and a ratio of 2 to 7 between
	# the size of the marginal axes and the main axes in both directions.
	# Also adjust the subplot parameters for a square plot.
	gs = fig.add_gridspec(2, 2,  width_ratios=(7, 2), height_ratios=(2, 7),
	                      left=0.1, right=0.9, bottom=0.1, top=0.9,
	                      wspace=0.05, hspace=0.05)

	ax = fig.add_subplot(gs[1, 0])
	ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
	ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

	# use the previously defined function
	res3 = y_hat3-rv_daily[index]
	res6 = y_hat6-rv_daily[index]

	scatter_hist(res3, res6, erv_daily[index], erv_daily[index], ax, ax_histx, ax_histy)

	plt.savefig('scatter_histogram.pdf')
	plt.show()

	_ = plt.hist(y_hat3-y_hat6, bins=20, color='black', alpha=0.5)
	plt.xlabel('Model 3 - Model 6 [m/s]')
	plt.savefig('model_comparison_histogram.pdf')
	plt.show()

	np.percentile(y_hat3-y_hat6, 50+34.1) - np.percentile(y_hat3-y_hat6, 50-34.1)

	# -----------------------
	# fwhm_bis with lags 
	# -----------------------
	for iday in range(5):
		day = iday+3
		
		bjd_daily = np.loadtxt('bjd_daily.txt')

		fwhm_bis = np.vstack((fwhm_daily, bis_daily)).T
		scaler = StandardScaler()
		fwhm_bis = scaler.fit_transform(fwhm_bis)
		efwhm_bis = np.vstack((efwhm_daily/np.std(fwhm_daily), ebis_daily/np.std(bis_daily))).T  

		x_fwhm_bis_int = np.zeros((len(bjd_int), 2))
		for i in range(2):
			f = interp1d(bjd_daily, fwhm_bis[:,i], fill_value='extrapolate')
			x_fwhm_bis_int[:,i] = f(bjd_int)

		x_fwhm_bis_int_lag = np.zeros((len(rv_daily_int[day:-day]), 2*(2*day+1)))
		for n in range(len(rv_daily_int[day:-day])): #new
			for i in range((2*day+1)):
				x_fwhm_bis_int_lag[n, (2*i):(2*i+2)] = x_fwhm_bis_int[n+i, 0:2]

		rv_daily_int_lag = rv_daily_int[day:-day]

		lasso = Lasso(alpha=alpha, max_iter=10000).fit(x_fwhm_bis_int_lag, rv_daily_int_lag, sample_weight=1/erv_daily_int[day:-day]**2)

		coeff_matrix = np.zeros((day*2+1, 2))
		for i in range(day*2+1):
			coeff_matrix[i, :] = lasso.coef_[(i*2):(i*2+2)]

		w_std = np.zeros(2)
		for i in range(2):
			_, w_std[i] = weighted_avg_and_std(x_fwhm_bis_int[day:-day,i], 1/erv_daily_int[day:-day]**2)

		variance_matrix = np.zeros(coeff_matrix.shape)
		for i in range(coeff_matrix.shape[0]):
			variance_matrix[i,:] = coeff_matrix[i,:]*w_std
		variance_matrix = variance_matrix**2
		variance_matrix = variance_matrix / variance_matrix.sum() * 100

		imshow_matrix(variance_matrix, file_name='fwhm_bis_coef') # ax.set_yticklabels(['FWHM', 'BIS'])
>>>>>>> 6d94324e434f9f91ee1b557918fc0182432b9c1d

	# -----------------------

	k_feature = 2
	short_variation = np.zeros((k_feature, len(bjd_daily)))
	long_variation 	= np.zeros((k_feature, len(bjd_daily)))

	for i in range(k_feature):
		_, short_variation[i,:], long_variation[i,:] = long_short_divide(
			x=bjd_daily, y=fwhm_bis[:,i], yerr=efwhm_bis[:,i], r=100)

	plot_all_but_corr(k_mode=4, t=bjd_daily,
		ind=feature_matrix.T, eind=np.hstack([efwhm_bis, efwhm_bis]).T, 
		ts_xlabel='BJD - 2400000 [d]', 
		height_ratio = 0.8,
		pe_xlabel='Period [days]',
		ind_yalbel=['S-FWHM', 'S-BIS', 'L-FWHM', 'L-BIS'],
		file_name='fwhm_bis_variation.pdf')
	plt.close()


# test the weighted averaged shift
rv_test = np.zeros(shift_spectrum.shape[0])
for i in range(shift_spectrum.shape[0]):
	rv_test[i] = sum(shift_function[i] * power_spectrum[i]) / np.sum(power_spectrum[i])

# for i in range(shift_spectrum.shape[0]):
# 	rv_test[i] = sum(shift_spectrum[i] * power_spectrum[i]) / np.sum(power_spectrum[i])

#----------------#
# RV time-series #
#----------------#
plt.rcParams.update({'font.size': 14})
alpha=0.5
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

fig, axes = plt.subplots(figsize=(15, 3))
plt.gcf().subplots_adjust(bottom=0.2)
# plt.errorbar(bjd, rv-np.mean(rv_daily), rv_err, c='black', marker='.', ls='none', alpha= 0.05, label='Unbinned RV')
plt.errorbar(bjd_daily, rv_daily-np.mean(rv_daily), erv_daily, c='purple', marker='o', ls='none', alpha= 0.5, label='Daily RV')
# plt.plot(bjd_int, rv_daily_int-np.mean(rv_daily), '.')
# plt.title('HARPS-N three years solar RV')
plt.legend()
plt.xlabel('BJD - 2400000 [d]')
plt.ylabel('RV [m/s]')
plt.savefig('rv_daily.png')
plt.show()



from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

def periodogram_overlap(x=bjd_daily, y=power_spectrum, dy=err_power_spectrum,
				plot_min_t=1, study_min_t=5, max_f=1, spp=100, title='',
				label = '',
				file_name='Periodogram.png'):


	if y.ndim == 1:
		N = 1
	else:
		N = y.shape[1]

	time_span = (max(x) - min(x))
	min_f   = 1/time_span

	plt.rcParams.update({'font.size': 12})
	plt.subplots(figsize=(15, 2))
	plt.gcf().subplots_adjust(bottom=0.25)
	plt.title(title)

	for i in range(N):
		plt.subplot(N,1,i+1)
		if i == 0:
			plt.title(title)

		if y.ndim == 1:
			frequency, power = LombScargle(x, y, dy).autopower(minimum_frequency=min_f,
															   maximum_frequency=max_f,
															   samples_per_peak=spp)
		else:
			frequency, power = LombScargle(x, y[:, i], dy[:, i]).autopower(minimum_frequency=min_f,
																		   maximum_frequency=max_f,
																		   samples_per_peak=spp)

		plot_x = 1/frequency
		idxx = plot_x>study_min_t
		height = max(power[plot_x>study_min_t])*0.5
		plt.plot(plot_x, power, ls='-', alpha=0.8, label=label)
		plt.legend()
		peaks, _ = find_peaks(power[idxx], height=height*0.5)
		plt.plot(plot_x[peaks], power[peaks], "o")

		for n in range(len(plot_x[peaks])):
			plt.text(plot_x[peaks][n], power[peaks][n], '%.1f' % plot_x[peaks][n], fontsize=12)

		plt.xlim([plot_min_t,time_span])
		plt.ylim([0,2.5*height])
		plt.xscale('log')
		plt.ylabel('Power')
		plt.xlabel('Period [d]')

	plt.savefig(file_name)
	plt.show()

periodogram_overlap(x=bjd, y=df["fwhm"], dy=df["fwhm_err"], title='', label='FWHM', file_name = 'fwhm_Periodogram.png') # height*0.6

periodogram_overlap(x=bjd, y=df["bis_span"], dy=df["bis_span_err"], title='', label='BIS SPAN', file_name = 'bis_Periodogram.png') # height*0.5

periodogram_overlap(x=bjd_daily, y=rv_daily, dy=erv_daily, file_name = 'rv_Periodogram.png') # height*0.5

########
# FWHM #
########
plt.rcParams.update({'font.size': 12})

fig, axes = plt.subplots(figsize=(15, 2.5))
plt.gcf().subplots_adjust(bottom=0.25)
plt.errorbar(bjd, df["fwhm"]/1000, df["fwhm_err"]/1000, marker='.', ls='none', alpha= 0.1, label='Unbinned')
plt.errorbar(bjd_daily, fwhm_daily/1000, efwhm_daily/1000, marker='o', ls='none', alpha= 0.5, label='Daily binned')
# plt.title('HARPS-N three years solar spectra FWHM')
plt.legend()
plt.xlabel('BJD - 2400000 [d]')
plt.ylabel('FWHM [km/s]')
plt.savefig('FWHM.png')
plt.show()

#######
# BIS #
#######
plt.rcParams.update({'font.size': 12})

fig, axes = plt.subplots(figsize=(15, 2.5))
plt.gcf().subplots_adjust(bottom=0.25)
plt.errorbar(bjd, df["bis_span"], df["fwhm_err"], marker='.', ls='none', alpha= 0.1, label='Unbinned')
plt.errorbar(bjd_daily, bis_daily, ebis_daily, marker='o', ls='none', alpha= 0.5, label='Daily binned')
# plt.title('HARPS-N three years solar spectra BIS SPAN')
plt.legend()
plt.xlabel('BJD - 2400000 [d]')
plt.ylabel('BIS SPAN [m/s]')
plt.savefig('BIS_SPAN.png')
plt.show()





fig, axes = plt.subplots(figsize=(20, 4))
plt.errorbar(bjd_daily, rv_raw_daily - np.mean(rv_raw_daily), erv_daily, marker='o', ls='none', alpha= 0.5, label='rv_raw_daily')
plt.errorbar(bjd_daily, RV_gauss-np.mean(RV_gauss), erv_daily, marker='o', ls='none', alpha= 0.5, label='RV_Gauss')
plt.title('HARPS-N 3 yrs solar RV')
plt.legend()
plt.xlabel('bjd')
plt.ylabel('rv [m/s]')
plt.savefig('rv_raw_daily.png')
plt.show()

# I do not understand why rv_test is not 0 as expected
if 0:
	clean_rv = rv_daily - rv_test

	plt.plot(rv_daily, clean_rv, '.')
	plt.show()

	plt.plot(bjd_daily, rv_daily-np.mean(rv_daily), '.', label='rv_daily')
	plt.plot(bjd_daily, 4*(rv_test-np.mean(rv_test)), '.', label='rv_test')
	plt.xlabel('bjd')
	plt.ylabel('rv [m/s]')
	plt.legend()
	# plt.savefig('rv_daily-rv_test.png')
	plt.show()

	plt.plot(bjd_daily, rv_daily - rv_test, '.', label='rv_daily - rv_test')
	plt.xlabel('bjd')
	plt.ylabel('rv [m/s]')
	plt.legend()
	# plt.savefig('rv_daily-rv_test.png')
	plt.show()



	np.std(rv_daily)
	np.std(rv_daily - rv_test)

	np.std(rv_test)

	from astropy.timeseries import LombScargle
	from scipy.signal import find_peaks

	time_span = (max(bjd_daily) - min(bjd_daily))
	min_f   = 1/time_span
	max_f   = 1
	spp     = 100  # spp=1000 will take a while; for quick results set spp = 10
	xc 		= 365/2

	frequency0, power0 = LombScargle(bjd_daily, rv_daily, erv_daily).autopower(minimum_frequency=min_f,
													   maximum_frequency=max_f,
													   samples_per_peak=spp)

	plot_x = 1 / frequency0
	plt.plot(plot_x, power0, ls='-', label='rv_daily', alpha=0.8)
	idxx = plot_x > 5
	height = max(power0[plot_x > 2]) * 0.5
	plt.axvline(x=xc, color='k', linestyle='--', alpha=0.75)
	peaks, _ = find_peaks(power0[idxx], height=height)
	plt.plot(plot_x[peaks], power0[peaks], "x")
	for n in range(len(plot_x[peaks])):
		plt.text(plot_x[peaks][n], power0[peaks][n], '%.1f' % plot_x[peaks][n], fontsize=12)

	frequency0, power0 = LombScargle(bjd_daily, rv_test, erv_daily).autopower(minimum_frequency=min_f,
													   maximum_frequency=max_f,
													   samples_per_peak=spp)
	plt.plot(plot_x, -power0, ls='-', label='rv_test', alpha=0.8)

	idxx = plot_x > 5
	height = max(power0[plot_x > 2]) * 0.5
	plt.axvline(x=xc, color='k', linestyle='--', alpha=0.75)
	plt.axvline(x=13.5, color='k', linestyle='--', alpha=0.75)
	plt.axvline(x=27, color='k', linestyle='--', alpha=0.75)
	peaks, _ = find_peaks(power0[idxx], height=height)
	plt.plot(plot_x[peaks], -power0[peaks], "x")
	for n in range(len(plot_x[peaks])):
		plt.text(plot_x[peaks][n], -power0[peaks][n], '%.1f' % plot_x[peaks][n], fontsize=12)
	plt.xlim([1, time_span])
	plt.ylim([-2 * height, 2 * height])
	plt.xscale('log')
	plt.ylabel('Power%d' % (i + 1))
	plt.legend()
	plt.savefig('FIESTA_periodogram_comparison.png')
	plt.show()


if 0:
	# compare the rvs
	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, rv_raw-np.mean(rv_raw), '.', alpha=0.2, label='HARPS-N')
	plt.plot(bjd_daily, RV_gauss-np.mean(RV_gauss), '.', alpha=0.5, label='Gaussian fit of CCF_daily_corrected')
	plt.plot(bjd_daily, rv_daily*1000-np.mean(rv_daily*1000), '.', alpha=0.5, label='rv_daily')
	plt.xlabel('bjd')
	plt.ylabel('rv [m/s]')
	plt.legend()
	# plt.savefig('rv_raw_comparison.png')
	plt.show()

	for n in range(632):
		plt.plot(V_grid, CCF_daily_corrected[:,n] - CCF_daily[:,n], alpha = 0.2)
	plt.savefig('CCF_diff.png')	
	plt.close()

	for n in range(632):
		plt.plot(V_grid, CCF_daily[:,n] - np.mean(CCF_daily, axis=1), alpha = 0.2)
	# plt.savefig('CCF_diff.png')	
	# plt.close()
	for n in range(632):
		plt.plot(V_grid, CCF_daily_corrected[:,n] - np.mean(CCF_daily_corrected, axis=1), alpha = 0.2)

	for n in range(632):
		plt.plot(V_grid, CCF_daily_corrected[:,n] - CCF_daily[:,0], alpha = 0.2)

# Haven't figured it out yet!!!


#==============================================================================
# Plots 
#==============================================================================
plt.rcParams.update({'font.size': 12})

def time_series(x=bjd_daily, y=power_spectrum, dy=err_power_spectrum, N=None,
				ylabel='k=',
				title='Time series',
				file_name='Time_series.png'):
	if N==None:
		N = y.shape[1]
	plt.subplots(figsize=(12, N))

	for i in range(N):
		ax = plt.subplot(N, 1, i+1)
		if i == 0:
			plt.title(title)
		plt.errorbar(x, y[:, i], dy[:, i], marker='.', ls='none', alpha=0.5, ms=5)
		plt.ylabel(ylabel+str(i+1))
		if i != N-1:
			ax.set_xticks([])
		else:
			plt.xlabel('BJD - 2400000 [d]')
	plt.savefig(file_name)
	plt.show()


def periodogram(x=bjd_daily, y=power_spectrum, dy=err_power_spectrum, N=None,
				plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=None,
				ylabel=None,
				title = 'Periodogram',
				file_name='Periodogram.png'):

	if N==None:
		N = y.shape[1]
	time_span = (max(x) - min(x))
	min_f   = 1/time_span

	plt.subplots(figsize=(12, N))

	for i in range(N):
		ax = plt.subplot(N,1,i+1)
		if i == 0:
			plt.title(title)

		frequency, power = LombScargle(x, y[:, i], dy[:, i]).autopower(minimum_frequency=min_f,
													   maximum_frequency=max_f,
													   samples_per_peak=spp)

		plot_x = 1/frequency
		idxx = plot_x>study_min_t
		height = max(power[plot_x>study_min_t])*0.5
		plt.plot(plot_x, power, ls='-', label=r'$\xi$'+str(i+1), alpha=0.8)
		peaks, _ = find_peaks(power[idxx], height=height)
		plt.plot(plot_x[peaks], power[peaks], "o")
		if xc != None:
			plt.axvline(x=xc, color='k', linestyle='--', alpha = 0.5)

		for n in range(len(plot_x[peaks])):
			plt.text(plot_x[peaks][n], power[peaks][n], '%.1f' % plot_x[peaks][n], fontsize=10)

		plt.xlim([plot_min_t,time_span])
		plt.ylim([0,2.5*height])
		plt.xscale('log')
		plt.ylabel(ylabel+str(i+1))

		if i != N-1:
			ax.set_xticks([])
		else:
			plt.xlabel('Period [d]')

	plt.savefig(file_name)
	plt.show()


def bgls_periodogram(x=bjd_daily, y=power_spectrum, dy=err_power_spectrum, N=None,
				plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=None,
				ylabel=None,
				title = 'Periodogram',
				file_name='Periodogram.png'):

	if N==None:
		N = y.shape[1]
	time_span = (max(x) - min(x))
	min_f   = 1/time_span

	plt.subplots(figsize=(12, N))

	for i in range(N):
		ax = plt.subplot(N,1,i+1)
		if i == 0:
			plt.title(title)

		periods, probs = bgls(x, y[:, i], dy[:, i], plow=2, phigh=500, ofac=100)
		# plt.plot(periods, probs, '-'); 
		# # plt.xlim([0,20]); 
		# plt.show()
		plot_x = periods
		power = np.log(probs)
		idxx = (plot_x>=study_min_t)
		height = max(power[plot_x>study_min_t])*0.5
		plt.plot(plot_x, power, ls='-', label=r'$\xi$'+str(i+1), alpha=0.8)
		peaks, _ = find_peaks(power[idxx], height=height)
		plt.plot(plot_x[peaks], power[peaks], "o")
		if xc != None:
			plt.axvline(x=xc, color='k', linestyle='--', alpha = 0.5)

		for n in range(len(plot_x[peaks])):
			plt.text(plot_x[peaks][n], power[peaks][n], '%.1f' % plot_x[peaks][n], fontsize=10)

		plt.xlim([plot_min_t,time_span])
		plt.ylim([0,2.5*height])
		plt.xscale('log')
		plt.ylabel(ylabel+str(i+1))

		if i != N-1:
			ax.set_xticks([])
		else:
			plt.xlabel('Period [d]')

	plt.savefig(file_name)
	plt.show()	


time_series(x=bjd_daily, y=power_spectrum, dy=err_power_spectrum, N=None,
				title='$A_k$ time series',
				file_name='FIESTA_amplitude_time_series.png')

periodogram(x=bjd_daily, y=power_spectrum, dy=err_power_spectrum, N=None,
			plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=365/2,
			ylabel='k=',
			title='$A_k$ Periodogram',
			file_name='FIESTA_amplitude_periodogram.png')

time_series(x=bjd_daily, y=shift_function-np.mean(shift_function, axis=0), dy=err_shift_spectrum, ylabel='k=', N=None,
				title=r'$\Delta RV_k$ time series [m/s]',
				file_name='FIESTA_shift_time_series.png')

periodogram(x=bjd_daily, y=shift_function, dy=err_shift_spectrum, N=None,
			plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=365/2,
			ylabel='k=',
			title=r'$\Delta RV_k$ Periodogram',
			file_name='FIESTA_shift_periodogram.png')





#----------------------------
# New time-series
#----------------------------
def periodogram(x, y, dy, N=None,
				plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=None,
				ylabel=None,
				title = 'Periodogram',
				file_name='Periodogram.png'):
	
	from scipy.signal import find_peaks
	from astropy.timeseries import LombScargle

	if N==None:
		N = y.shape[1]
	time_span = (max(x) - min(x))
	min_f   = 1/time_span

	plt.subplots(figsize=(10, N))

	for i in range(N):
		ax = plt.subplot(N,1,i+1)
		if i == 0:
			plt.title(title)

		frequency, power = LombScargle(x, y[:, i], dy[:, i]).autopower(minimum_frequency=min_f,
													   maximum_frequency=max_f,
													   samples_per_peak=spp)

		plot_x = 1/frequency
		idxx = (plot_x>plot_min_t) & (plot_x<time_span/2)
		height = max(power[idxx])*0.4
		plt.plot(plot_x[idxx], power[idxx], 'k-', label=r'$\xi$'+str(i+1), alpha=0.5)
		peaks, _ = find_peaks(power[idxx], height=height)
		plt.plot(plot_x[idxx][peaks], power[idxx][peaks], "ro")
		if xc != None:
			plt.axvline(x=xc, color='k', linestyle='--', alpha = 0.5)

		for n in range(len(plot_x[idxx][peaks])):
			plt.text(plot_x[idxx][peaks][n], power[idxx][peaks][n], '%.1f' % plot_x[idxx][peaks][n], fontsize=10)

		plt.xlim([plot_min_t,time_span/2])
		plt.ylim([0, 3*height])
		plt.xscale('log')
		if i==0:
			plt.ylabel('HARPS-N')
		if i>0:
			plt.ylabel(ylabel+str(i))

		if i != N-1:
			ax.set_xticks([])
		else:
			plt.xlabel('Period [days]')

	plt.savefig(file_name)
	plt.show()


periodogram(x=bjd_daily, y=np.vstack((rv_daily, power_spectrum)).T, dy=np.vstack((erv_daily, err_power_spectrum)).T,
			plot_min_t=2, study_min_t=5, max_f=1, spp=100,
			ylabel='k=',
			file_name='FIESTA_amplitude_periodogram.png')

periodogram(x=bjd_daily, y=np.vstack((rv_daily, shift_function)).T, dy=np.vstack((erv_daily, err_shift_spectrum)).T,
			plot_min_t=2, study_min_t=5, max_f=1, spp=100,
			ylabel='k=',
			file_name='FIESTA_shift_periodogram.png')




#---------------------------------
# New time-series and periodogram 
#---------------------------------


names = ['sepal length [cm]', 'sepal width [cm]',
         'petal length [cm]', 'petal width [cm]']

fig, axes = scatterplotmatrix(np.vstack((rv_daily, shift_function)).T, figsize=(20, 20), alpha=0.5)
# fig, axes = scatterplotmatrix(X[y==1], fig_axes=(fig, axes), alpha=0.5)
# fig, axes = scatterplotmatrix(X[y==2], fig_axes=(fig, axes), alpha=0.5, names=names)

# np.vstack((rv_daily, shift_function))

plt.tight_layout()
plt.show()














# CCF_daily = CCF_daily[:, 0::5]
# eCCF_daily = eCCF_daily[:, 0::5]

# bjd_daily = bjd_daily[0::5]
# rv_daily = rv_daily[0::5]
# rv_raw_daily = rv_raw_daily[0::5]
# erv_daily = erv_daily[0::5]

# power_spectrum 		= power_spectrum[0::5, :]
# err_power_spectrum 	= err_power_spectrum[0::5, :]
# shift_spectrum 		= shift_spectrum[0::5, :]
# err_shift_spectrum 	= err_shift_spectrum[0::5, :]
# shift_function 		= shift_function[0::5]


#----------------------------------
# PCA 
#----------------------------------
if 0: # old

	X = shift_function

	from sklearn.preprocessing import StandardScaler
	X = StandardScaler().fit_transform(X)

	# Compute the Eigenvectors and Eigenvalues
	covariance_matrix = np.cov(X.T)
	eigen_values, eigen_vectors = np.linalg.eig(covariance_matrix)

	# Singular Value Decomposition (SVD)
	eigen_vec_svd, _, _= np.linalg.svd(X.T)

	variance_explained = [(i / sum(eigen_values)) * 100 for i in eigen_values]

	projection_matrix = (eigen_vectors.T)[:6].T
	# array = np.array([0,1,2,3,5,6,10,11,12,13])
	# projection_matrix = (eigen_vectors.T)[array].T

	X_pca = X.dot(projection_matrix)

	fig, axes = plt.subplots(figsize=(12, 10))
	for i in range(6):
		ax = plt.subplot(6,1,i+1)
		plt.plot(bjd_daily, X_pca[:, i], '.', alpha=0.5)
		plt.ylabel('PCA%d' %(i+1))
		if i != 6-1:
			ax.set_xticks([])
		else:
			plt.xlabel('date_bjd')		
	plt.savefig('PCA.png')
	plt.show()
	plt.close()


	# Method 2 
	n_pca = 6
	from sklearn.decomposition import PCA
	pca = PCA(n_components=n_pca)
	x_pca = pca.fit_transform(X)

	fig, axes = plt.subplots(figsize=(12, 18))
	for i in range(n_pca):
		ax = plt.subplot(n_pca,1,i+1)
		plt.plot(bjd_daily, x_pca[:, i], '.', alpha=0.8)
		# plt.plot(bjd_daily[19:613-1], moving_average(x_pca[:, i], 40), '-', alpha=1)
		plt.ylabel('PCA%d' %(i+1))
		if i != n_pca-1:
			ax.set_xticks([])
		else:
			plt.xlabel('date_bjd')		
	plt.savefig('skl_PCA.png')
	plt.show()
	plt.close()

	print(pca.explained_variance_ratio_)

	cumulative_variance_explained = np.cumsum(pca.explained_variance_ratio_) * 100

	plt.plot(np.arange(n_pca)+1, cumulative_variance_explained, '.')
	plt.xlabel("Number of components")
	plt.ylabel("Cumulative explained variance")
	plt.title("Explained variance vs Number of components")
	plt.savefig('cumulative_variance_explained.png')
	plt.show()


#----------------------------------
# weighted pca
#----------------------------------

# X.shape 
# n_samples x n_features

'''
- Ludovic -- Weighted principal component analysis: a weighted covariance eigendecomposition approach https://doi.org/10.1093/mnras/stu2219
- wPCA python package
		https://github.com/jakevdp/wpca
'''
if 0: # using the wPCA package

	from wpca import PCA, WPCA, EMPCA

	ThisPCA = WPCA
	X = shift_function
	weights = 1. / err_shift_spectrum**2

	if weights is None:
	    kwds = {}
	else:
	    kwds = {'weights': weights}

	# Compute the PCA vectors & variance
	pca = ThisPCA()

	pca.fit(X, **kwds)
	pca_score = pca.transform(X, **kwds)

	# shouldn't be doing the wpca package for the scores because 
	# the wpca package does not consider the non-diagonal terms 
	for i in range(err_shift_spectrum.shape[0]):
		w_matrix = np.diag(1. / err_shift_spectrum[i,0:10]**2)
		err_score = np.linalg.inv(pca.components_.T @ w_matrix @ pca.components_)

	# PCA variance
	pca.explained_variance_ratio_[0:10]

	print(sum(pca.explained_variance_ratio_[0:9]))


# for now, only consider the diagonal terms 
# X 		= power_spectrum
# weights = 1 / err_power_spectrum  # may need to include the Fourier power later
X 		= shift_function
weights = 1 / err_shift_spectrum # may need to include the Fourier power later

np.savetxt('X.txt', X_new)
np.savetxt('X_err.txt', err_shift_spectrum)
np.savetxt('W.txt', weights)


# weights[:, 6:11] = weights[:, 6:11] * 10


time_series(x=bjd_daily, y=my_pca_score, dy=my_err_pca_score, N=6,
			ylabel='PCA',
			title=r'$\Delta RV_k$ PCA Scores Time Series',
			file_name='FIESTA_RV_my_pca_score_time_series.png')




k_mode 	= 8
alpha1, alpha2 = [0.5,0.2]
widths 	= [8,1]
heights = [1]*(k_mode+1)
gs_kw 	= dict(width_ratios=widths, height_ratios=heights)
plt.rcParams.update({'font.size': 10})
fig6, f6_axes = plt.subplots(figsize=(10, k_mode+1), ncols=2, nrows=k_mode+1, constrained_layout=True,
                             gridspec_kw=gs_kw)
for r, row in enumerate(f6_axes):
	for c, ax in enumerate(row):		
		if c==0:
			if r==0:
				ax.errorbar(bjd_daily, rv_daily-np.mean(rv_daily), erv_daily, marker='.', color='black', ls='none', alpha=alpha1)
				ax.set_title('Time-series')
				ax.set_ylabel('$RV_{HARPS}$')
			else:				
				ax.errorbar(bjd_daily, my_pca_score.T[r-1,:], my_err_pca_score.T[r-1,:],  marker='.', color='black', ls='none', alpha=alpha1)
				ax.set_ylabel('PCA%d' %(r))
			if r!=k_mode:
				ax.set_xticks([])
			else:
				ax.set_xlabel('BJD - 2400000 [d]')
		if c==1:
			if r==0:
				reg = LinearRegression().fit(rv.reshape(-1, 1), rv.reshape(-1, 1))
				score = reg.score(rv.reshape(-1, 1), rv.reshape(-1, 1))
				ax.set_title('score = {:.2f}'.format(score))
				ax.plot(rv_daily-np.mean(rv_daily), rv_daily-np.mean(rv_daily), 'k.', alpha = alpha2)				
			if r>0:
				reg = LinearRegression().fit(rv_daily.reshape(-1, 1), my_pca_score.T[r-1,:].reshape(-1, 1))
				score = reg.score(rv_daily.reshape(-1, 1), my_pca_score.T[r-1,:].reshape(-1, 1))
				ax.set_title('score = {:.2f}'.format(score))
				ax.plot(rv_daily-np.mean(rv_daily), my_pca_score.T[r-1,:], 'k.', alpha = alpha2)
			if r!=k_mode:
				ax.set_xticks([])
			else:
				ax.set_xlabel('$RV_{HARPS}$')
			ax.yaxis.tick_right()
plt.savefig('time-series_and_PCA_shift_correlation.png')
# plt.savefig('time-series_and_PCA_A_correlation.png')			
plt.show()

periodogram(x=bjd_daily, y=np.vstack((rv_daily, my_pca_score.T)).T, dy=np.vstack((erv_daily, my_err_pca_score.T)).T, N=7,
			plot_min_t=2, study_min_t=5, max_f=1, spp=100, xc=365/2,
			ylabel='PCA',
			title=r'$\Delta RV_k$ PCA Scores Periodogram',
			file_name='FIESTA_RV_PCA_periodogram.png')

'''
periodogram(x=bjd_daily, y=my_pca_score, dy=my_err_pca_score, N=6,
			plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=365/2,
			ylabel='PCA',
			title=r'$\Delta RV_k$ PCA Scores Periodogram',
			file_name='FIESTA_RV_PCA_periodogram.png')
'''


# time_series(x=bjd_daily, y=my_pca_score, dy=my_err_pca_score, N =6,
# 			ylabel='PCA',
# 			title=r'$A_k$ PCA Scores Time Series',
# 			file_name='FIESTA_amplitudes_my_pca_score_time_series.png')

periodogram(x=bjd_daily, y=np.vstack((rv_daily, my_pca_score.T)).T, dy=np.vstack((erv_daily, my_err_pca_score.T)).T, N=7,
			plot_min_t=2, study_min_t=5, max_f=1, spp=100, xc=365/2,
			ylabel='PCA',
			title='$A_k$ PCA Scores Periodogram',
			file_name='FIESTA_amplitude_PCA_periodogram.png')
'''
periodogram(x=bjd_daily, y=my_pca_score, dy=my_err_pca_score, N=6,
			plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=365/2,
			ylabel='PCA',
			title='$A_k$ PCA Scores Periodogram',
			file_name='FIESTA_amplitude_PCA_periodogram.png')
'''

my_pca_score1=my_pca_score
my_err_pca_score1 = my_err_pca_score
my_pca_score2=my_pca_score
my_err_pca_score2 = my_err_pca_score

plt.rcParams.update({'font.size': 12})
N = my_pca_score.shape[1]
plt.subplots(figsize=(12, N))

for i in range(N):
	ax = plt.subplot(N, 1, i+1)
	if i == 0:
		plt.title('PCA scores')
	plt.plot(bjd_daily, pca_score[:,i], '.', alpha=0.2)
	# plt.plot(bjd_daily, pca_score3[:,i], '.', alpha=0.2)
	plt.ylabel(str(i+1))
	if i != N-1:
		ax.set_xticks([])
	else:
		plt.xlabel('date_bjd')
# plt.savefig('pca_score_comparison4.png')
plt.show()



pca_score = np.loadtxt('pca_score.txt', delimiter=',')
pca_score2 = np.loadtxt('pca_score2.txt', delimiter=',')
pca_score3 = np.loadtxt('pca_score3.txt', delimiter=',')

time_series(x=bjd_daily, y=my_pca_score, dy=my_err_pca_score, N = 6,
			ylabel='PCA',
			title='PCA Scores Time Series',
			file_name='FIESTA_my_pca_score_time_series.png')

time_series(x=bjd_daily, y=C.T, dy=err_C.T*0, N=n_pca,
			ylabel='PCA',
			title='C time series',
			file_name='FIESTA_C_time_series.png')

time_series(x=np.arange(11), y=P, dy=np.zeros(P.shape),
			ylabel='PCA',
			title='my_pca_score time series',
			file_name='FIESTA_P_time_series.png')
# my_pca_score is consistent with the C calculated below

# now not working with the WPCA(n_components=n_pca).fit(X_new, **kwds)


# determine how many pca scores are needed



for i in range(n_pca):
	np.savetxt('C' + str(i+1) + '.txt', C[i,:])
	np.savetxt('err_C' + str(i+1) + '.txt', err_C[i, :])


if 0:

	# Compute the PCA vectors & variance
	pca = WPCA(n_components=12)

	pca.fit(X_new, **kwds)
	pca_score = pca.transform(X_new, **kwds)
	print(pca_score.shape) # (632, 17)
	P = pca.components_ #(12, 17)
	print(P.shape)


#----------------------------------
# run mcmc on X
#----------------------------------
P, pca_score = WPCA_results(X=shift_function, weights=1/err_shift_spectrum)
P2, pca_score2 = WPCA_results2(X=shift_function, weights=1/err_shift_spectrum)

time_series(x=bjd_daily, y=pca_score, dy=err_C.T, N=n_pca,
			ylabel='PCA',
			title='pca_score time series',
			file_name='FIESTA_pca_score_time_series.png')

periodogram(x=bjd_daily, y=my_pca_score, dy=my_err_pca_score, N=n_pca,
			plot_min_t=1, study_min_t=5, max_f=1, spp=100,
			title='PCA scores Periodogram',
			file_name='FIESTA_pca_score_periodogram.png')


time_series(x=bjd_daily, y=pca_score2, dy=err_C.T, N=n_pca,
			ylabel='PCA',
			title='pca_score2 time series',
			file_name='FIESTA_pca_score2_time_series.png')



# get the PCA components and scores
def WPCA_results(X, weights=None, ncomp=None):

	if weights is None:
	    kwds = {}
	else:
	    kwds = {'weights': weights}	

	if ncomp is None:
		ncomp = X.shape[1]

	X_new = np.zeros(X.shape)

	for i in range(X.shape[1]):
		X_new[:,i] = X[:,i] - np.average(X[:,i], weights=weights[:,i])

	pca = WPCA(n_components = ncomp)
	pca.fit(X_new, **kwds)
	pca_score = pca.transform(X_new, **kwds)
	P = pca.components_	

	return P, pca_score


# def WPCA_results2(X, weights=None, ncomp=None):
# 
# 	if weights is None:
# 	    kwds = {}
# 	else:
# 	    kwds = {'weights': weights}
# 
# 	if ncomp is None:
# 		ncomp = X.shape[1]
# 
# 	X_new = np.zeros(X.shape)
# 
# 	for i in range(X.shape[1]):
# 		X_new[:,i] = X[:,i] - np.average(X[:,i], weights=weights[:,i])
# 
# 	pca_score = pca.transform(X_new, **kwds)
# 	P = pca.components_
# 
# 	return P, pca_score



pca_score_array = np.zeros((X.shape[0], 12, 1000))

for k in range(1000):

	X_noise = np.random.normal(X, err_shift_spectrum)
	_, pca_score_array[:,:,k] = WPCA_results2(X_noise, weights=weights, ncomp=12)

err_score_mcmc = np.std(pca_score_array, axis=2)


for k in range(pca_score_array.shape[1]):
	plt.hist(pca_score_array[0,k,:], bins = 20)
	plt.savefig('histogram_pca_' + str(k+1) + '.png')
	plt.close()

# compare the wpca results and the one from the paper

# C(:,i) = inv(P'*w*P) * (P'*w*X(:,i))

X_wpca = X_new.T
X_wpca.shape #(17, 632)

C 		= np.zeros((12, X.shape[0]))
err_C 	= np.zeros(C.shape)

W = weights.T

P = pca.components_.T
P.shape # (17, 12)

# W = W[:, 0:10]

for i in range(X_wpca.shape[1]):

	w = np.diag(W[:,i])**2
	C[:,i] 	= np.linalg.inv(P.T @ w @ P) @ (P.T @ w @ X_wpca[:,i])
	
	# the error is discribed by a covariance matrix of C
	Cov_C 	= np.linalg.inv(P.T @ w @ P) 

	# inv(P'*w*P) * P' * w * cov(X(:,i)) * w * P * inv(P'*w*P)
	# Cov_C = np.linalg.inv(P.T @ w @ P) @ P.T @ w @ np.diag(err_shift_spectrum[i,:])**2 @ w @ P @ np.linalg.inv(P.T @ w @ P)

	diag_C  = np.diag(Cov_C)
	err_C[:, i] = diag_C**0.5


for i in range(5):
	np.savetxt('C' + str(i+1) + '.txt', C[i,:])
	np.savetxt('err_C' + str(i+1) + '.txt', err_C[i, :])

#----------------------------------
# plot the pca scores
#----------------------------------

plt.rcParams.update({'font.size': 16})
fig, axes = plt.subplots(figsize=(18, n_pca*1.5))
for i in range(n_pca):
	ax = plt.subplot(n_pca,1,i+1)
	# 
	# plt.plot(bjd_daily, C[i,:], '.', alpha=0.3)
	# plt.plot(bjd_daily, pca_score[:, i], '.', alpha=0.3)
	plt.errorbar(bjd_daily, C[i, :], err_C[i, :], marker='.', ls='none', alpha= 0.5)
	# plt.errorbar(bjd_daily, pca_score[:, i], err_score_mcmc[:, i], marker='.', ls='none', alpha= 0.2)
	# plt.plot(bjd_daily[19:613-1], moving_average(x_pca[:, i], 40), '-', alpha=1)
	# plt.plot(bjd_daily, err_C[i, :], '.', alpha=0.3)
	# plt.plot(bjd_daily, err_score_mcmc[:, i], '.', alpha=0.3)
	plt.ylabel('PCA%d' %(i+1))
	if i != n_pca-1:
		ax.set_xticks([])
	else:
		plt.xlabel('date_bjd')		
# plt.savefig('wPCA.png')
# plt.savefig('wPCA_err.png')
plt.savefig('C.png')
# plt.savefig('C_err.png')
# plt.savefig('Compare_score.png')
# plt.savefig('Compare_errors2.png')
plt.show()

# or
time_series(x=bjd_daily, y=C.T, dy=err_C.T, N=n_pca,
			ylabel='PCA',
			title=r'$\Delta RV$ PCA time series',
			file_name='FIESTA_shift_PCA_time_series.png')


# plt.plot(np.mean(eCCF_daily, axis=0), err_score_mcmc[:, 0], '.')
# plt.show()

#----------------------------------
# pca scores periodogram
#----------------------------------
periodogram(x=bjd_daily, y=C.T, dy=err_C.T, N=n_pca,
			plot_min_t=1, study_min_t=5, max_f=1, spp=100, xc=398.98/2,
			title=r'$\Delta RV$ PCA Periodogram',
			file_name='FIESTA_shift_PCA_periodogram.png')


if 1:
	#----------------------------------
	# *principal component* periodogram
	#----------------------------------
	from astropy.timeseries import LombScargle
	from scipy.signal import find_peaks
	fig, axes = plt.subplots(n_pca, 1)

	time_span = (max(bjd_daily) - min(bjd_daily))
	min_f   = 1/time_span
	max_f   = 1
	spp     = 100  # spp=1000 will take a while; for quick results set spp = 10
	xc 		= 365/2

	fig, axes = plt.subplots(figsize=(8, 12))
	plt.title('Periodogram')
	for i in range(n_pca):
		ax = plt.subplot(n_pca,1,i+1)
		frequency0, power0 = LombScargle(bjd_daily, x_pca[:, i]).autopower(minimum_frequency=min_f,
	                                                   maximum_frequency=max_f,
	                                                   samples_per_peak=spp)

		plot_x = 1/frequency0
		plt.plot(plot_x, power0, ls='-', label=r'$\xi$'+str(i+1), alpha=0.8)
		plt.axvline(x=xc, color='k', linestyle='--', alpha = 0.75)
		peaks, _ = find_peaks(power0, height=max(power0[plot_x>2])*0.5)
		plt.plot(plot_x[peaks], power0[peaks], "x")
		for n in range(len(plot_x[peaks])):
			plt.text(plot_x[peaks][n], power0[peaks][n], '%.1f' % plot_x[peaks][n], fontsize=12)
		plt.xlim([1,time_span])
		plt.ylim([0,max(power0)])
		plt.xscale('log')
		plt.ylabel('PCA%d' %(i+1))
		# ax.set_xticks([])
	plt.savefig('x_pca_periodogram.png')


	#----------------------------------
	# RV periodogram
	#----------------------------------

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.title('Periodogram')
	frequency0, power0 = LombScargle(bjd, rv, rv_err).autopower(minimum_frequency=min_f,
	                                               maximum_frequency=max_f,
	                                               samples_per_peak=spp)
	frequency1, power1 = LombScargle(bjd_daily, rv_daily).autopower(minimum_frequency=min_f,
	                                               maximum_frequency=max_f,
	                                               samples_per_peak=spp)

	plot_x = 1/frequency0
	plt.plot(plot_x, power0, ls='-', alpha=0.5, label='raw')
	plt.plot(1/frequency1, power1, ls='-', alpha=0.5, label='daily')
	plt.axvline(x=xc, color='k', linestyle='--', alpha = 0.75)
	peaks, _ = find_peaks(power0, height=max(power0[plot_x>2])*0.15)
	plt.plot(plot_x[peaks], power0[peaks], "x")
	for n in range(len(plot_x[peaks])):
		plt.text(plot_x[peaks][n], power0[peaks][n], '%.1f' % plot_x[peaks][n], fontsize=12)
	plt.xlim([1,time_span])
	plt.ylim([0,max(power0)])
	plt.legend()
	plt.xscale('log')
	plt.ylabel('power')
	plt.savefig('rv_periodogram.png')



	XX = sigma_by_order
	XX = StandardScaler().fit_transform(XX)

	# Compute the Eigenvectors and Eigenvalues
	covariance_matrix = np.cov(XX.T)
	eigen_values, eigen_vectors = np.linalg.eig(covariance_matrix)

	# Singular Value Decomposition (SVD)
	# eigen_vec_svd, _, _= np.linalg.svd(XX.T)

	variance_explained = [(i / sum(eigen_values)) * 100 for i in eigen_values]

	projection_matrix = (eigen_vectors.T)[:2].T

	X_pca = X.dot(projection_matrix)

	fig, axes = plt.subplots(figsize=(12, 4))
	for i in range(2):
		ax = plt.subplot(2,1,i+1)
		plt.plot(bjd_daily, X_pca[:, i], '.', alpha=0.5)
		plt.ylabel('PCA%d' %(i+1))
		if i != 2-1:
			ax.set_xticks([])
		else:
			plt.xlabel('date_bjd')		
	plt.savefig('PCA.png')
	plt.show()
	plt.close()


	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, rv, '.', alpha=0.2, label='rv_sid')
	plt.plot(bjd_daily, rv_daily, '.', alpha=0.5, label='rv_sid_daily')
	plt.plot(bjd_daily, rv_daily)
	plt.legend()
	plt.ylabel('rv [m/s]')
	plt.savefig('rv_sid.png')	
	plt.close()
	# plt.show()

	fig, axes = plt.subplots(figsize=(12, 3))
	plt.plot(bjd, shift_spectrum[:, 0] - df['rv_raw'], '.', alpha=0.5)
	# plt.show()

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, df['rhk'], '.', alpha=0.2)
	plt.ylabel('rhk')
	plt.savefig('rhk.png')	
	# plt.show()
	plt.close()

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd, df['bis_span'], '.', alpha=0.2)
	plt.ylabel('bis_span')
	plt.savefig('bis_span.png')	
	# plt.show()

	# plt.plot(df['bis_span'], df['rv'], '.', alpha=0.2)
	# plt.show()

	fig, axes = plt.subplots(figsize=(12, 4))
	plt.plot(bjd[idx1], rv[idx1], '.', alpha=0.2, label='rv_sid')
	plt.ylabel('rv [m/s]')
	plt.savefig('rv_sid_s2.png')	
	plt.close()


	# Take a detailed look at 2015-07-29: file[9] to file[76]
	# observing from 9:01 to 15:24
	# group by hour 

	idx 				= (-10<V_grid) & (V_grid<10)
	CCF_2015_07_29 		= np.zeros((49, 76-9+1))
	CCF_2015_07_29_nor 	= np.zeros((49, 76-9+1))
	CCF_2015_07_29_daily= np.zeros(49)
	# date0 				= df['filename'][i][8:18]

	for i in range(9,76+1):
		hdulist     = fits.open(valid_ccf[i])
		data        = hdulist[1].data		
		CCF_2015_07_29[:, i-9] 		= data[69,:]
		CCF_2015_07_29_nor[:, i-9] 	= 1 - data[69,:] / np.mean(data[69,~idx])
		CCF_2015_07_29_daily 		+= data[69,:]
	CCF_2015_07_29_daily = 1 - CCF_2015_07_29_daily / np.mean(CCF_2015_07_29_daily[~idx])

	plt.plot(V_grid, CCF_2015_07_29_daily, label='CCF_2015_07_29_daily')
	plt.savefig('CCF_2015_07_29_daily.png')
	plt.xlabel('V_grid')
	plt.ylabel('CCF')
	plt.legend()
	plt.close()

	for i in range(68):
		plt.plot(V_grid, CCF_2015_07_29_nor[:,i] - CCF_2015_07_29_daily)
	plt.xlabel('V_grid')
	plt.ylabel('CCF')
	# plt.legend()
	plt.savefig('CCF_2015_07_29.png')
	plt.close()

	sigma_2015_07_29 = np.zeros(68)
	for i in range(68):
		popt, pcov 	= curve_fit(gaussian, V_grid, CCF_2015_07_29_nor[:,i], p0=[0.5, (max(V_grid)+min(V_grid))/2, 1, 0])
		sigma_2015_07_29[i] = popt[2]

	popt, pcov 	= curve_fit(gaussian, V_grid, CCF_2015_07_29_daily, p0=[0.5, (max(V_grid)+min(V_grid))/2, 1, 0])
	sigma_2015_07_29_daily = popt[2]



# compare the hyperparameters

hp = np.zeros((6,19))
dir = '/Users/az/Documents/GitHub/GLOM_RV_Example/hyperparameters/'
for i in range(6):
	hp[i,:] = np.loadtxt(dir+'HARPS_N_fit1_total_hyperparameters_'+str(i)+'.txt')

for i in range(19):
	plt.plot(np.arange(6)+1, hp[:,i], '-')
	plt.title('hyperparameter_'+str(i+1))
	plt.savefig(dir+'hyperparameter_'+str(i+1)+'.png')
	plt.close()


