'''
	Test how mlr performs for the first trunks of RV_HARPS
'''
bjd_daily 		= np.loadtxt('bjd_daily.txt')
idx_bjd 		= bjd_daily<58100
bjd_daily_part1 = bjd_daily[idx_bjd]

t1_min = min(bjd_daily_part1) 	#= 57233.05440362564
t1_max = max(bjd_daily_part1)	#= 58068.0739790099

bjd_daily_lag1 	= bjd_daily_part1[(bjd_daily_part1>=t1_min+day-0.3) & (bjd_daily_part1<=t1_max-day+0.3)]

bjd_daily_lag 	= bjd_daily_lag1

feature_matrix_int_lag 	= np.zeros((len(bjd_daily_lag), k_feature2*(day+1)))

for i in range((2*day+1)):	
	for j in range(k_feature):
		f = interp1d(bjd_daily, feature_matrix[:,j], fill_value='extrapolate')
		feature_matrix_int_lag[:, k_feature*i+j] = f(bjd_daily_lag-day+i)

		# cs = CubicSpline(bjd_daily, feature_matrix[:,j], extrapolate=True)
		# feature_matrix_int_lag[:, k_feature*i+j] = cs(bjd_daily_lag-day+i)

index = ((bjd_daily>=t1_min+day-0.3) & (bjd_daily<=t1_max-day+0.3))

feature_matrix_int_lag[:,k_feature*(2*day+1):] = feature_matrix[index, k_feature:]

y_hat3, w_std_all3, w_rms, score, variance_matrix3 = mlr(feature_matrix_int_lag, target_vector=rv_daily[index], etarget_vector=erv_daily[index], feature_matrix2=feature_matrix[index])
imshow_matrix(variance_matrix3, file_name='fiesta_multi_coef')

#---------------------------------------------------------------------------------#
# Multiple Regression Model with multidays for short-term variabilities 
#---------------------------------------------------------------------------------#
'''
	Interpolate RV_HARPS and indicators to the same evenly spaced V_grid
	Not recommended because it changes the original data
'''
bjd_daily 		= np.loadtxt('bjd_daily.txt')
idx_bjd 		= bjd_daily<58100
bjd_daily_part 	= bjd_daily[idx_bjd]
decimal  		= np.median(bjd_daily_part - [int(bjd_daily_part[i]) for i in np.arange(len(bjd_daily_part))])
bjd_int 		= np.arange(int(min(bjd_daily_part)), int(max(bjd_daily_part))) + decimal
f 				= interp1d(bjd_daily, rv_daily, fill_value='extrapolate')
rv_daily_int 	= f(bjd_int)
f 				= interp1d(bjd_daily, erv_daily, fill_value='extrapolate')
erv_daily_int	= f(bjd_int)

for iday in range(5):
	day 		= iday+3

	k_feature2 	= feature_matrix.shape[1]
	k_feature 	= int(k_feature2/2)
	
	feature_matrix_int = np.zeros((len(bjd_int), k_feature2))
	for i in range(k_feature2):
		f = interp1d(bjd_daily, feature_matrix[:,i],fill_value='extrapolate')
		feature_matrix_int[:,i] = f(bjd_int)
	
	rv_daily_int_lag = rv_daily_int[day:-day]

	feature_matrix_int_lag = np.zeros((len(rv_daily_int[day:-day]), k_feature2*(day+1)))

	for n in range(len(rv_daily_int[day:-day])):
		for i in range((2*day+1)):			
			feature_matrix_int_lag[n, (k_feature*i):(k_feature*i+k_feature)] = feature_matrix_int[n+i, 0:k_feature] #(835,6)
	feature_matrix_int_lag[:,k_feature*(2*day+1):] = feature_matrix_int[day:-day, k_feature:]

	y_hat3, w_std_all3, w_rms, score, variance_matrix3 = mlr(feature_matrix_int_lag, target_vector=rv_daily_int_lag, etarget_vector=erv_daily_int[day:-day], feature_matrix2=feature_matrix_int[day:-day, :])
	imshow_matrix(variance_matrix3, file_name='fiesta_multi_coef')

	# FWHM-BIS L and S
	y_hat6, w_std_all6, w_rms, score, variance_matrix6 = mlr(feature_matrix_int_lag, target_vector=rv_daily_int_lag, etarget_vector=erv_daily_int[day:-day], feature_matrix2=feature_matrix)
	imshow_matrix(variance_matrix6, file_name='fwhm_bis_multi_coef') # not working? 
		