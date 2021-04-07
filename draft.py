
for N_FIESTA in np.arange(5,12):

	print('-------------------------------')
	print('-----------N={:d}-----------------'.format(N_FIESTA))
	print('-------------------------------')

	shift_spectrum, err_shift_spectrum, power_spectrum, err_power_spectrum, RV_gauss = FIESTA(V_grid, CCF_daily, eCCF_daily, out=N_FIESTA)
	shift_spectrum 		= shift_spectrum * 1000
	err_shift_spectrum 	= err_shift_spectrum * 1000


	shift_function = np.zeros(shift_spectrum.shape)
	for i in range(shift_spectrum.shape[1]):
		shift_function[:,i] = shift_spectrum[:,i] - rv_daily

	N_FIESTA_freq 	= shift_spectrum.shape[1]

	# WPCA 
	from wpca import PCA, WPCA, EMPCA
	X 		= shift_function
	weights = 1 / err_shift_spectrum  # may need to include the Fourier power later
	# X 		= power_spectrum
	# weights = 1 / err_power_spectrum  # may need to include the Fourier power later

	mean 	= np.zeros(X.shape[1])
	std 	= np.zeros(X.shape[1])
	for i in range(X.shape[1]):
		mean[i] 		= np.average(X[:,i], weights=weights[:,i])
		std[i] 			= np.average((X[:,i]-mean[i])**2, weights=weights[:,i])**0.5
		
	
	if weights is None:
	    kwds = {}
	else:
	    kwds = {'weights': weights}

	print(X.shape) # (632, 17)

	# subtract the weighted mean for each measurement type (dimension) of X
	X_new = np.zeros(X.shape)

	if 0: # or?
		from sklearn.preprocessing import StandardScaler	
		X_new = StandardScaler().fit_transform(X)

	for i in range(X.shape[1]):
		X_new[:,i] = (X[:,i] - mean[i]) / std[i]
		weights[:,i] 	= 1 / (err_shift_spectrum[:,i] / std[i]) # may need to include the Fourier power later
		# weights[:, i] = 1 / (err_power_spectrum[:, i] / std[i])

	n_pca = X.shape[1]
	# n_pca = 2
	# Compute the PCA vectors & variance
	pca = WPCA(n_components=n_pca)

	pca.fit(X_new, **kwds)
	pca_score = pca.transform(X_new, **kwds)
	# print(pca_score.shape)  # (632, 17)
	P = pca.components_  # (12, 17)
	# print(P.shape)

	X_wpca = X_new.T
	X_wpca.shape  # (17, 632)

	C = np.zeros((n_pca, X.shape[0]))
	err_C = np.zeros(C.shape)

	W = weights.T

	P = pca.components_.T
	# P.shape  # (17, 12)

	for i in range(X_wpca.shape[1]):
		w = np.diag(W[:, i]) ** 2
		C[:, i] = np.linalg.inv(P.T @ w @ P) @ (P.T @ w @ X_wpca[:, i])

		# the error is discribed by a covariance matrix of C
		Cov_C = np.linalg.inv(P.T @ w @ P)

		# inv(P'*w*P) * P' * w * cov(X(:,i)) * w * P * inv(P'*w*P)
		# Cov_C = np.linalg.inv(P.T @ w @ P) @ P.T @ w @ np.diag(err_shift_spectrum[i,:])**2 @ w @ P @ np.linalg.inv(P.T @ w @ P)

		diag_C = np.diag(Cov_C)
		err_C[:, i] = diag_C ** 0.5

	# determine how many pca scores are needed

	cumulative_variance_explained = np.cumsum(pca.explained_variance_ratio_) * 100 # look again the difference calculated from the other way
	print(cumulative_variance_explained)

	for i in range(len(cumulative_variance_explained)):
		if cumulative_variance_explained[i] < 90:
			n_pca = i
	n_pca += 2
	if cumulative_variance_explained[0]>90:
		n_pca = 1

	print('{:d} pca scores account for {:.2f}% variance explained'.format(n_pca, cumulative_variance_explained[n_pca-1]))

	print('std(C) and midean(err_C) are\n',
		np.around(np.std(C[0:n_pca,:], axis=1), decimals=1), '\n',
		np.around(np.median(err_C[0:n_pca,:], axis=1), decimals=1))

	for i in range(n_pca):
		np.savetxt(str(N_FIESTA)+ '_C' + str(i+1) + '.txt', C[i,:])
		np.savetxt(str(N_FIESTA)+ '_err_C' + str(i+1) + '.txt', err_C[i, :])


	plt.rcParams.update({'font.size': 16})
	fig, axes = plt.subplots(figsize=(18, n_pca*1.5))
	for i in range(n_pca):
		ax = plt.subplot(n_pca,1,i+1)
		if i==0:
			plt.title('N_FIESTA = ' + str(N_FIESTA))
		plt.errorbar(bjd_daily, C[i, :], err_C[i, :], marker='.', ls='none', alpha= 0.5)
		plt.ylabel('PCA%d' %(i+1))
		if i != n_pca-1:
			ax.set_xticks([])
		else:
			plt.xlabel('date_bjd')		
	plt.savefig(str(N_FIESTA)+ '_C.png')
	plt.show()		

