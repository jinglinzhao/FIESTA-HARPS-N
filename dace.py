# from dace.sun import Sun
# import pandas as pd
# import matplotlib.pyplot as plt
# from astropy.io import fits
# import os
# import numpy as np
# # %matplotlib

# timeseries = Sun.get_timeseries()
# data = pd.DataFrame(timeseries)

# # Show the different timeseries available
# data.keys()

# # Plot the raw radial velocities, in the Solar System barycenter rest frame
# plt.figure()
# plt.plot(data.date_bjd, data.rv_raw, 'o')
# plt.xlabel('BJD - 2400000 [d] [m/s]')
# plt.ylabel('RV [m/s]')
# plt.tight_layout()