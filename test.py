import numpy as np
import skued
from scipy.optimize import curve_fit

# Load data from file first
# 2 x N array:
#   first row is time-delay
#   second row is diffracted intensity
block = np.load('docs\\tutorials\\tseries1.npy')
time, intensity = block[0, :], block[1, :]

# Compute initial guesses for this curve (optional)
initial_guesses = (0,                                   # time-zero
                    intensity.min() - intensity.max(),   # amplitude
                    1,                                   # time-constant
                    intensity.max())                     # offset

params, pcov = curve_fit(skued.exponential, time, intensity, 
                            p0 = initial_guesses)
                                                        
tzero, amplitude, tconst, offset = params
best_fit_curve = skued.exponential(time, *params)

import matplotlib.pyplot as plt
plt.plot(time, intensity, '.k')
plt.plot(time, best_fit_curve, '-r')
plt.show()