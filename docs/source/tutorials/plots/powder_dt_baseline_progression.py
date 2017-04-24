import matplotlib.pyplot as plt
import numpy as np
from skued import gaussian
from skued.baseline import baseline_dt

s, intensity = np.load('powder.npy')

# Double exponential inelastic background and substrate effects
background = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

signal = intensity + background + substrate1 + substrate2

for l in range(5):
	baseline = baseline_dt(signal, level = l, max_iter = 100)
	plt.plot(s, baseline, color = 'r', label = 'Baseline level {}'.format(l))

plt.plot(s, signal, 'k-', label = 'Diffraction')

plt.title('Diffraction pattern of rutile VO$_2$')
plt.xlabel('Scattering angle (1/$\AA$)')
plt.ylabel('Diffracted intensity (counts)')
plt.xlim([s.min(), s.max()])
plt.ylim([0, 100])
plt.show()