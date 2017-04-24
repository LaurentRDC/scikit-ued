import matplotlib.pyplot as plt
import numpy as np
from skued import gaussian

s, intensity = np.load('powder.npy')

# Double exponential inelastic background and substrate effects
background = 75 * np.exp(-7 * s) + 55 * np.exp(-2 * s)
substrate1 = 0.8 * gaussian(s, center = s.mean(), fwhm = s.mean()/4)
substrate2 = 0.9 * gaussian(s, center = s.mean()/2.5, fwhm = s.mean()/4)

plt.plot(s, intensity + background + substrate1 + substrate2, 'k-',
		 label = 'Diffraction')
plt.plot(s, background + substrate1 + substrate2, 'r-', 
		 label = 'Baseline')

plt.legend()

plt.title('Diffraction pattern of rutile VO$_2$')
plt.xlabel('Scattering angle (1/$\AA$)')
plt.ylabel('Diffracted intensity (counts)')
plt.xlim([s.min(), s.max()])
plt.ylim([0, 100])
plt.show()