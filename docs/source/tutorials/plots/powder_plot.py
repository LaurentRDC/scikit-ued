import matplotlib.pyplot as plt
import numpy as np

s, intensity = np.load('powder.npy')
plt.plot(s, intensity, 'k-')
plt.title('Background-subtracted diffraction pattern of rutile VO$_2$')
plt.xlabel('Scattering angle (1/$\AA$)')
plt.ylabel('Diffracted intensity (counts)')
plt.xlim([s.min(), s.max()])
plt.ylim([0, 20])
plt.show()