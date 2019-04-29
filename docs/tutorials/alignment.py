# This code used to be part of the documentation
# and rendered with sphinx
#
# However, the masked_register_translation function requires
# too much memory, and so readthedocs would kill the documentation
# build. Therefore, we render the image locally instead.

import matplotlib.pyplot as plt
import numpy as np
from skimage.feature import masked_register_translation
from skued import shift_image, diffread

ref = diffread("Cr_1.tif")
im = diffread("Cr_2.tif")

mask = np.ones_like(ref, dtype=np.bool)
mask[0:1250, 950:1250] = False

shift = masked_register_translation(im, ref, mask)
shifted = shift_image(im, -1 * shift)

fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(nrows=2, ncols=3, figsize=(9, 6))
ax1.imshow(ref, vmin=0, vmax=200)
ax2.imshow(im, vmin=0, vmax=200)
ax3.imshow(ref - im, cmap="RdBu_r", vmin=-100, vmax=100)
ax4.imshow(mask, vmin=0, vmax=1, cmap="binary")
ax5.imshow(shifted, vmin=0, vmax=200)
ax6.imshow(ref - shifted, cmap="RdBu_r", vmin=-100, vmax=100)

for ax in (ax1, ax2, ax3, ax4, ax5, ax6):
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)

ax1.set_title("Reference")
ax2.set_title("Data")
ax3.set_title("Difference")
ax4.set_title("Mask")
ax5.set_title("Aligned data")
ax6.set_title("Diff. after shift")

plt.tight_layout()
plt.savefig("cr_alignment.png")
