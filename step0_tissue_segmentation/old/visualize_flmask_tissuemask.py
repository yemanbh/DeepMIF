from skimage import io
import matplotlib.pyplot as plt
from matplotlib import colors
import os
import numpy as np
import seaborn as sns
# sns.set(palette='bright')

composite_im_dir = r'Y:\yhagos\Vectra_deep_learning\Raw-data\Panels\Panel1-Immune-Tcell\14897\3_[51072,8635]_composite_image.tif'
flmask_im_dir = r'Y:\yhagos\Vectra_deep_learning\Data\20200721_Analysis\FL_segmentation\Panel1-Immune-Tcell\Corrected_Segmentation\14897\[51072,8635].jpg'
tissue_mask_im_dir = r'Y:\yhagos\Vectra_deep_learning\Data\20200721_Analysis\bcg_segmentation\Panel1-Immune-Tcell\14897\[51072,8635].jpg'

t_mask_im = io.imread(tissue_mask_im_dir)
composite_im = io.imread(composite_im_dir)
fl_mask_im = io.imread(flmask_im_dir)

fig, ax = plt.subplots(2, 3)
ax[0, 0].imshow(composite_im, cmap='gray')
ax[0, 0].set_title('composite')
ax[0, 1].imshow(t_mask_im, interpolation='nearest')
ax[0, 1].set_title('background')

ax[0, 2].imshow(fl_mask_im>0, interpolation='nearest')
ax[0, 2].set_title('Inside FL (original)')

# inside and outside follicule after removing background pixels
t_mask_im = (t_mask_im > 0) * 1
t_inside_fl = ((fl_mask_im > 0) - t_mask_im) > 0
ax[1, 1].imshow(t_inside_fl>0, interpolation='nearest')
ax[1, 1].set_title('Inside FL (final)')

t_outside_fl = ((fl_mask_im == 0) - t_mask_im) > 0
seg_image = np.zeros(t_inside_fl.shape, dtype='uint8')
seg_image[np.where(t_inside_fl)] = 10
seg_image[np.where(t_outside_fl)] = 20

# ax[1, 0].imshow(fl_mask_im == 0, interpolation='nearest')
cmap = plt.get_cmap("Dark2")
cmap = cmap([0, 1, 2])
cmap_name = 'my_list'
n_bin = 3
cm = colors.LinearSegmentedColormap.from_list(
        cmap_name, cmap, N=n_bin)
# norm = colors.BoundaryNorm(bounds, len(cmap))

ax[1, 0].imshow(seg_image, interpolation='nearest', cmap=cm)
ax[1, 0].set_title('Outside FL (original)')

ax[1, 2].imshow(t_outside_fl>0, interpolation='nearest')
ax[1, 2].set_title('Outside FL (final)')
for i in range(2):
	for j in range(3):
		ax[i, j].axis('off')

plt.tight_layout()
plt.savefig(os.path.join('output\sample_segmentation.jpg'), dpi=600)
plt.show()