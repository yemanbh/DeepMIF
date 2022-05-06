from skimage import io
import matplotlib.pyplot as plt
from matplotlib import colors
import os
import numpy as np
import seaborn as sns
# sns.set(palette='bright')
msi_id = '[54735,16839]'
composite_im_dir = r'Y:\yhagos\Vectra_deep_learning\Raw-data\Panels\Panel1-Immune-Tcell\14897\3_{}_composite_image.tif'.format(msi_id)
flmask_im_dir = r'Y:\yhagos\Vectra_deep_learning\Data\20200721_Analysis\FL_segmentation\Panel1-Immune-Tcell\Corrected_Segmentation\14897\{}.jpg'.format(msi_id)
tissue_mask_im_dir = r'Y:\yhagos\Vectra_deep_learning\Data\20200721_Analysis\bcg_segmentation_final\Panel1-Immune-Tcell\14897\{}.jpg'.format(msi_id)

# color map
t_mask_im = io.imread(tissue_mask_im_dir)
composite_im = io.imread(composite_im_dir)
fl_mask_im = io.imread(flmask_im_dir)

fig, ax = plt.subplots(2, 2)
ax[0, 0].imshow(composite_im, cmap='gray')
ax[0, 0].set_title('composite image')
cmap = plt.get_cmap("tab20c")
cmap = cmap([10, 1])
cm = colors.LinearSegmentedColormap.from_list(name='my_list', colors=cmap, N=2)

ax[0, 1].imshow(t_mask_im, interpolation='nearest', cmap=cm)
ax[0, 1].set_title('background')
cmap = plt.get_cmap("tab20c")
cmap = cmap([14, 5])
cm = colors.LinearSegmentedColormap.from_list(name='my_list', colors=cmap, N=2)
ax[1, 0].imshow(fl_mask_im>0, interpolation='nearest', cmap=cm)
ax[1, 0].set_title('Manual FL segmentation')

# inside and outside follicule after removing background pixels
t_mask_im = (t_mask_im > 0) * 1
t_inside_fl = ((fl_mask_im > 0) - t_mask_im) > 0


t_outside_fl = ((fl_mask_im == 0) - t_mask_im) > 0
seg_image = np.zeros(t_inside_fl.shape, dtype='uint8')
seg_image[np.where(t_inside_fl)] = 10
seg_image[np.where(t_outside_fl)] = 20

cmap = plt.get_cmap("tab20c")
cmap = cmap([1, 5, 10])
cm = colors.LinearSegmentedColormap.from_list(name='my_list', colors=cmap, N=3)
ax[1, 1].imshow(seg_image, interpolation='nearest', cmap=cm)
ax[1, 1].set_title('Outside FL (original)')

for i in range(2):
	for j in range(2):
		ax[i, j].axis('off')

plt.tight_layout()
plt.savefig('output\{}.pdf'.format(msi_id))
plt.show()