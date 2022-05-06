# -*- coding: utf-8 -*-
"""
Created on 30/09/2020

@author: yhagos
"""
from skimage import io
from skimage.color import rgb2gray
import os
from skimage.morphology import dilation, erosion, disk
from skimage import measure
import numpy as np
import pandas as pd

def tissue_seg(input_dir, output_dir):

	# these are parameters optimized from intensity profile or noisy binary object
	intensity_threshold = 0.007
	area_threshold = 5000
	
	tissue_area_file = os.path.join(output_dir, 'tissue_area.csv')

	print(f"--> processing {input_dir}")
	os.makedirs(output_dir, exist_ok=True)
	
	tissue_area_df = pd.DataFrame(columns=['PatientID', 'TissueArea'])
	
	# initialize vars
	tissue_area = 0
	patient_id = os.path.basename(input_dir)
	images_file_names = [file_name for file_name in os.listdir(input_dir) if 'composite_image'  in file_name]
	for file_name in images_file_names:
		# msi_id = os.path.splitext(csv_file_name)[0]
		msi_id = file_name[file_name.find('[') : file_name.find(']') + 1]
		
		
		# get output file path
		output_file_path = os.path.join(output_dir, msi_id + '.jpg')
		if os.path.isfile(output_file_path):
			print(msi_id, 'exists in output path')
			continue

		file_path = os.path.join(input_dir, file_name)
		im = io.imread(file_path)
		
		# rgb to gray
		im_gray = rgb2gray(im)
		
		# binary
		bcg = (im_gray <= intensity_threshold) * 1
		
		# remove noise using morphological operation
		labeled_im = measure.label(bcg)
		if len(np.unique(labeled_im)) > 1:
			props = measure.regionprops_table(label_image=labeled_im, properties=('label', 'area'))
			bcg_label = props['label'][props['area'] > area_threshold]
			bcg_clean = np.isin(labeled_im, bcg_label)
			bcg_clean = dilation(bcg_clean, selem=disk(2))
			bcg_clean = erosion(bcg_clean, selem=disk(2))
			
		else:
			bcg_clean = bcg
		
		# save file
		tissue = (1 - bcg_clean)
		tissue_area += np.sum(tissue == 1)
		io.imsave(output_file_path, 255*tissue.astype('uint8'))
	

	tissue_area_df.loc[tissue_area_df.__len__()] = [patient_id, tissue_area]
	tissue_area_df.to_csv(tissue_area_file, index=False)
	
	