# -*- coding: utf-8 -*-
"""
Created on 30/09/2020

@author: yhagos
"""
from skimage import io
from skimage.color import rgb2gray
import matplotlib.pyplot as plt
import os
import re
from skimage.morphology import dilation, erosion, disk
from skimage import measure
import numpy as np

from configs import seg_replacement, panel_names, combined_csv_dir, area_threshold, \
	intensity_threshold, output_dir, raw_data_dir



for panel in [panel_names[0]]:
	for patient_id in os.listdir(os.path.join(combined_csv_dir, panel, 'Combined_csv')):
		print(panel, patient_id)
		csv_dir = os.path.join(combined_csv_dir, panel, 'Combined_csv', patient_id)
		csv_file_names = [file_name for file_name in os.listdir(csv_dir) if 'csv' in file_name]
		
		# images dir
		images_dir = os.path.join(raw_data_dir, panel, patient_id)
		images_file_names = [file_name for file_name in os.listdir(images_dir) if seg_replacement[panel]  in file_name]
		for csv_file_name in csv_file_names:
			msi_id = os.path.splitext(csv_file_name)[0]
			
			# output file name
			output_path = os.path.join(output_dir, panel, patient_id)
			os.makedirs(output_path, exist_ok=True)
			output_file_path = os.path.join(output_path, msi_id + '.jpg')
			if os.path.isfile(output_file_path):
				print(msi_id, 'exists in output path')
				continue
			
			dapi_file_name_list = [x for x in images_file_names if msi_id in x]
			
			if len(dapi_file_name_list) > 1 :
				print(panel, patient_id, msi_id)
				raise Exception('multiple dpi files are found for single msi')
			if len(dapi_file_name_list) == 0:
				print(panel, patient_id, msi_id)
				raise Exception('no dpi files is found')

			
			dapi_file_name = dapi_file_name_list[0]
			dapi_file_name = dapi_file_name.replace(' ', '')
			file_path = os.path.join(images_dir, dapi_file_name)
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
				bcg_clean = (255 * bcg_clean).astype('uint8')
			else:
				bcg_clean = bcg
			
			# save file
			io.imsave(output_file_path, bcg_clean)

print('done')


