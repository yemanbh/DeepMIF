# -*- coding: utf-8 -*-
"""
Created on 30/09/2020

@author: yhagos
"""
import multiprocessing as mp
import os
from skimage.morphology import dilation, erosion, disk
from skimage import io
from skimage.color import rgb2gray
from skimage import measure
import numpy as np
import pandas as pd

# import seaborn as sns
# from skimage.filters import threshold_otsu
# import matplotlib.pyplot as plt

class TissueSeg:
	def __init__(self, input_dir, output_dir, num_processes=1, intensity_threshold=0.843, area_threshold=5000):
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.num_processes = num_processes
		self.intensity_threshold = intensity_threshold
		self.area_threshold = area_threshold
	def tissue_seg(self, images_file_names, p_n):
		# these are parameters optimized from intensity profile () or noisy binary object
		# intensity_threshold = 0.843
		
		
		tissue_area_file = os.path.join(self.output_dir, 'tissue_area.csv')
		
		print(f"--> processing {self.input_dir}")
		os.makedirs(self.output_dir, exist_ok=True)
		
		tissue_area_df = pd.DataFrame(columns=['PatientID', 'TissueArea'])
		
		# initialize vars
		# tissue_area = 0
		# patient_id = os.path.basename(self.input_dir)
		
		
		for file_name in images_file_names:
			# msi_id = os.path.splitext(csv_file_name)[0]
			msi_id = file_name[file_name.find('['): file_name.find(']') + 1]
			
			# get output file path
			output_file_path = os.path.join(self.output_dir, msi_id + '.jpg')
			if os.path.isfile(output_file_path):
				print(msi_id, 'exists in output path')
				continue
			
			file_path = os.path.join(self.input_dir, file_name)
			im = io.imread(file_path)
			
			# rgb to gray
			im_gray = rgb2gray(im)
			# intensity_threshold = threshold_otsu(im_gray)
			
			# binary
			tissue = (im_gray <= self.intensity_threshold) * 1
			
			tissue = dilation(tissue, selem=disk(15))
			tissue = dilation(tissue, selem=disk(15))
			# remove noise using morphological operation
			labeled_im = measure.label(tissue)
			if len(np.unique(labeled_im)) > 1:
				props = measure.regionprops_table(label_image=labeled_im, properties=('label', 'area'))
				tissue_labels = props['label'][props['area'] > self.area_threshold]
				tissue = 1 * np.isin(labeled_im, tissue_labels)
	
			
			# save file
			# tissue = (1 - bcg_clean)
			# tissue_area += np.sum(tissue == 1)
			io.imsave(output_file_path, 255 * tissue.astype('uint8'))
			# fig, ax = plt.subplots(1, 4)
			# ax[0].imshow(im)
			# ax[1].imshow(im_gray)
			# ax[2].imshow(tissue)
			# ax[2].set_title(f"threshold={round(intensity_threshold, 3)}")
			# sns.kdeplot(im_gray.flatten(), ax=ax[3])
			# plt.savefig(file_name)
			# plt.show()
			# plt.show()
			# plt.close()
		
		# tissue_area_df.loc[tissue_area_df.__len__()] = [patient_id, tissue_area]
		# tissue_area_df.to_csv(tissue_area_file, index=False)
		
	def distribute_data(self, file_names, n, num_elem_per_process):
		file_names_list_list = []
		for i in range(n):
			start_ = i * num_elem_per_process
			if i < n - 1:
				files = file_names[start_: start_ + num_elem_per_process]
			else:
				files = file_names[start_:]
			
			file_names_list_list.append(files)
		
		# the last batch might be more
		n_0 = file_names_list_list[0].__len__()
		n_last = file_names_list_list[-1].__len__()
		
		if n_last != n_0:
			print('--> found difference in data distribution. --> making data equally distributed across process')
			diff = n_last - n_0
			to_distribute = file_names_list_list[-1][-diff:]
			
			# distribute them
			for k, file_name in enumerate(to_distribute):
				file_names_list_list[k].append(file_name)
			
			# update last batch
			file_names_list_list[-1] = file_names_list_list[-1][:-diff]
		
		return file_names_list_list
	
	def run_multi_process(self, img_files_names_list):
		if len(img_files_names_list) < self.num_processes:
			n = len(img_files_names_list)
		else:
			n = self.num_processes
		
		num_elem_per_process = int(len(img_files_names_list) // n)
		
		print(
			f"--> there are {img_files_names_list.__len__()} DAPI image --> divided to {n} processes --> each {num_elem_per_process} images ")
		
		file_names_list_list = self.distribute_data(img_files_names_list, n, num_elem_per_process)
		for i in range(n):
			print(f"{i} -- > {file_names_list_list[i]}")
		# create list of processes
		processes = [
			mp.Process(target=self.tissue_seg, args=(file_names_list_list[process_num], process_num))
			for process_num in range(n)
		]
		
		print('{} processes created'.format(n))
		
		# Run processes
		for p in processes:
			p.start()
		
		# Exit the completed processes
		for p in processes:
			p.join()
		print('All Processes finished!!!')
	
	def run(self):
		# print('Processing images in folder : {}'.format(self.subdir_name))
		# markers = self.nucleus_like_proteins + self.other_proteins
		# img_files_list_all_ = os.listdir(os.path.join(self.input_dir, self.subdir_name))
		img_files_list_all = [file_name for file_name in os.listdir(self.input_dir) if 'DAPI' in file_name]
		# img_files_list_all = [file_name for file_name in img_files_list_all_ if 'DAPI' in file_name]
		print(f"images: {img_files_list_all}")
		img_files_names_list = img_files_list_all
		
		if self.num_processes > 1:
			self.run_multi_process( img_files_names_list=img_files_names_list)
		else:
			# print(f"proteins:{self.img_name_pattern}")
			# print(img_files_names_list)
			self.tissue_seg( img_files_names_list, 0)


if __name__ == '__main__':
	pass
	# tissue_seg(input_dir=r'F:\Projects\Vectra_deep_learning\Lorenzo\example\test',
	#            output_dir=r'F:\Projects\Vectra_deep_learning\Lorenzo\profile_example\test')

