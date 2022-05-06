import os
import pandas as pd
from skimage import io
import numpy as np

class CellDensity:
	def __init__(self, input_dir, output_dir, cell_names, tissue_seg_dir):
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.cell_names = cell_names
		self.tissue_seg_dir = tissue_seg_dir
		
		os.makedirs(self.output_dir, exist_ok=True)
	
	@staticmethod
	def add_dicts(dict1, dict2):
		sum_dict = dict()
		for k in dict1.keys():
			if k in dict2.keys():
				sum_dict[k] = dict1[k] + dict2[k]
			else:
				sum_dict[k] = dict1[k]
		return sum_dict
	def get_cell_density_data(self):
		cell_density_data = {f"{cell_name}":0 for cell_name in self.cell_names}
		patient_id = os.path.basename(self.input_dir)
		cell_density_data['PatientID'] = patient_id
		# tissue_area = pd.read_csv(self.tissue_area_file)
		
		tissue_area = 0
		for msi_file in os.listdir(self.input_dir):
			df = pd.read_csv(os.path.join(self.input_dir, msi_file))
			# count cells
			temp_cells = df.CellPhenotype.value_counts().to_dict()
			cell_density_data = self.add_dicts(cell_density_data, temp_cells)
			
			# get tissue area
		
			msi_id = msi_file[msi_file.find('['): msi_file.find(']') + 1]
			# read tissue mask image
			tissue_mask_image = io.imread(
				os.path.join(self.tissue_seg_dir, msi_id + '.jpg'))
			tissue_area += np.sum(tissue_mask_image > 220)
		
		cell_density_data['TissueAreaPixels^2'] = tissue_area
		# change dict to csv format
		df = pd.DataFrame(cell_density_data, index=[0])
		print(df)
		df.to_csv(os.path.join(self.output_dir, patient_id + '.csv'), index=False)

if __name__ == '__main__':
	output_dir = r'F:\Projects\Vectra_deep_learning\Lorenzo\output_test2'
	xx = "CD4+/CD8+#CD4+/FOXP3+#CD4+/PD1+#CD8+/FOXP3+#CD8+/PD1+#CD4+/CD8-/FOXP3-/PD1-#CD4-/CD8+/FOXP3-/PD1-#CD4-/CD8-/PD1+"
	subdir_name = 'UH13-16987 TIFF'
	stat_params = dict(input_dir=os.path.join(output_dir, 'cd_cc', 'Co_expression', subdir_name),
	                   output_dir=os.path.join(output_dir, 'stat_data'),
	                   cell_names=xx.split('#'),
	                   tissue_seg_dir=os.path.join(output_dir, 'tissue_seg', subdir_name))
	stat = CellDensity(**stat_params)
	stat.get_cell_density_data()

