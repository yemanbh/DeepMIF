import os
import pandas as pd

def combine_stat_data(input_dir, output_dir):
	df_combined = pd.DataFrame()
	for file_name in os.listdir(input_dir):
		df = pd.read_csv(os.path.join(input_dir, file_name))
		
		df_combined = pd.concat([df_combined, df], axis=0, sort=False, ignore_index=True)
	
	# save files
	os.makedirs(output_dir, exist_ok=True)
	df_combined.to_csv(os.path.join(output_dir, 'combined_stat_files.csv'), index=False)