# -*- coding: utf-8 -*-
"""
Created on 15/05/2020

@author: yhagos
"""
import pandas as pd
import os
import numpy as np
import itertools
from scipy.spatial.distance import cdist
import multiprocessing as mp
pd.options.mode.chained_assignment = None

class IdentifyMarkersCoExpression:
	def __init__(self, combined_cell_pos_dir, patient_id, output_dir, threshold, coexpression_proteins, num_processes=1):
		self.combined_cell_pos_dir = combined_cell_pos_dir
		self.patient_id = patient_id
		self.output_dir = output_dir
		self.threshold = threshold
		self.coexpression_proteins = coexpression_proteins
		self.num_processes = num_processes


	def identify_co_expressing_cells(self, file_names_list, process_num):
		# create output file names and check if they exist in the path
		save_dir = os.path.join(self.output_dir, self.patient_id)
		os.makedirs(save_dir, exist_ok=True)

		for n_, file_name in enumerate(file_names_list):
			print('Process:{}, Patient Id:{}, File name:{}, {}/{}'.format(process_num + 1, self.patient_id, file_name,
																		  n_ + 1, len(file_names_list)))
			msi_name = os.path.splitext(file_name)[0]
			output_csv = os.path.join(save_dir, msi_name + '_co_exp.csv')

			# if os.path.isfile(output_csv) :
			# 	print('{} already exists'.format(output_csv))
			# 	continue

			cell_data_df = pd.read_csv(os.path.join(self.combined_cell_pos_dir, self.patient_id, file_name))
			col_names = cell_data_df.columns

			cell_data_df_copy = cell_data_df.copy()
			overlap_df_all = pd.DataFrame(columns=col_names)
			# overlap_markes_pos = dict()

			for co_exp in self.coexpression_proteins:
				
				co_expression_available = True
				coexp_protein_list = co_exp.split('+')
				# if len(coexp_protein_list) < 3:
				# 	print('stop here')
				# 	continue

				# empty index
				coexp_markers_index_database = {f"{protein}": [] for protein in coexp_protein_list}

				protein_pairs = list(itertools.combinations(coexp_protein_list, 2))
				for protein_pair in protein_pairs:
					# print(protein_pair)
					# if protein_pair[0] == 'B' and protein_pair[1] == 'C':
					# 	print('stop here')
					# protein 1 data frame
					protein_1_data = cell_data_df.loc[cell_data_df['Component'] == protein_pair[0]].reset_index()
					protein_1_data = protein_1_data.rename(columns={'index': 'INDEX_'})

					# for more than 2 markers expression if there is data from previous computation; consider it
					if coexp_markers_index_database[protein_pair[0]].__len__() != 0:
						protein_1_data = protein_1_data.loc[protein_1_data['INDEX_'].isin(coexp_markers_index_database[protein_pair[0]]), :]
					else:
						pass


					# protein 2 data frame
					protein_2_data = cell_data_df.loc[cell_data_df['Component'] == protein_pair[1]].reset_index()
					protein_2_data = protein_2_data.rename(columns={'index': 'INDEX_'})

					if coexp_markers_index_database[protein_pair[1]].__len__() != 0:
						protein_2_data = protein_2_data.loc[protein_2_data['INDEX_'].isin(coexp_markers_index_database[protein_pair[1]]), :]
					else:
						pass
					

					overlap_index_input1, overlap_index_input2 = self.get_co_exp_cells_detail(protein_1_data, protein_2_data)
					if overlap_index_input1.__len__() == 0:
						co_expression_available = False
						break
					indexs_dict = dict()
					indexs_dict[protein_pair[0]] = overlap_index_input1
					indexs_dict[protein_pair[1]] = overlap_index_input2
					coexp_markers_index_database = self.update_coexpression_database(coexp_markers_index_database, indexs_dict)

				# update which is overlapping and not
				if co_expression_available:
					overlapping_indices = self.get_index_co_expressing_markers_position(coexp_markers_index_database)
					cell_data_df_copy.loc[overlapping_indices, 'Component'] = 'co_exp'
	
					# get overlap data
					overlap_df = self.get_overlap_data(coexp_database=coexp_markers_index_database, data=cell_data_df_copy.copy())
					overlap_df['Component'] = co_exp
					overlap_df_all = pd.concat([overlap_df_all, overlap_df], ignore_index=True, axis=0, sort=False)
					
				else:
					pass

			cell_data_df_copy.drop(columns=['Class'], inplace=True)
			overlap_df_all.drop(columns=['Class'], inplace=True)

			# drop all cells co-expressing different markers from cell_data_df_copy
			# cell_data_df_copy.drop(cell_data_df_copy.index[cell_data_df_copy['Component'] == 'co_exp'],	 inplace=True)
			non_overlap_df_data = cell_data_df_copy.loc[cell_data_df_copy['Component'] != 'co_exp', :]

			# concatenate single marker expressing cells and co-expressing cells
			# combined_df_all = pd.concat([overlap_df_all, cell_data_df_copy], ignore_index=True, axis=0, sort=False)
			combined_df_all = pd.concat([overlap_df_all, non_overlap_df_data], ignore_index=True, axis=0, sort=False)
			combined_df_all.to_csv(output_csv, index=False)

	def get_overlap_data(self, coexp_database: dict, data: pd.DataFrame) -> pd.DataFrame:
		df_overlap = pd.DataFrame()
		for protein, index_values in coexp_database.items():
			df = data.iloc[index_values, :].reset_index(drop=True)
			if df_overlap.__len__() == 0:
				df_overlap = df
			else:
				# summation
				df_overlap[['X', 'Y']] = df_overlap[['X', 'Y']] + df[['X', 'Y']]

		# average
		df_overlap[['X', 'Y']] = df_overlap[['X', 'Y']] // coexp_database.__len__()

		return df_overlap

	def get_index_co_expressing_markers_position(self, coexp_database: dict) -> list:
		index_list = []
		for _, val in coexp_database.items():
			index_list += val

		return index_list

	def update_coexpression_database(self, database, current_computation):
		updated_database = database.copy()
		if self.__database_empty(database):

			for protein, values in current_computation.items():
					updated_database[protein] = current_computation[protein]
		else:

			# do the update using the non-empty field
			for protein, values in current_computation.items():
				if database[protein].__len__() != 0:
					# get index of common values
					common_values_index = self.__get_common_values_index(values1=database[protein],
																		 values2=current_computation[protein])
					# update dictionary
					updated_database = self.__do_updates(database, common_values_index)
					break
			# the above loop only will look only for one protein: so make sure you update here
			for protein, values in current_computation.items():
				updated_database[protein] = values

		return updated_database

	def __do_updates(self, database: dict, index_values_included: list) -> dict:
		new_db = dict()
		for protein, values in database.items():
			if values.__len__() != 0:
				new_db[protein] = [values[index_val] for index_val in index_values_included]
			else:
				new_db[protein] = []

		return new_db

	def __database_empty(self, database):
		vals = []
		for _, val in database.items():
			vals += val
		return vals.__len__() == 0

	def __get_common_values_index(self, values1: list, values2: list) -> list:
		common_values_index = []
		for val in values2:
			common_values_index.append(values1.index(val))
		
		return common_values_index

	def get_co_exp_cells_detail(self, df1, df2):

		# if either of cell_data_df_p1 or cell_data_df_p2 are empty
		if len(df1) == 0 or len(df2) == 0:
			return [], []
		else:
			# compute euclidean distance between the cell pos
			euclidean_dist = cdist(df1[['X', 'Y']].to_numpy(), df2[['X', 'Y']].to_numpy(), metric='euclidean')
			# find the location where in 2d space, which is the minimum distance
			arg_min_dist = np.argmin(euclidean_dist, axis=1)
			is_argmin = np.full(euclidean_dist.shape, fill_value=False, dtype=bool)
			is_argmin[(np.array(np.arange(euclidean_dist.shape[0])), np.array(arg_min_dist))] = True

			# distance and threshold
			is_overlap = euclidean_dist <= self.threshold

			# masking: identify overlapping expression level
			is_overlap = np.logical_and(is_overlap, is_argmin)

			df1_index, df2_index = np.where(is_overlap)

			# THESE ARE INDICES IN THE ORIGINAL CSV FILES
			overlap_index_input1 = df1.iloc[df1_index, :]['INDEX_'].to_list()
			overlap_index_input2 = df2.iloc[df2_index, :]['INDEX_'].to_list()

			# return overlap_df, overlap_pos_dict
			return overlap_index_input1, overlap_index_input2

	def run_co_exp_analysis(self):

		file_lists = os.listdir(os.path.join(self.combined_cell_pos_dir, self.patient_id))

		if self.num_processes > 1:
			n = len(file_lists)
			if n < self.num_processes:
				num_processes = n
			num_elem_per_process = int(np.ceil(n / self.num_processes))

			file_names_list_list = []

			for i in range(self.num_processes):
				start_ = i * num_elem_per_process
				x = file_lists[start_: start_ + num_elem_per_process]
				file_names_list_list.append(x)

			print('{} processes created.'.format(self.num_processes))
			# create list of processes
			processes = [
				mp.Process(target=self.identify_co_expressing_cells,
						   args=(file_names_list_list[process_num],
								 process_num)) for process_num in range(self.num_processes)]

			print('processes created')

			# Run processes
			for p in processes:
				p.start()

			# Exit the completed processes
			for p in processes:
				p.join()
			print('All Processes finished!!!')
		else:
			self.identify_co_expressing_cells(file_lists, 0)


if __name__ == '__main__':
	pass