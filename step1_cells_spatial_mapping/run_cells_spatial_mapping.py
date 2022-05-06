# -*- coding: utf-8 -*-
"""
Created on 14/04/2020

@author: yhagos
"""
import os
import json

from step1_cells_spatial_mapping.CellSpatialMapping import DetectCells

def run_cd_cc(data_dir, subdir_name, output_dir, proteins, nucleus_like_proteins, other_proteins, num_cpu=1, do_classification=True, image_scale=1):
	cd_model_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models',  'CD_Top1', 'best_model.h5')
	nuclear_cc_model_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models',  'nuclear', 'best_model.h5')
	others_cc_model_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'models',  'non-nuclear', 'best_model.h5')

	#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	# predProbT0.8_DistT11_AreaT18, best performing param combinations for vectra cell detection
	# these params should be optimized for your dataset
	pred_prob_threshold = 0.5
	distance_threshold = 12
	area_threshold = 5
	# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	params = dict(
		cell_detection_model_dir=cd_model_dir,
		others_cc_model_dir=others_cc_model_dir,
		nuclear_cc_model_dir=nuclear_cc_model_dir,
		input_dir=data_dir,
		output_dir=output_dir,
		num_processes=num_cpu,
		scale=image_scale,
		subdir_name=subdir_name,
		split_area_threshold=0.95,
		pred_probability_threshold=pred_prob_threshold,
		postprocess=True,
		cell_names_ordered=('Neg', 'Pos'),
		cell_label_text=os.path.join(os.path.dirname(os.path.abspath(__file__)),
		                             'utils', 'cell_label.txt'),
		do_classification=do_classification,
		prediction_prob_area_threshold=area_threshold,
		distance_threshold=distance_threshold,
		count_loss_used=True,
		img_name_pattern=proteins,
		save_annotated_image=True,
		save_detection_prob_map=True,
		nucleus_like_proteins=nucleus_like_proteins,
		other_proteins=other_proteins,
	)
	
	# save hyper parameters
	os.makedirs(output_dir, exist_ok=True)
	with open(os.path.join(output_dir, '_params_.json'), 'w') as fp:
		json.dump(params, fp, indent=5)
	
	# run
	obj = DetectCells(**params)
	obj.run()
