
import os
import sys
import time
from parse_arguments import get_parsed_arguments
from step1_cells_spatial_mapping.run_cells_spatial_mapping import run_cd_cc
from step2_coexpression_analysis.run_coexpression_analysis import run_coexp
from step3_mark_composite_images.AnnotateMifImage import AnnotateMifImage
from step0_tissue_segmentation.tissue_segmentation_V1 import TissueSeg
from step4_cell_aboundance_stat.CellDensityData import CellDensity
from utils import combine_stat_files

if __name__ == '__main__':
	
	print('arguments passed', sys.argv)
	args = get_parsed_arguments()
	
	input_dir = "/input"
	output_dir = "/output"
	print(args.nuclear_markers)
	print(args.non_nuclear_markers)
	nuclear_markers = str(args.nuclear_markers).split('#')
	other_markers = str(args.non_nuclear_markers).split('#')
	co_expression_phenotype = str(args.cell_phenotypes).split('#')
	co_exp_dist_threshold = int(args.distance)
	subdir_name = int(args.subdir_name)
	
	# these are constant variables
	#ToDo: this applies in windows only; need to be updated for Ubuntu
	print('*'*40)
	print(f'Docker is using {os.cpu_count()} cpu')
	print('*'*40)
	
	num_cpu = os.cpu_count()
	spatial_analysis = False
	mask_dir = None

	print('*'*40)
	print("Running cell detection and classification with these parameters:")

	print(f'input dir:{input_dir} \n '
		  f'output dir: {output_dir} \n '
		  f'nuclea markers: {nuclear_markers} \n non-nuclear markers: {other_markers} \n '
		  f'cell phenotypes: {co_expression_phenotype} \n co-expression distance: {co_exp_dist_threshold}')
	print('*'*40)

	# data_dir = os.path.dirname(data_dir)
	# subdir_name = os.path.basename(data_dir)
	data_dir = input_dir
	slides_list = os.listdir(input_dir)
	# for k, subdir_name in enumerate(slides_list):
	print(f'evaluating slide:{subdir_name}')
	

	t0 = time.time()
	# tissue segmentation
	seg_params = dict(input_dir=os.path.join(data_dir, subdir_name),
	           num_processes=num_cpu,
	           output_dir=os.path.join(output_dir, 'tissue_seg', subdir_name))
	seg = TissueSeg(**seg_params)
	seg.run()
	
	# cell detection and classification on de-convoluted images
	run_cd_cc(output_dir=os.path.join(output_dir, 'cd_cc'),
			  subdir_name=subdir_name,
			  data_dir=data_dir,
			  nucleus_like_proteins=nuclear_markers,
			  other_proteins=other_markers,
			  proteins=other_markers + nuclear_markers,
			  num_cpu=num_cpu)

	# co-expression analysis
	run_coexp(input_dir=os.path.join(output_dir, 'cd_cc', 'AnnotatedCellsCoord'),
			  output_dir=os.path.join(output_dir, 'cd_cc'),
			  patient_id=subdir_name,
			  co_exp_threshold=co_exp_dist_threshold,
			  proteins=other_markers + nuclear_markers,
			  co_exp=co_expression_phenotype,
			  num_processes=num_cpu)

	# annotation on mIF images
	annotation_params = dict(
		data_dir=data_dir,
		subdir_name=subdir_name,
		# cell_phenotypes=nuclear_markers + other_markers + co_expression_phenotype,
		cell_phenotypes=co_expression_phenotype,
		output_dir=os.path.join(output_dir, 'cd_cc', 'annotated_composite'),
		csv_dir=os.path.join(output_dir, 'cd_cc', 'Co_expression'),
		num_cpu=num_cpu)

	annotation_obj = AnnotateMifImage(**annotation_params)
	annotation_obj.run()
	time_elapsed = (time.time() - t0) / 60
	with open(os.path.join(output_dir, 'Time_elapsed.txt'), 'a') as fp:
		fp.write(f"{subdir_name}: {time_elapsed} minutes\n")
	
	# generate stat data
	stat_params = dict(input_dir=os.path.join(output_dir, 'cd_cc', 'Co_expression', subdir_name),
	                  output_dir=os.path.join(output_dir, 'stat_data'),
	                  cell_names=co_expression_phenotype,
	                  tissue_seg_dir=os.path.join(output_dir, 'tissue_seg', subdir_name))
	stat = CellDensity(**stat_params)
	stat.get_cell_density_data()
	
	# what percent of cases are analysed
	percent_covered = (k + 1) / slides_list.__len__()
	sys.stdout.write(f"Total complete: {int(percent_covered * 100)}%\n")

	combine_stat_files.combine_stat_data(input_dir=os.path.join(output_dir, 'stat_data'),
	                                     output_dir=os.path.join(output_dir, 'stat_data_combined'))
