import os
# import configs
from step2_coexpression_analysis.package.S2_combine_positive_cells import combine_positive_cells
from step2_coexpression_analysis.package.S4_identify_coexpressing_cells import IdentifyMarkersCoExpression

def run_coexp(input_dir, output_dir, patient_id, proteins, co_exp, co_exp_threshold=3, num_processes=1):
	print('*'*50)
	print('Combine positive markers')
	print('*'*50)
	# combine annotations
	params1 = dict(
	input_dir=input_dir,
	output_dir=os.path.join(output_dir, 'Combined_csv'),
	markers=proteins,
	patient_id=patient_id
	)
	combine_positive_cells(**params1)
	# FOR CO-EXPRESSION, THERE ARE MSI WITH MISSING COMPONENT AND THEY SHOULD BE IGNORED
	# DURING FURTHER ANALYSIS
	print('*'*50)
	print('CO-EXPRESSION ANALYSIS')
	print('*'*50)
	
	# co-expression
	params3 = dict(
		combined_cell_pos_dir=os.path.join(output_dir,  'Combined_csv'),
		output_dir=os.path.join(output_dir,  'Co_expression'),
		threshold=co_exp_threshold,
		coexpression_proteins=co_exp,
		num_processes=num_processes,
		patient_id=patient_id
	)
	coexp = IdentifyMarkersCoExpression(**params3)
	coexp.run_co_exp_analysis()
	
	print('DONE!')

if __name__ == '__main__':
	params = dict(input_dir=r'', # this doen't matter
				  output_dir=r'C:\Users\yhagos\Dropbox (ICR)\Yeman\Projects\DeepMIFUI\co_exp_testing\simulation\data',
				  patient_id='PatientX',
				  proteins=('M1', 'M2', 'M3', 'M4'),
				  co_exp = ('M1+M2', 'M1+M3', 'M1+M2+M3+M4', 'M1+M2+M3', 'M1+M4'),
				  co_exp_threshold=9,
				  num_processes=1)
	run_coexp(**params)
