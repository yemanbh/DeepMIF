import argparse

def get_parsed_arguments():
    parser = argparse.ArgumentParser()
    # ['-distance', value, '-input_dir', value, '-output_dir', value, '-cell_phenotypes', value '-nuclear_m', value, '-non_nuclear_markers', value]
    parser.add_argument('-distance', dest='distance')
    parser.add_argument('-cell_phenotypes', dest='cell_phenotypes')
    parser.add_argument('-nuclear_m', dest='nuclear_markers')
    parser.add_argument('-non_nuclear_m', dest='non_nuclear_markers')
	parser.add_argument('-subdir_name', dest='folder name')
    parser.add_argument('-n', '--num_cpu', dest='num_cpu', help='number of cpu for multiprocessing', type=int, default=1)
    parser.add_argument('-s', '--spatial', dest='spatial_analysis', help='apply spatial analysis or not', default=False, action='store_true')
    args = parser.parse_args()

    return args
