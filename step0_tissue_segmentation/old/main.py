import os
from tissue_segmentation import tissue_seg
from parse_arguments import get_parsed_arguments

args = get_parsed_arguments()

tissue_seg(input_dir=args.input_dir, output_dir=args.output)