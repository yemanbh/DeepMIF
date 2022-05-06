import argparse

def get_parsed_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--input_dir', dest='input_dir', help='directory of the input data')
    parser.add_argument('-o', '--output', dest='output', help='directory of output')
    args = parser.parse_args()

    return args
