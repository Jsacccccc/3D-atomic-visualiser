import argparse

parser = argparse.ArgumentParser(
    prog = "Single_piece_model",
    description = "Creates a basic, colourable, single piece model")
parser.add_argument("--input", help = "this is the input file name")
parser.add_argument("--output", help = "this is the output file name")
parser.add_argument("-scale", help = "sets the scale of the model (default: 1.0)")
parser.add_argument("-cutoff", help = "sets the cutoff length for bonds (default: 1.2)")
args = parser.parse_args()
input_filename = args.input
output_filename = args.output
scale = args.scale
print(input_filename)
print(output_filename)
print(scale)