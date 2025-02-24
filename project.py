import openscad as o
from ase.neighborlist import neighbor_list
import numpy as np
import math
import ase.io
from properties import *
import subprocess
from datetime import datetime
import argparse

time_date = datetime.now()
time_date = time_date.strftime("_" + "%Y_" + "%m_" + "%d_" + "%H_" + "%M")

scale = 1.0
cutoff = 1.2

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
# these two have to be floats otherwise they can't be multiplied which breaks the script
scale = float(args.scale)
cutoff = float(args.cutoff)

molecule = ase.io.read(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/real_test_files_xyz/"+ str(input_filename))
#molecule = ase.io.read(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/clean_XYZs/"+"d"+str(input_filename)+".xyz")
#molecule = ase.io.read(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/real_test_files_xyz/Slab_with_dislocation_for3Dprint.xyz")
a = molecule
a.set_pbc(False)
b = a.get_positions()
species_uniq = np.unique(a.get_chemical_symbols())
species = a.get_chemical_symbols()
cutoffs = {(str(s1), str(s2)): (bond_radii[s1] + bond_radii[s2]) * cutoff for s1 in species_uniq for s2 in species_uniq}
size_multiplier = [((atom_size[atom]/31)*0.6*scale) for atom in species]
species_colour = [atom_colours[atom] for atom in species]
print(species_colour)

for x,y in enumerate(range(0,len(size_multiplier))):
    if size_multiplier[x] > 2*scale:
        size_multiplier[x] = 2



i, j, d, D = neighbor_list('ijdD', a, cutoffs)
b *= 10*scale
count = 0
for x in b:
    o.colour_sphere(species_colour[count],(4.5*size_multiplier[count]))
    print(4.5*size_multiplier[count])
    count += 1
# the * before x means it is iterative and uses the x,y and z coordinates in x iteratively
    o.translate(*x)

#cylindrical coordinates, make sure x is zero
#z rotate begins in x axis


def z_rot():
    z_rot.rho = math.sqrt(u_[0] ** 2 + u_[1] ** 2)
    z_rot.angle = math.atan2(u_[1], u_[0]) * (180 / math.pi)

def y_rot():
    y_rot.base = math.sqrt(u_[0] ** 2 + u_[1] ** 2)
    if u_[2] >= 0:
        y_rot.angle = 90 - (math.acos(y_rot.base) * (180 / math.pi))
    elif u_[2] < 0:
        y_rot.angle = 90 + (math.acos(y_rot.base) * (180 / math.pi))


u = D / d[:, None]
#u is the length vector in x y and z directions of each atom
counter = 0
for i_, j_, d_, D_, u_ in zip(i, j, d, D, u):
    if i_ < j_:
        print(i_, j_, d_, D_, u_)
        o.colour_cylinder("grey", (5 * scale), (d_ * 10 * (1.2 * scale)))
        z_rot()
        y_rot()
        o.rotate(0, y_rot.angle, z_rot.angle)
        o.translate(*b[i_])

o.output(o.result())
# this subprocess runs a command line command which renders the scad file generated in this script to be exported to a 3mf file.
subprocess.run(["/home/isaac/Code/colorscad/colorscad.sh", "-i", "/home/isaac/PycharmProjects/3D-atomic-visualiser/project.scad", "-o", "/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/single_body_export_files/" + str(output_filename) + str(time_date) + ".3mf"])
