from solid import *
import openscad as o
import ase.build
from ase.neighborlist import neighbor_list
import numpy as np
import math
import test as t
from properties import bond_radii
from properties import atom_size
from ase.data import covalent_radii
import subprocess
import os


a = []
b = []
scale = 1.0
def asemolecule():
    a = ase.build.molecule('CH4')
    return a
def databasemolecule():
    a = t.molecule
    return a
def bulkgen():
    a = ase.build.bulk("C",cubic = True)
    a*=(2,2,2)
    a.set_pbc(False)
#   si.write("C:\\Users\\Isaac\\Documents\\Si.xyz")
    return a

a = bulkgen()
b = a.get_positions()
species_uniq = np.unique(a.get_chemical_symbols())
species = a.get_chemical_symbols()
cutoffs = {(str(s1), str(s2)): (bond_radii[s1] + bond_radii[s2]) * 1.1 for s1 in species_uniq for s2 in species_uniq}
size_multiplier = [((atom_size[atom]/31)*0.6*scale) for atom in species]

#print(a)
i, j, d, D = neighbor_list('ijdD', a, cutoffs)
#print(i,j,d,D)
b *= 10*scale
#print(len(b))
count = 0
for x in b:
    o.colour_sphere("red",(4.5*size_multiplier[count]*scale))
    count += 1
# the * before x means it is iterative and uses the x,y and z coordinates in x iteratively
    o.translate(*x*scale)
    print(x)

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
counter = 0;
for i_, j_, d_, D_, u_ in zip(i, j, d, D, u):
    if i_ < j_:
        print(i_, j_, d_, D_, u_)
#        o.cylinder((5 * scale), d_ * 10 * (1.2 * scale))
        o.colour_cylinder("grey", (5 * scale), (d_ * 10 * (1.2 * scale)))
        z_rot()
        y_rot()
        o.rotate(0, y_rot.angle, z_rot.angle)
        o.translate(*b[i_])

o.output(o.result())
#model = import_scad("C:\\Users\\Isaac\\Documents\\python\\working_example\\project.scad")

subprocess.run(["/usr/bin/openscad", "-o", "/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/project.3mf", "project.scad"])
#os.chdir("C:\\Users\\Isaac\\Documents\\python\\working_example\\colorscad-0.5.2")
#subprocess.call([".\\colorscad.sh", "-i", "C:\\Users\\Isaac\\Documents\\python\\working_example\\project.scad",
#                 "-o", "C:\\Users\\Isaac\\Documents\\python\\working_example\\exported models\\colortest.3mf"], shell=True)