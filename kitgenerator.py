import openscad as o
from solid import *
import ase.build
from ase.neighborlist import neighbor_list
import numpy as np
import math
import test as t
from properties import bond_radii
from properties import atom_size
import subprocess

# a contains the structure and b contains the positions of atoms in the structure
a = []
b = []
scale = 1.0


#these three definitions define the source used to generate structures.
def asemolecule():
    a = ase.build.molecule('CH3CH2OCH3')
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


a = asemolecule()
b = a.get_positions()
# species_uniq make a tuple of 1 instance of each type of atom
species_uniq = np.unique(a.get_chemical_symbols())
# species makes a tuple of the species of every atom in the structure
species = a.get_chemical_symbols()
cutoffs = {(str(s1), str(s2)): (bond_radii[s1] + bond_radii[s2]) * 1.1 for s1 in species_uniq for s2 in species_uniq}
size_multiplier = [((atom_size[atom]/31)*0.6*scale) for atom in species]

i, j, d, D = neighbor_list('ijdD', a, cutoffs)
b *= 10*scale

place = 0
for x in b:
    #this spawns an atom of a given colour and size based on its element
    o.colour_sphere("red",(4.5*size_multiplier[place]*scale))
    x_value = (1+(place//6))*15*scale
    y_value = (place%6)*15*scale
    z_value = 4.5*size_multiplier[place]*scale
    o.translate(x_value,y_value,z_value)
    place += 1

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
#u is unit length of
counter = 0;
#atomlist = np.arange(5*len(u)).reshape(len(u),3)
atomlist = []
for i_, j_, d_, D_, u_ in zip(i, j, d, D, u):
#   print(i_, j_, d_, D_, u_)
#   o.cylinder((5 * scale), d_ * 10 * (1.2 * scale))
#   o.colour_cylinder("grey", (5 * scale), (d_ * 10 * (1.2 * scale)))
    z_rot()
    y_rot()
#   o.rotate(0, y_rot.angle, z_rot.angle)
#   o.translate(*b[i_]*scale)
    counter = + 1
    atomlist = np.append(atomlist,[i_,j_,d_,z_rot.angle,y_rot.angle], axis=(counter-1))

check = []
counter = 0
#for i in range(len(u)):
#
#    o.difference(
#        if check[]
#        o.sphere()
#    )
#    counter += 1

test = difference()(
    sphere(8),
    translate([-4,-4,2])(cube(8))
)
scad_render_to_file(test, 'solidpytest.scad')
print(atomlist)
print(len(u))
o.output(o.result())