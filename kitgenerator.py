from solid import *
import openscad as o
import ase.build
from ase.neighborlist import neighbor_list
import numpy as np
import math
import test as t
from properties import bond_radii
from properties import atom_size
import subprocess
import copy

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
count = 1
z_rot_total = [0] * len(i)
y_rot_total = [0] * len(i)
i_total = [0] * len(i)
j_total = [0] * len(i)
for i_, j_, d_, D_, u_ in zip(i, j, d, D, u):
#    print(i_,j_,d_,D_,u_)
    z_rot()
    z_rot_total[count-1] = z_rot.angle
    y_rot()
    y_rot_total[count-1] = y_rot.angle
    i_total[count-1] = int(i_)
    j_total[count-1] = int(j_)
#    print(z_rot.angle)
#    print(y_rot.angle)
    count += 1
print(i_total)
print(j_total)

atom_index = []
x_tot = []
y_tot = []
z_tot = []
place_main = 0
place_sub = 0
for place_main,x in enumerate(b):
    fresh = []
    atom_number = i_total.count(place_main)
    for y in range(0,atom_number):
        fresh.append(j_total[place_sub])
        place_sub += 1
    atom_index.append(fresh)
    x_value = (1 + (place_main // 6)) * 15 * scale
    y_value = (place_main % 6) * 15 * scale
    z_value = 4.5 * size_multiplier[place_main] * scale
    x_tot.append(x_value)
    y_tot.append(y_value)
    z_tot.append(z_value)



count = 0
positive = 0
negative = 0
total = 0
place_sub = 0
for count,x in enumerate(atom_index):
    positive = sphere(4.5*size_multiplier[count]*scale, segments=50)
    for place_sub,y in enumerate(atom_index[count]):
        negative += rotate([0, y_rot_total[place_sub], z_rot_total[place_sub]])(cylinder(r=2,h=5,segments=50))
        negative += rotate([0, y_rot_total[place_sub], z_rot_total[place_sub]])(translate([-20,-20,0.8*4.5*size_multiplier[count]*scale])(cube(40)))
        positive = positive - negative
    total = total + translate([x_tot[count],y_tot[count],z_tot[count]])(positive)

scad_render_to_file(total, '/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/total.scad')


#counter = 0;
#atomlist = np.arange(5*len(u)).reshape(len(u),3)
#atomlist = []
#for i_, j_, d_, D_, u_ in zip(i, j, d, D, u):
#   print(i_, j_, d_, D_, u_)
#   o.cylinder((5 * scale), d_ * 10 * (1.2 * scale))
#   o.colour_cylinder("grey", (5 * scale), (d_ * 10 * (1.2 * scale)))
#    z_rot()
#    y_rot()
#   o.rotate(0, y_rot.angle, z_rot.angle)
#   o.translate(*b[i_]*scale)
#    counter = + 1
#    atomlist = np.append(atomlist,[i_,j_,d_,z_rot.angle,y_rot.angle], axis=(counter-1))

#check = []
#counter = 0
#for i in range(len(u)):
#
#    o.difference(
#        if check[]
#        o.sphere()
#    )
#    counter += 1



#print(atomlist)
#print(len(u))
#o.output(o.result())