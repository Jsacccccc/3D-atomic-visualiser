from solid import *
import openscad as o
import ase.build
from ase.neighborlist import neighbor_list
import numpy as np
import math
import test as t
from properties import *
import subprocess
import argparse

# a contains the structure and b contains the positions of atoms in the structure
a = []
b = []
scale = 0.5


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
def structure_from_file():
    a = t.real_molecule
    a.set_pbc(False)
    return a

a = asemolecule()
b = a.get_positions()
# species_uniq make a tuple of 1 instance of each type of atom
species_uniq = np.unique(a.get_chemical_symbols())
# species makes a tuple of the species of every atom in the structure
species = a.get_chemical_symbols()
cutoffs = {(str(s1), str(s2)): (bond_radii[s1] + bond_radii[s2]) * 1.1 for s1 in species_uniq for s2 in species_uniq}
size_multiplier = [((atom_size[atom]/31)*0.6) for atom in species]
species_colour = [atom_colours[atom] for atom in species]
species_colour_uniq = [atom_colours[atom] for atom in species_uniq]
species_colour_uniq.sort(reverse = False)

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
    x_value = ((place_main // 6)) * 15 * scale
    y_value = (place_main % 6) * 15 * scale
    z_value = 4.5 * size_multiplier[place_main] * scale
    x_tot.append(x_value)
    y_tot.append(y_value)
    z_tot.append(z_value)

# the purpose of this loop is to sort out how many atoms of each colour are being printed
# since ase automatically collects groups of the same species together in the species list, all this does is count them in order
atom_colour_temp = []
atom_colour = []
y = 0
# this loop counts the number of atoms of each colour/species
for x in species_colour_uniq:
    atom_colour_temp.append(species_colour.count(x))
# they are then sorted from least to most common colour/species, as this is the order they are generated in
atom_colour_temp.sort(reverse = False)
temp = 0
for y,z in enumerate(atom_colour_temp):
    temp += atom_colour_temp[y]
    atom_colour.append(temp)


count = 0
positive = 0
negative = 0
atom_total = 0
place_main = 0
place_sub = 0
mini_count = 0
count_sub = 0
size_total = []
# this first loop runs for every species/colour present in the structure where z is equal to the number of atoms of each colour

# this is the loop that prints the spheres for the kit
# this first loop runs for every atom in the structure
print(atom_colour[mini_count])
for count,x in enumerate(atom_index):
    sphere_size = 4.5*size_multiplier[count]*scale
    size_total.append(sphere_size)
    positive = sphere(sphere_size, segments=50)
    negative = 0
# this second loop runs for every bond connection for every atom in the structure, so place_sub is needed to keep track
    for y in atom_index[count]:
# 'negative' is the material being removed from each atom where 'positive' is the atom itself
        negative += rotate([0, y_rot_total[place_sub], z_rot_total[place_sub]])(cylinder(r1 = (1.8 * scale),r2 = (2.5 * scale),h = (10 * scale), segments = 50))
        negative += rotate([0, y_rot_total[place_sub], z_rot_total[place_sub]])(translate([-20*scale,-20*scale,0.8*4.5*size_multiplier[count]*scale])(cube(40 * scale)))
        positive = positive - negative
        place_sub += 1
    count_sub += 1
    atom_total = atom_total + translate([x_tot[count_sub - 1] + max(size_total),y_tot[count_sub - 1] + max(size_total),z_tot[count]])(positive)
    if (count + 1) == atom_colour[mini_count]:
        atom_total = atom_total + cube([x_tot[count] + (2 * max(size_total)),y_tot[count] + (2 * max(size_total)),0.3])
        scad_render_to_file(atom_total, '/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/atoms_to_print'+ str(count) +'.scad')
        subprocess.run(['/usr/bin/openscad', '-o',
                        '/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/kit_export_models/atoms_to_print.3mf',
                        '/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/atoms_to_print.scad'])
        mini_count += 1
        atom_total = 0
        count_sub = 0


x_tot = []
y_tot = []
counter = 0
for counter,any in enumerate(range(0,place_sub)):
    x_value = (0 + (counter // 6)) * 10 * scale
    y_value = (counter % 6) * 10 * scale
    x_tot.append(x_value)
    y_tot.append(y_value)

print(x_tot)
print(y_tot)

counter = 0
bond_total = 0
for i_,j_,d_ in zip(i,j,d):
    if i_ < j_:
#       print(np.float64(d_)*12)
        bond = 0
        bond_num = 0
        bond_len_base = ((np.float64(d_)*12)-(0.85 * 4.5 * (size_multiplier[i_] + size_multiplier[j_]))) * scale
        bond = rotate([0,0,0])(cylinder(r = (2.5 * scale),h = bond_len_base, segments = 50))
        bond += rotate([0,0,0])(translate([0,0,(-2 * scale)])(cylinder(r = (2 * scale) - 0.03,h = ((4 * scale) + bond_len_base),segments = 50)))
        bond_total += translate([x_tot[counter] + (2.5 * scale),y_tot[counter] + (2.5 * scale),(2 * scale)])(bond)
        counter += 1
bond_total += cube([x_tot[counter] + (5 * scale), y_tot[counter] + (5 * scale),0.3])
scad_render_to_file(bond_total,'/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/bonds_to_print.scad' )
subprocess.run(['/usr/bin/openscad', '-o', '/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/kit_export_models/bonds_to_print.3mf', '/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/bonds_to_print.scad'])

#bond_num = rotate([0, 0, 0])(linear_extrude(15)(translate([0, 0, 0])(text(, size = 10))))
#scad_render_to_file(bond_num,'/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/bonds_num_test.scad' )
