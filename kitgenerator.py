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
from tqdm import tqdm
import time
from datetime import datetime

parser = argparse.ArgumentParser(
    prog = "Single_piece_model",
    description = "Creates a basic, colourable, single piece model")
parser.add_argument("--input", help = "this is the input file name")
#parser.add_argument("--output", help = "this is the output file name")
parser.add_argument("-scale", help = "sets the scale of the model (default: 1.0)")
parser.add_argument("-cutoff", help = "sets the cutoff length for bonds (default: 1.2)")
args = parser.parse_args()
input_filename = args.input
#output_filename = args.output
# these two have to be floats otherwise they can't be multiplied which breaks the script
scale = float(args.scale)
cutoff = float(args.cutoff)

# generates the time and date to be concatenated into output file names
time_date = datetime.now()
time_date = time_date.strftime("_" + "%Y_" + "%m_" + "%d_" + "%H_" + "%M_")
print(time_date)

# reads the user's file and assigns the molecule to a
molecule = ase.io.read(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/real_test_files_xyz/"+ str(input_filename))
#molecule = ase.io.read(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/clean_XYZs/"+"d"+str(input_filename)+".xyz")
a = molecule
# needed to sort atoms into groups of the same species, otherwise the kit tries to put different species of atoms together.
a = a[a.numbers.argsort()]
# b grabs the coordinates of every atom in the structure
b = a.get_positions()
# species_uniq make a tuple of 1 instance of each type of atom
species_uniq = np.unique(a.get_chemical_symbols())
# species makes a tuple of the species of every atom in the structure
species = a.get_chemical_symbols()
print(species)
cutoffs = {(str(s1), str(s2)): (bond_radii[s1] + bond_radii[s2]) * 1.1 for s1 in species_uniq for s2 in species_uniq}
size_multiplier = [((atom_size[atom]/31)*0.6) for atom in species]
species_colour = [atom_colours[atom] for atom in species]
species_colour_uniq = [atom_colours[atom] for atom in species_uniq]
species_colour_uniq.sort(reverse = False)
#neighbour list used to list all
i, j, d, D = neighbor_list('ijdD', a, cutoffs)
b *= 10*scale
subprocess.run(['mkdir','/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/kit_export_models/'+'printing_files_'+str(input_filename)+str(time_date)])

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
atom_colour_temp = [0] * len(species_colour_uniq)
atom_colour = []
y = 0
sub_loop = 0
# this loop/section counts the number of atoms of each colour/species and first takes the number of each species then sorts them from most to least common and compounds them as this is needed for indexing while generating atoms.
atom_colour_temp[len(species_colour_uniq) - 1] = species_colour[len(species_colour_uniq) - 1].count(species_colour[len(species_colour_uniq) - 1])

for x in range(1,len(species)):
    print(atom_colour_temp)
    atom_colour_temp[sub_loop] += 1
    if species_colour[x] != species_colour[x - 1]:
        sub_loop += 1
#extra addition needed to the last one because the range starts at 1 so it misses the last element out

# this goes through atom_colour_temp and adds all of the previous values to the current one, so the final element is equal to the number of atoms.
temp = 0
for y,z in enumerate(atom_colour_temp):
    temp += atom_colour_temp[y]
    atom_colour.append(temp)
print(atom_colour)

# this is ugly but all of these need to be set to zero before atoms can be generated.
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
progress_bar = tqdm(total = len(atom_index), desc = "Atom rendering progress: ")
print(y_tot)
for count,x in enumerate(atom_index):
    sphere_size = 4.5*size_multiplier[count]*scale
    size_total.append(sphere_size)
    positive = sphere(sphere_size, segments=50)
    negative = 0
# this second loop runs for every bond connection for every atom in the structure, so place_sub is needed to keep track
    for y in atom_index[count]:
# 'negative' is the material being removed from each atom where 'positive' is the atom
# joint slot is recessed by 2.1*scale into the atom since the pin on the bond is 2*scale units and there needs to be some leeway for it to fit into.
# multiplied by 0.8 since that is where the flat surface formed by the sphere is - specifically recessed from that rather than the radius of the atom as calculated above
        negative += rotate([0, y_rot_total[place_sub], z_rot_total[place_sub]])(translate([0,0,(0.8*4.5*size_multiplier[count]*scale)-(2.1*scale)])(cylinder(r1 = (2.3 * scale),r2 = (2.5 * scale),h = (2.11 * scale), segments = 50)))
        negative += rotate([0, y_rot_total[place_sub], z_rot_total[place_sub]])(translate([-20*scale,-20*scale,0.8*4.5*size_multiplier[count]*scale])(cube(40 * scale)))
        positive = positive - negative
        place_sub += 1
    count_sub += 1
    positive += rotate([0,0,90])(translate([-5 * scale,-5 * scale,- z_tot[count]])(linear_extrude(height = 1)(text(size = 5, text = str(count)))))
    atom_total += translate([x_tot[count_sub - 1] + max(size_total),y_tot[count_sub - 1] + max(size_total),z_tot[count]])(positive)
    time.sleep(0.1)
    progress_bar.update(1)
    if (count + 1) == atom_colour[mini_count]:
        print(mini_count)
        print(count_sub)
        if count_sub >= 6:
            atom_total = atom_total + cube([x_tot[count_sub - 1] + (2 * max(size_total)), (75.0 * scale) + (2 * max(size_total)), 0.2])
        else:
            atom_total = atom_total + cube([x_tot[count_sub - 1] + (2 * max(size_total)),y_tot[count_sub - 1] + (2 * max(size_total)),0.2])
        atom_total = atom_total + translate([x_tot[count_sub - 1] + (2 * max(size_total)),0,0])(cube([10,10,0.3])) + translate([x_tot[count_sub - 1] + (2 * max(size_total)) + 1.25,1.25,0])(linear_extrude(height = 1)(text(str(species[count]), 7.5)))
        scad_render_to_file(atom_total, '/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/atoms_to_print'+ str(time_date) + str(mini_count) + '_' + '.scad')
        subprocess.run(['/usr/bin/openscad', '-o',
                        '/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/kit_export_models/'+'printing_files_'+str(input_filename)+str(time_date)+'/'+'atoms_to_print'+ str(time_date) + str(mini_count) + '_' +'.3mf',
                        '/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/atoms_to_print'+ str(time_date) + str(mini_count) + '_' + '.scad'])
        mini_count += 1
        atom_total = 0
        count_sub = 0


# coordinate values need to be reset and recalculated for the bonds
x_tot = []
y_tot = []
counter = 0
for counter,any in enumerate(range(0,place_sub)):
    x_value = (0 + (counter // 6)) * 10 * scale
    y_value = (counter % 6) * 10 * scale
    x_tot.append(x_value)
    y_tot.append(y_value)


# this loop generates the bonds and renders/exports them - i_,j_,d_ calculated again for simplicity sake rather than referencing them
counter = 0
bond_total = 0
bond_progress = tqdm(total = (len(i)/2), desc = "Bond rendering progress: ")
# i lists the number of an atom once for every bond it has, e.g carbons show up as 4 of their index number

for i_,j_,d_ in zip(i,j,d):
    if i_ < j_:
        print(species[i_],species[j_])
#       print(np.float64(d_)*12)
        bond = 0
        bond_num = 0
        bond_len_base = ((np.float64(d_)*12)-(0.85 * 4.5 * (size_multiplier[i_] + size_multiplier[j_]))) * scale
        bond = rotate([0,0,0])(cylinder(r = (2.5 * scale),h = bond_len_base, segments = 50))
#       bond += rotate([90,0,90])(translate([-5,y_tot[counter],x_tot[counter]/2])(text(size = scale * 5, text = str(i_))))
# pin on bond is translated by -2*scale meaning that it gives the joint a unit length of 2
        bond += rotate([0,0,0])(translate([0,0,(-2 * scale)])(cylinder(r = (2 * scale) - 0.03,h = ((4 * scale) + bond_len_base),segments = 50)))
        bond += rotate([0, 0, 90])(translate([-2, -6 * scale, 0 - (2 * scale)])(linear_extrude(height = 1)(text(size=5, text=str(i_)+"-"+str(j_)))))
        bond_total += translate([x_tot[counter] + (2.5 * scale),y_tot[counter] + (2.5 * scale),(2 * scale)])(bond)
        counter += 1
        bond_progress.update(1)
if counter >= 6:
    bond_total += cube([x_tot[counter - 1] + (10 * scale), (60.0 * scale) + (5 * scale), 0.3])
else:
    bond_total += cube([x_tot[counter - 1] + (10 * scale),y_tot[counter] + (5 * scale),0.3])
#bond_total += cube([x_tot[counter] + (5 * scale), y_tot[counter] + (5 * scale),0.3])
scad_render_to_file(bond_total,'/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/bonds_to_print'+ str(time_date) +'.scad' )
subprocess.run(['/usr/bin/openscad', '-o', '/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/kit_export_models/'+'printing_files_'+str(input_filename)+str(time_date)+'/'+'bonds_to_print'+ str(time_date) +'.3mf', '/home/isaac/PycharmProjects/3D-atomic-visualiser/scad_files/bonds_to_print'+ str(time_date) +'.scad'])


