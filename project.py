import openscad as o
import ase.build
from ase.neighborlist import neighbor_list
import numpy
import math
import test as t
from properties import bond_radii
from ase.data import covalent_radii


#a = ase.build.molecule('HCCl3')
#b = (a.get_positions())

a = t.molecule
b = (t.molecule.get_positions())
#a = ase.build.bulk("Si",cubic = True)
#a*=(2,2,2)
#a.set_pbc(False)
#b = a.get_positions()
#si.write("C:\\Users\\Isaac\\Documents\\Si.xyz")
species = numpy.unique(a.get_chemical_symbols())
cutoffs = {(str(s1),str(s2)):(bond_radii[s1]+bond_radii[s2])*1.1 for s1 in species for s2 in species}
print(cutoffs)
print(covalent_radii)
i,j,d,D = neighbor_list('ijdD',a,cutoffs)
#i,j,d,D = neighbor_list('ijdD',a,{("C","H"):1.2,("C","Cl"):1.9})
#print(i,j,d,D)
b *= 10
print(b)
for x in b:

    o.sphere(9)
    o.translate(*x)
    
#cylindrical coordinates, make sure x is zero
#z rotate begins in x axis

    
def z_rot():
    z_rot.rho = math.sqrt(u_[0]**2 + u_[1]**2)
    z_rot.angle = math.atan2(u_[1],u_[0])*(180/math.pi)
#    opposite = u_[1]/z_rot.rho
#    opposite = min(opposite,1.0)
#    opposite = max(opposite,-1.0)
#    if u_[0] >= 0:
#        z_rot.angle = math.asin(opposite)*(180/math.pi)
#    elif u_[0] < 0 and u_[1] >= 0:
#        z_rot.angle = (math.pi - math.asin(opposite))*(180/math.pi)
#    elif u_[0] < 0 and u_[1] < 0:
#        z_rot.angle = (math.pi - math.asin(opposite))*(180/math.pi)
    
def y_rot():
    y_rot.base = math.sqrt(u_[0]**2 + u_[1]**2)
    if u_[2] >= 0:
        y_rot.angle = 90 - (math.acos(y_rot.base)*(180/math.pi))
    elif u_[2] < 0:
        y_rot.angle = 90 + (math.acos(y_rot.base)*(180/math.pi))

u = D/d[:,None]
print(len(u))
counter = 0;
for i_,j_,d_,D_,u_ in zip(i,j,d,D,u):
    if i_ < j_: 
        print(i_,j_,d_,D_,u_)
        counter += 1
        o.cylinder(5,d_*10)
        z_rot()
        y_rot()
        o.rotate(0,y_rot.angle,z_rot.angle)
        o.translate(*b[i_])
    
    
    
    



    










o.output(o.result()) 
