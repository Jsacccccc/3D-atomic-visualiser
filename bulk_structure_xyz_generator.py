import ase.build

bulk_species = "Al"

a = ase.build.bulk(str(bulk_species), cubic=True)
a *= (2,2,2)
a.set_pbc(False)
a.write(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/real_test_files_xyz/" + str(bulk_species) + "_lattice.xyz")
