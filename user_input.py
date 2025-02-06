import ase.build
import ase.io
import test as t

def asemolecule(self):
    a = ase.build.molecule(self)
    return a
def databasemolecule():
    a = molecule
    return a
def bulkgen():
    a = ase.build.bulk("C",cubic = True)
    a*=(2,2,2)
    a.set_pbc(False)
#   si.write("C:\\Users\\Isaac\\Documents\\Si.xyz")
    return a
def structure_from_file():
    a = real_molecule
    a.set_pbc(False)
    return a

#This reads the custom file that the user has been instructed to place into this folder. Since the names are not pre-determined, this is the only part that the user will have to change manually
real_molecule = ase.io.read(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/real_test_files_xyz/He_cluster_in_W.xyz")

#This part prints the following text and prompts the user to pick which script to run
print("Please choose a mode to generate your model:\n1. One piece model - this model generates as one body and is recolourable in a slicer"
      "\n2. Kit generator - this model generates an assemblable kit where species are separated into separate recolourable files")
model_selection = input("Please type either 1 or 2 to select your option: ")

#This checks that the user has entered a valid response and exits the script if not with a brief description why
if float(model_selection) != 1.0 and 2.0:
    exit("invalid response >:(")

#This part prompts the user to pick which source to generate their model from
print("\nPlease choose a source to generate your model from:\n1. A model from the built in ASE database"
      "\n2. A model from the built in database containing 133,885 models"
      "\n3. A bulk structure"
      "\n4. A custom user submitted file (check readme for more information)")
source_selection = input("Please choose a number between 1 and 4: ")

#This if statement list determines what to do for each input and checks that the user has entered a valid input, terminating the script if not with a brief description why
if float(source_selection) == 1.0:
    pick_molecule = input("Write the formula of the molecule you would like to model (use all caps, e.g. CH4): ")
    a = asemolecule(pick_molecule)

elif float(source_selection) == 2.0:
    pick_database = input("Which file would you like to choose from the database? \nPick a number between 1 and 133,885: ")
    molecule = ase.io.read(r"/home/isaac/PycharmProjects/3D-atomic-visualiser/clean_XYZs/d"+str(pick_database)+".xyz", format="xyz")
    a = databasemolecule()

elif float(source_selection) == 3.0:
    a = bulkgen()

elif float(source_selection) == 4.0:
    a = structure_from_file()
else:
    exit("invalid response >:(")

#this first if statement opens, reads and executes the code from the project.py file if option 1 is selected in the first choice
if float(model_selection) == 1.0:
    with open("project.py") as file:
        exec(file.read())
#this second elif statement opens, reads and executes the kitgenerator code if option 2 is selected first
elif float(model_selection) == 2.0:
    with open("kitgenerator.py") as file:
        exec(file.read())
#both of these options use the a variable calculated using the user's choices