3D atomic visualiser program,  2025-02-19

System Requirements:
-------------------
- Linux required due to ColorSCAD dependency
- 16 GB ram or more for rendering structures of more than ~50 atoms


Installation Requirements: READ CAREFULLY
-------------------------
Install modules in using python console:
Type: pip install <module>

- ASE
- solidpython
- numpy
- tqdm

Need to install OpenSCAD from here (with guide): https://openscad.org/downloads.html
Need to install ColorSCAD from here (with guide): https://github.com/jschobben/colorscad

Some directories also need to be updated:
In project.py: IMPORTANT
-------------
- Use ctrl+f to find "/home/isaac/PycharmProjects/3D-atomic-visualiser/clean_XYZs/" and replace /home/isaac/PycharmProjects with the directory you have saved the program to.

- Go to the last line of the script which reads: "subprocess.run(["/home/isaac/Code/colorscad/colorscad.sh", "-i", "/home/isaac/PycharmProjects/3D-atomic-visualiser/project.scad", "-o", "/home/isaac/PycharmProjects/3D-atomic-visualiser/export_models/single_body_export_files/" + str(output_filename) + str(time_date) + ".3mf"])"
- Change the first directory to the directory that you have installed ColorSCAD to.
- For the last two directories, change the /home/isaac/PycharmProjects to the directory that you are running the program from.


In kitgenerator.py: IMPORTANT
------------------
- Use ctrl+f to find every instance of "/home/isaac/PycharmProjects" and replace this with the directory of the program on your system. There should be 7 instances to replace.
- Go to the last line and replace "/usr/bin/openscad" and replace this with the directory where you have installed OpenSCAD on your system.




General Usage:
-------------

For the project.py script, the arguments are:
--input       for the name of your input file including the file extension (xyz, pdb etc).

--output      for specifying the name of the 3mf output file, EXCLUDING the file extension.

-cutoff       this specifies the multiplier for the cutoff of the bonds and should only be
              changed from 1.2 if bonds aren't generating when they should be.

-scale        This specifies the scale factor, where 1 is default (used over slicer scaling
              so that numbers and raft thickness remain constant).

Example command: python project.py -scale 1.0 -cutoff 1.2 --input hello.xyz --output goodbye

For the kitgenerator.py script, the arguments are the same apart from the --output command is not used, since the number of output files varies depending on how many species are present in the chemical structure. The atom files will be named "atoms_to_print" followed by the date, time and number in the structure and the bonds will be named "bonds_to_print" and both of these will be found in the "export_models/kit_export_models" folder in the program files.


Colouring structures from the project.py script:

1. Firstly, generate a model using the project.py script and open the .3mf file generated using Bambu Studio (or similar slicer program)
2. On the left hand side, change the global selection to object. This shows each species as a separate body.
3. Add filaments using the + button and click on the coloured box next to each body to assign a filament to it.
4. Slice and print your model.

Printing Structures from the kitgenerator.py script:

1. Firstly, generate your model using the command line with the kitgenerator.py script and open the folder "/export_models/kit_export_models" as your exported models will be here.
2. Open all of the atoms_to_print and bonds_to_print files with the current date and time at the end of the filename in separate instances of Bambu slicer.
3. Print each set of atoms/bonds in different filaments but at standard printer settings.
4. Once printed, note the identifying numbers on the plates for the atoms/bonds.
5. Carefully remove atoms one at a time, remove their corresponding bonds and join them together - they should have a good friction fit if printed correctly.
