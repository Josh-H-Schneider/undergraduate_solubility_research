This read_me is for the simulator.py program and will explain how to use it.

run python simulator_pack.py and after a few seconds a prompt will appear and tell you what to do


below is outdated

# This read_me file contains the information required to succesfully perform pt 2 of the fragment simulation started in part 1. 
# Part 2 will combine the fragment and the receptor into one gro and top file, as well as adding the itp file into the top file.
# After this it minimize the energy of the system, equilibrate the system, simulate the system, and then center the simulated system so it can be correctly viewed.


# Steps:
# Firstoff, you will need the top, gro, and itp files from the make_fragment_topologies.py program.
# then in the frag_box folder make another folder that is just the name of the fragment and make another folder named the same in the combo_folder folder.
# then put all the previously made fragment files except the itp file in the frag_box fragment folder and place the itp file into the fragment folder in combo_folder.
# Take the gro file for the receptor, name it what you want and just place it into the receptor_box.
# take the itp files for the protein and place them in the fragment folder in the combo_folder folder.
# take the top file for the receptor system and place it in the fragment folder in combo_folder and name it combo.top ***Make copy of it as it will change directly***
# then make a make a copy of all the loose files and the tip3 forcefield folder (If you are using it) and place the copies into the fragment folder in the combo_folder.
# Name all the newly placed files according to the conventions listed in the files in the MRP folder in combo_folder. 
# ***Use these files as a reference if you have questions as they are all properly named and formatted***
# You will need to change this section in the NVT file:

# tcoupl                  = V-rescale                     ; modified Berendsen thermostat
# tc-grps                 = Protein WAA Water_and_ions    ; two coupling groups - more accurate
# tau_t                   = 0.1   0.1  0.1                 ; time constant, in ps
# ref_t                   = 300   300  300                  ; reference temperature, one for each group, in K
# ; Pressure coupling

# then fix it so the WAA is the name of your fragment and copy and paste lines 21-23 from this fixed file to the md.mdp file.

# After this all steps are prepped and you can move on to running the command.

# The command to run will be from the fragment_builder directory and is:
# simulator.py <Fragment name> <Receptor name> <number of copies of fragment you want in solution>.
                                                     
                                                    
# This command should work UP UNTIL the energy minimization step at which point it will error out. 
# The error message will say that the atoms numbers in the gro and top files do not match and it will give you a number for the atoms counts of both the gro and top.
# subtract the gro # from the top # and divide it by 3 (Only works if your solis water if else, divide by the number of atoms in one SOL molecule). 
# after this if the number is positive you will subtract that from the number next to SOL at the bottom of the .TOP document ***DO NOT TOUCH ANYTHING ELSE***.
# if it is negative add.
                                                     
# After all these changes have been made run the program again with the ***EXACT SAME*** command.

# The program should now run until completion and spit out all the correctly named docs into the fragment folder in the combo_folder folder.
                                                     
