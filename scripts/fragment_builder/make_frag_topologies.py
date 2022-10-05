import os
import sys
import mdtraj as md
sys.path.append('/media/tun37047/HD/tun37047/sampl9-voelzlab/scripts/fragment_builder')

from Frag_Builder_System import Fragment_Build
from pathlib import Path
import shutil


# Functions

def convert_topfile(topfile, itpfile=None):
    """Reads in a GROMACS *.top file and converts the uppercase atomtypes (e.g. 'C1') to
    lower case ('c1') atomtypes, to avoid conflict with AMBER ff atomtypes.

    INPUTS
        topfile

    PARAMETERS
        itpfile        If specified, a corresponding *.itp file will be written by removing
                       the header and footer.

    WARNING            This function will overwrite the topfile.
    """

    # Parse the directives of the topfile
    pathcorrected_top = Path(topfile)
    topology = open(pathcorrected_top, 'r')

    lines = topology.readlines()
    topology.close()

    def is_directive(line):
        fields = (line.strip()).split(' ')
        if (fields[0] == '[') and (fields[2] == ']'):
            return True
        else:
            return False

    if lines[-1] != '\n':
        lines.append('\n')

    header = [ lines.pop(0) ]
    directives = {}
    keys_in_order = []

    i = 0
    while i < len(lines):

        if is_directive(lines[i]):
            key = lines[i]
            keys_in_order.append(key)
            directives[ lines[i] ] = []
            i +=1
            while lines[i] != '\n':
                directives[ key ].append( lines[i] )
                i +=1
            directives[ key ].append( lines[i] )
            i +=1

        else:
            header.append( lines[i] )
            i +=1

    # Convert uppercase atomtypes to lower in the '[ atomtypes ]' directive
    for key in keys_in_order:

        if key.count('atomtypes') == 1:
            for i in range(len(directives[key])):
                if directives[key][i][0] != ';':
                    directives[key][i] = directives[key][i].lower()

        if key.count('atoms') == 1:
            for i in range(len(directives[key])):
                if directives[key][i][0] != ';':
                    # Assume that the atom name occurs in the first 22 characters
                    directives[key][i] = directives[key][i][0:22].lower() + directives[key][i][22:]



    with open(pathcorrected_top, 'w') as new_top:
        for line in header:
            new_top.write(line)
        for key in keys_in_order:
            new_top.write(key)
            for line in directives[ key ]:
                new_top.write(line)
    new_top.close()

    # If the itpfile is specified, write a corresponing *.itp file with the header and foot removed 
    if itpfile != None:
        pathcorrected_itp = Path(itpfile)

        with open(pathcorrected_itp, 'w') as new_itp:
            for line in header:
                new_itp.write(line)
            for key in keys_in_order:
                if key.count('defaults') + key.count('atomtypes') + key.count('system') + key.count('molecules') == 0:
                    new_itp.write(key)
                    for line in directives[ key ]:
                        new_itp.write(line)

        new_itp.close()

def fix_unl(grofile, topfile, itpfile, lig_prefix):
    """This will get rid of UNL strings in these files."""

    # Get rid of UNL strings in the *.gro file
    with open(grofile, 'r') as file :
        filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('UNL', lig_prefix)
        filedata = filedata.replace('UNK', lig_prefix)
        with open(grofile, 'w') as file:
            file.write(filedata)


    # Get rid of UNL strings in the *.top file    
    with open(topfile, 'r') as file :
        filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('UNL', lig_prefix)
        filedata = filedata.replace('UNK', lig_prefix)
        # Write the file out again
        with open(topfile, 'w') as file:
            file.write(filedata)

    # Get rid of UNL strings in the *.itp file  
    with open(itpfile, 'r') as file :
        filedata = file.read()
        # Replace the target string
        filedata = filedata.replace('UNL', lig_prefix)
        filedata = filedata.replace('UNK', lig_prefix)
        # Write the file out again
        with open(itpfile, 'w') as file:
            file.write(filedata)

#def fix_naming (topfile, grofile, built_topfile, built_grofile) :
    #os.rename(built_topfile, topfile)
    #os.rename(built_grofile, grofile)


# Main

if __name__ == "__main__":

    usage = """
Usage:   python make_frag_topologies.py LIGPREFIX [fragmentdir]

    This script will create GROMACS *.gro, *.top and *.itp files for small molecules from 

    INPUTS 
        LIGPREFIX       This prefix (e.g. WAA) is the naming convention for *.sdf and/or *.pdb files
                        (e.g. WAA.sdf and WAA.pdb) found in the fragment directory (Default: ./frag_box)

        fragmentdir     OPTIONAL.   Use this as the fragment directory instead of ./frag_box

    OUTPUTS

        Will create these files in the fragment directory: LIGPREFIX.gro, LIGPREFIX.top, and LIGPREFIX.itp

    Examples:
        $ python make_frag_topologies.py AAB
        or
        $ python make_frag_topologies.py AAB my_fragments

    Help:        python make_frag_topologies.py -h 

"""

    if len(sys.argv) < 2:
        print(usage)
        sys.exit(1)

    lig_prefix = sys.argv[1]

    if lig_prefix == '-h':
        print(usage)
        sys.exit(1)

    if len(sys.argv) >= 2:
        project_directory = sys.argv[2]
    else:
        project_directory = './frag_box/'     

 
    s = Fragment_Build(project_directory, ionic_strength=0, temperature=298.15,
            small_molecule_forcefield='openff-2.0.0', hydrogen_mass=1, timestep=2)
            
    s.build_fragment(ligand_prefix=f'{project_directory}/{lig_prefix}', small_mol_rec=False)
    
    #built_grofile, built_topfile = s.build_fragment(ligand_prefix=f'{project_directory}/{lig_prefix}', small_mol_rec=False)

    os.rename( './frag_box/none.gro', f'./frag_box/{lig_prefix}.gro')
    #os.replace(f'./frag_box/{lig_prefix}.gro', f'{project_directory}.gro')
    shutil.move(f'./frag_box/{lig_prefix}.gro', f'{project_directory}')

    os.rename( './frag_box/none.top', f'{project_directory}.top')
    #os.replace(f'./frag_box/{lig_prefix}.top', f'{project_directory}.top')
    shutil.move(f'./frag_box/{lig_prefix}.top', f'{project_directory}')

    grofile = os.path.join(project_directory, f'{lig_prefix}.gro')
    topfile = os.path.join(project_directory, f'{lig_prefix}.top')
   
 

    # The Fragment_Build() creates files None.* ???
    #built_grofile = os.path.join(project_directory, 'None.gro')
   # built_topfile = os.path.join(project_directory, 'None.top')
    
    #fix_naming(topfile, grofile, built_topfile, built_grofile)

    
    # Convert the topology file to lowercase atomtypes and write an itp file
    itpfile = os.path.join(project_directory, f'{lig_prefix}.itp')
    convert_topfile(topfile, itpfile=itpfile)
    # convert_to_itp(topfile, itpfile)
    
    # change the naming residue UNL to the lig_prefix, if needed
    
    fix_unl(grofile, topfile, itpfile, lig_prefix)


     

