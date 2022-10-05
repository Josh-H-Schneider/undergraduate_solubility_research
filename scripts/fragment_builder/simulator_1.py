import os
import sys
import fileinput
import os.path
from os import path
import MDAnalysis as mda
import gromacs

frag_prefix = sys.argv[1]
receptor_prefix = sys.argv[2]
numol = int(sys.argv[3])

# NOTE the receptor top file must be copied into the combo_folder and named combo.top

fragitp = '"' + frag_prefix + '.itp"'
count = 0
project_frag_directory = f'./frag_box/{frag_prefix}/'
project_receptor_directory = './receptor_box/'
project_combo_directory = f'./combo_folder/{frag_prefix}/'

receptor =  os.path.join(project_receptor_directory, f'{receptor_prefix}.gro')
fragment =  os.path.join(project_frag_directory, f'{frag_prefix}.gro')
itpfile = os.path.join(project_combo_directory, f'{frag_prefix}.itp')
grocombo =  os.path.join(project_combo_directory, f'combo{frag_prefix}{numol}.gro')
topcombo =  os.path.join(project_combo_directory, f'combo{frag_prefix}{numol}.top')
emcombo =  os.path.join(project_combo_directory, f'combo_em_{frag_prefix}{numol}')
nvtcombo =  os.path.join(project_combo_directory, f'combo_nvt_{frag_prefix}{numol}')
mdcombo =  os.path.join(project_combo_directory, f'combo_md_{frag_prefix}{numol}')


file_gro_exists = path.exists(grocombo)

if file_gro_exists == True:
    pass
else:
    insert = gromacs.insert_molecules(f=receptor, ci=fragment, o=grocombo, nmol=f'{numol}', replace="12") #use 12 for rashad system

file_top_exists = path.exists(topcombo)

#####################################################################

# The below 6-7 lines are to prevent the files with different nmols from being named the same

if os.path.exists(os.path.join(project_combo_directory, 'combo.top')):
    os.rename(os.path.join(project_combo_directory, 'combo.top'), topcombo)

########################################################################



#The below top file itp inserter has tested on it's own and works
if file_top_exists == True:
    pass
else:
    for line in fileinput.FileInput(topcombo,inplace=1):
        if count < 1:
            if "#include" in line:
                line=line.replace(line,line+"\n; Include fragment topology\n#include " + fragitp + "\n")
                count += 1
        print (line, end='')


    with open(topcombo, 'a') as filedata:
        filedata.write(frag_prefix + "           "  + str(numol))
        filedata.close()

########################################################################
em_mdpfile = os.path.join(project_combo_directory, f'em_{frag_prefix}.mdp')
#tprfile = os.path.join(project_combo_directory, 'em.tpr')

file_tpr_exists = path.exists(f'{emcombo}.tpr')

if file_tpr_exists == True:
    pass
else:
    gromacs.grompp(f=em_mdpfile, c=grocombo, p=topcombo, o= f'{emcombo}.tpr', maxwarn="3")
    gromacs.mdrun(s= f'{emcombo}.tpr', v=True, deffnm= emcombo)

############################################################################


nvt_file = os.path.join(project_combo_directory, f'nvt_{frag_prefix}.mdp')
nvt_tpr_file = os.path.join(project_combo_directory, f'nvt_{frag_prefix}{numol}.tpr')

#npt_file = os.path.join(project_combo_directory, 'npt.mdp')
#npt_tpr_file = os.path.join(project_combo_directory, 'npt.tpr')

file_nvt_tpr_exists = path.exists(nvt_tpr_file)
#file_npt_tpr_exists = path.exists(npt_tpr_file)

if file_nvt_tpr_exists == True:
    pass
else:
    gromacs.grompp(f= nvt_file, c= emcombo, p= topcombo, o= nvt_tpr_file, maxwarn="3")
    gromacs.mdrun(s= nvt_tpr_file, deffnm= nvtcombo)

#if file_npt_tpr_exists == True:
    #pass
#else:
    #gromacs.grompp(f= npt_file, c="nvt.gro", p=topcombo, o=npt_tpr_file, maxwarn="3")
    #gromacs.mdrun(s= npt_tpr_file, deffnm="npt")




##########################################################################
md_mdpfile = os.path.join(project_combo_directory, f'md_{frag_prefix}.mdp')
md_tprfile = os.path.join(project_combo_directory, f'md_{frag_prefix}{numol}.tpr')

file_mdtpr_exists = path.exists(md_tprfile)

if file_mdtpr_exists == True:
    pass
else:
    gromacs.grompp(f= md_mdpfile, c=nvt_tpr_file, p= topcombo, o= md_tprfile)
    gromacs.mdrun(s= md_tprfile, v=True, deffnm= mdcombo)

###################################################################################
md_xtc_file = os.path.join(project_combo_directory, f'md_{frag_prefix}{numol}.xtc')
center_md = os.path.join(project_combo_directory, f'md_center_{frag_prefix}{numol}.xtc')
gromacs.trjconv(s= md_tprfile, f= md_xtc_file, o= center_md, center= "yes", pbc= "mol", ur= "compact", input=(f'{frag_prefix}', 'system'))
#chose protein or ligand for option 1 and for option 2 chose system
