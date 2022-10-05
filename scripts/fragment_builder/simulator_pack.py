import os
import sys
import fileinput
import os.path
from os import path
import MDAnalysis as mda
import gromacs
import shutil

frag_prefix = sys.argv[1]  #prefix for the fragment
receptor_prefix = sys.argv[2] #prefix for the receptor
numol = int(sys.argv[3]) #number of copies of fragment in solution
# NOTE the receptor top file must be copied into the combo_folder and named combo.top

fragitp = '"' + frag_prefix + '.itp"'
count = 0
project_frag_itp = sys.argv [4]#itp file for the fragment
project_top = sys.argv [5] #receptor top file
pack_gro = sys.argv [6] #combined gro file
project_combo_directory = sys.argv [7] #outputs
mdp_loc = sys.argv [8] # mdp file locations

endfol =  os.path.join(project_combo_directory, f'{frag_prefix}{numol}')
os.mkdir(endfol)



receptor =  os.path.join(project_receptor_directory, f'{receptor_prefix}.gro')
fragment =  os.path.join(project_frag_directory, f'{frag_prefix}.gro')
itpfile = os.path.join(endfol, f'{frag_prefix}.itp')
grocombo =  os.path.join(endfol, f'combo{frag_prefix}{numol}.gro')
topcombo =  os.path.join(endfol, f'combo{receptor_prefix}{numol}.top')
emcombo =  os.path.join(endfol, f'combo_em_{frag_prefix}{numol}')
emcombogro =  os.path.join(endfol, f'combo_em_{frag_prefix}{numol}.gro')
nvtcombo =  os.path.join(endfol, f'combo_nvt_{frag_prefix}{numol}')
nvtcombogro =  os.path.join(endfol, f'combo_nvt_{frag_prefix}{numol}.gro')
nptcombo =  os.path.join(endfol, f'combo_npt_{frag_prefix}{numol}')
nptcombogro =  os.path.join(endfol, f'combo_npt_{frag_prefix}{numol}.gro')
mdcombo =  os.path.join(endfol, f'combo_md_{frag_prefix}{numol}')
#ndx =  os.path.join(project_frag_directory, f'{frag_prefix}.ndx')



file_gro_exists = path.exists(grocombo)
####################################################################


with open(project_top) as f:
    with open(topcombo) as f1:
        for line in f:
            f1.write(line)

with open(pack_gro) as f2:
    with open(grocombo) as f3:
        for line in f2:
            f3.write(line)
            
with open(project_frag_itp) as f4:
    with open(itpfile) as f5:
        for line in f4:
            f5.write(line)
            
         

#####################################################################

# The below 6-7 lines are to prevent the files with different nmols from being named the same

#if os.path.exists(os.path.join(project_combo_directory, 'combo.top')):
    #os.rename(os.path.join(project_combo_directory, 'combo.top'), topcombo)

#combo_topfile = os.path.join(project_combo_directory, 'combo.top')
#if os.path.exists(combo_topfile):
    #  Don't rename, make a copy!
    #shutil.copyfile(combo_topfile, topcombo)

########################################################################

#file_top_exists = path.exists(topcombo)



#The below top file itp inserter has tested on it's own and works

def topfix (topcombo, fragitp, frag_prefix):
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
em_mdpfile = os.path.join(mdp_loc, 'em.mdp')
#tprfile = os.path.join(project_combo_directory, 'em.tpr')

file_tpr_exists = path.exists(f'{emcombo}.tpr')

def em (em_mdpfile, grocombo, topcombo, emcombo):
    gromacs.grompp(f=em_mdpfile, c=grocombo, p=topcombo, o= f'{emcombo}.tpr', maxwarn="3")
    gromacs.mdrun(s= f'{emcombo}.tpr', v=True, deffnm= emcombo)

############################################################################


nvt_file = os.path.join(mdp_loc, 'nvt.mdp')
nvt_tpr_file = os.path.join(endfol, f'nvt_{frag_prefix}{numol}.tpr')


file_nvt_tpr_exists = path.exists(nvt_tpr_file)

def nvt (nvt_file, emcombogro, topcombo, nvt_tpr_file, nvtcombo):
    gromacs.grompp(f= nvt_file, c= emcombogro, p= topcombo, o= nvt_tpr_file, maxwarn="3")
    gromacs.mdrun(s= nvt_tpr_file, deffnm= nvtcombo)

#######################################################################


npt_file = os.path.join(mdp_loc, f'npt_{frag_prefix}.mdp')
npt_tpr_file = os.path.join(endloc, f'npt_{frag_prefix}{numol}.tpr')

file_npt_tpr_exists = path.exists(npt_tpr_file)

def npt (npt_file, nvtcombogro, topcombo, npt_tpr_file, nptcombo):
   gromacs.grompp(f= npt_file, c= nvtcombogro, p= topcombo, o= npt_tpr_file, maxwarn="3")
   gromacs.mdrun(s= npt_tpr_file, deffnm= nptcombo)


##########################################################################
md_mdpfile = os.path.join(mdp_loc, f'md_{frag_prefix}.mdp')
md_tprfile = os.path.join(endfol, f'md_{frag_prefix}{numol}.tpr')

file_mdtpr_exists = path.exists(md_tprfile)

def md (md_mdpfile, nptcombogro, topcombo, md_tprfile, mdcombo):
    gromacs.grompp(f= md_mdpfile, c= nptcombogro, p= topcombo, o= md_tprfile)
    gromacs.mdrun(s= md_tprfile, v=True, deffnm= mdcombo)

###################################################################################
#md_xtc_file = os.path.join(project_combo_directory, f'md_{frag_prefix}{numol}.xtc')
#center_md = os.path.join(project_combo_directory, f'md_center_{frag_prefix}{numol}.xtc')
#gromacs.trjconv(s= md_tprfile, f= md_xtc_file, o= center_md, center= "yes", pbc= "mol", ur= "compact", input=('MRP', 'system'))
#chose protein or ligand for option 1 and for option 2 chose system

#############################################################
#############################################################


topfix(topcombo, fragitp, frag_prefix)

if file_tpr_exists == True:
    pass
else:
    em(em_mdpfile, grocombo, topcombo, emcombo)

if file_nvt_tpr_exists == True:
    pass
else:
    nvt(nvt_file, emcombogro, topcombo, nvt_tpr_file, nvtcombo)
    
if file_npt_tpr_exists == True:
    pass
else:
    npt (npt_file, nvtcombogro, topcombo, npt_tpr_file, nptcombo)
    
if file_mdtpr_exists == True:
    pass
else:
    md (md_mdpfile, nptcombogro, topcombo, md_tprfile, mdcombo)