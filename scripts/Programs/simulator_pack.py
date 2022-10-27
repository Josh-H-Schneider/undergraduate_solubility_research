import os
import sys
import fileinput
import os.path
from os import path
import MDAnalysis as mda
import gromacs
#import gromacs.run
import shutil
import subprocess


#class MDrunnerMPI(gromacs.run.MDrunner):
   # """Manage running :program:`mdrun` as an MPI multiprocessor job."""
   # mdrun = "gmx_mpi mdrun"
   # mpiexec = "/opt/local/bin/mpiexec"

try:
    frag_prefix = sys.argv[1]  #prefix for the fragment
    receptor_prefix = sys.argv[2] #prefix for the receptor
    numol = int(sys.argv[3]) #number of copies of fragment in solution
    # NOTE the receptor top file must be copied into the combo_folder and named combo.top

    fragitp = '"' + frag_prefix + '.itp"'
    project_frag_itp = sys.argv [4]#itp file for the fragment
    project_top = sys.argv [5] #receptor top file
    pack_gro = sys.argv [6] #combined gro file
    project_combo_directory = sys.argv [7] #outputs
    mdp_loc = sys.argv [8] # mdp file locations
    amber_files = sys.argv [9] # file path for the amber force field folders
   # cores = sys.argv [10] # of cores being used
except:
    print ("")
    print ("")
    print( "There may have been an error with your inputs")
    print( "This program requires 8 total system inputs to function and they are listed below")
    print ("1: The prefix or name of the frament, i.e. MRP or ILE")
    print ("2: The prefix or name of the receptor for the systems, i.e. RAS or BASE")
    print ("3: The integer value for the amount of molecules of fragment contained within the system")
    print ("4: The file path for the itp file for the fragment")
    print ("5: The file path for the top file for the receptor of the system")
    print ("6: The file path for the gro file of the combined receptor and fragment system")
    print ("7: The file path for where you would like all of the outputs to go")
    print ("8: The file path to the location that contains all of MDP files needed for simulation")
    print ("9: The file path to the location that contains amber force field folder needed for simulation")
   # print ("10: The number of cores being used")
    print ("")
    print ("")

#mdrun = "gmx_mpi mdrun"

endfol =  os.path.join(project_combo_directory, f'{frag_prefix}{numol}')


#receptor =  os.path.join(project_receptor_directory, f'{receptor_prefix}.gro')
#fragment =  os.path.join(project_frag_directory, f'{frag_prefix}.gro')
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
em_mdpfile_orig = os.path.join(mdp_loc, 'em.mdp')
em_mdpfile = os.path.join(endfol, 'em.mdp')
tprfile = os.path.join(project_combo_directory, 'em.tpr')

nvt_file_orig = os.path.join(mdp_loc, 'nvt.mdp')
nvt_file = os.path.join(endfol, 'nvt.mdp')
nvt_tpr_file = os.path.join(endfol, f'nvt_{frag_prefix}{numol}.tpr')


npt_file_orig = os.path.join(mdp_loc, f'npt.mdp')
npt_file = os.path.join(endfol, f'npt.mdp')
npt_tpr_file = os.path.join(endfol, f'npt_{frag_prefix}{numol}.tpr')

md_mdpfile_orig = os.path.join(mdp_loc, f'md.mdp')
md_mdpfile = os.path.join(endfol, f'md.mdp')
md_tprfile = os.path.join(endfol, f'md_{frag_prefix}{numol}.tpr')
cptf = os.path.join(endfol, f'combo_md_{frag_prefix}{numol}.cpt')
file_mdtpr_exists = path.exists(md_tprfile)
####################################################################

destination_folder = f'{endfol}/amber14sb.ff'

file_endfol_exists = path.exists(endfol)

if file_endfol_exists == True:
    pass
else:
    os.mkdir(endfol)

    shutil.copytree(amber_files, destination_folder)

    with open(project_top) as f:
        with open(topcombo, 'w+') as f1:
            for line in f:
                f1.write(line)

    with open(pack_gro) as f2:
        with open(grocombo, 'w+') as f3:
            for line in f2:
                f3.write(line)
            
    with open(project_frag_itp) as f4:
        with open(itpfile, 'w+') as f5:
            for line in f4:
                f5.write(line)

    with open(nvt_file_orig) as f6:
        with open(nvt_file, 'w+') as f7:
            for line in f6:
                f7.write(line) 
           
    with open(em_mdpfile_orig) as f8:
        with open(em_mdpfile, 'w+') as f9:
            for line in f8:
                f9.write(line)
       
    with open(npt_file_orig) as f10:
        with open(npt_file, 'w+') as f11:
            for line in f10:
                f11.write(line)

    with open(md_mdpfile_orig) as f12:
        with open(md_mdpfile, 'w+') as f13:
            for line in f12:
                f13.write(line)
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

def topfix (topcombo, fragitp, frag_prefix, numol):
    count = 0
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
file_tpr_exists = path.exists(f'{emcombo}.tpr')

def em (em_mdpfile, grocombo, topcombo, emcombo):
    gromacs.grompp(f=em_mdpfile, c=grocombo, p=topcombo, o= f'{emcombo}.tpr', maxwarn="3")
    gromacs.mdrun(s= f'{emcombo}.tpr', v=True, deffnm= emcombo)

############################################################################



file_nvt_tpr_exists = path.exists(nvt_tpr_file)

def nvt (nvt_file, emcombogro, topcombo, nvt_tpr_file, nvtcombo):
    gromacs.grompp(f= nvt_file, c= emcombogro, p= topcombo, o= nvt_tpr_file, maxwarn="3")
    gromacs.mdrun(s= nvt_tpr_file, deffnm= nvtcombo)

#################################################################

file_npt_tpr_exists = path.exists(npt_tpr_file)

def npt (npt_file, nvtcombogro, topcombo, npt_tpr_file, nptcombo):
   gromacs.grompp(f= npt_file, c= nvtcombogro, p= topcombo, o= npt_tpr_file, maxwarn="3")
   gromacs.mdrun(s= npt_tpr_file, deffnm= nptcombo)


##########################################################################


def md (md_mdpfile, nptcombogro, topcombo, md_tprfile, mdcombo):
   gromacs.grompp(f= md_mdpfile, c= nptcombogro, p= topcombo, o= md_tprfile)
   # gromacs.mdrun(s= md_tprfile, cpi =cptf ,v=True, deffnm= mdcombo, append=True)
    #mdrun_mpi = MDrunnerMPI(s= md_tprfile, cpi =cptf ,v=True, deffnm= mdcombo, append=True)
    #mdrun_mpi.run(ncores=cores)   
   # subprocess.Popen(['mpirun', '-np', cores, 'gmx_mpi', 'mdrun',
   #                   '-s', md_tprfile, '-cpi', cptf, '-v', '-deffnm', mdcombo, '-append'],
   # :                stdin=subprocess.PIPE)
###################################################################################
#md_xtc_file = os.path.join(project_combo_directory, f'md_{frag_prefix}{numol}.xtc')
#center_md = os.path.join(project_combo_directory, f'md_center_{frag_prefix}{numol}.xtc')
#gromacs.trjconv(s= md_tprfile, f= md_xtc_file, o= center_md, center= "yes", pbc= "mol", ur= "compact", input=('MRP', 'system'))
#chose protein or ligand for option 1 and for option 2 chose system

#############################################################
#############################################################


topfix(topcombo, fragitp, frag_prefix, numol)

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
    
#if file_mdtpr_exists == True:
    #pass
#else:
md (md_mdpfile, nptcombogro, topcombo, md_tprfile, mdcombo)
