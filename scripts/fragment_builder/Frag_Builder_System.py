from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from simtk import openmm, unit
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
from openbabel import openbabel
import parmed
import os, re, sys, subprocess
import numpy as np
global GMX_BIN
try:
    GMX_BIN = os.environ['GMXBIN']   
except:
    GMX_BIN = '/usr/local/gromacs/bin'
if not os.path.exists(GMX_BIN):
    print('Cannot find the GROMACS installation!  Please set environment variable GMX_BIN')
    sys.exit(1)

usage = """Usage:   python make_frag_topologies.py LIGPREFIX

    The *frag prefix* should be the same name for both the sdf and pdb files that are in the './frag_box' directory  and an example of how it should be is 
    AAB.pdb and AAB.sdf are both in the frag box so the frag prefix put after the python command will just be AAB.

    Help:        python make_frag_topologies.py -h 

    Example:     python make_frag_topologies.py AAB  
"""

class Fragment_Build:  
    def __init__(self, project_directory, box_size= (4.5,4.5,4.5), solvent_padding= None, ionic_strength= 100,
            pressure= 1.0, collision_rate= 91.0/unit.picoseconds, temperature= 300.0,
            timestep= 2.0, nsteps_equil= 5000, protein_forcefield= 'amber14/protein.ff14SB.xml',
            small_molecule_forcefield = 'openff-2.0.0', solvation_forcefield = 'amber14/tip3p.xml',
            write_gromacs_files= True, write_openmm_files= False, hydrogen_mass= 1.01,
            removeCMMotion= False, rigidWater= False, minimize= False, equilibrate= False):
    
        """ Prepare a receptor only, ligand only, or receptor-ligand system for an expanded ensamble
            simultion using openMM
            Parameters
            ----------
            project_directory: string
                directory where receptor and ligand files can be found
            box_size: three number tuple
                box dimentions in Nanometer
                box size will supersede solvent padding
                default is 4.5 nm
            solvent_padding: float or unit.quantity()
                distance between bio-molecule and edge of box (nm)
                WARNING solvant padding is overwritten by box_size
                will be passed to modeller.addSolvent()
                default = 0.35 nm
            ionic_strength: float or unit.quantity()
                Concentration of NaCL ions (millimolar)
                will be passed to modeller.addSolvent()
                default = 100 millimolar
            pressure: float or unit.quantity()
                pressure for NPT simulaton (atmospheres)
                will be passed to openmm.MonteCarloBarostat()
                default = 1 atmosphere
            collision_rate: float or unit.quantity()
               collision rate for Langevin thermostat (1/ps)
               will be passed to openmm.LangevinIntegrator()
               default = 91 / picoseconds
            temperature: float or unit.quantity()
               temperature for simulation (kelvin)
               will be passed to openmm.MonteCarloBarostat()
               default = 300 kelvin
            timestep: float or unit.quantity()
               time between force calculations (fs)
               will be passed to openmm.LangevinIntegrator()
               default 4 fs
            nsteps_equil: int
               number of simulation steps taken for equilibration
               default = 5000
            protein_forcefield: openmm .xml file
               protein forcefield to be used
               default = amber14/protein.ff14SB.xml
            small_molecule_forcefield: openmm .xml file
               small molecule forcefield to be used
               default = openff-1.3.0
               another option: gaff-2.11
            solvation_force_field: openmm .xml file
               solvation forcefield to be used
               default = amber14/tip3p.xml
            write_gromacs_files: Boolean
               write output as Gromacs files
               default = True
            write_openmm_files: Boolean
               write output as OpenMM files
               default = False
            hydrogen_mass: float or unit.quantity()
               use heavy hydrogen to increase timestep
               default = 4
            removeCMMotion: Boolean
               remove center of mass motion
               default = False
            rigidWater: Boolean
               use rigid or flexable water
               default = False
        """ # It turns out that box_size cannot be formatted this way. A plain-old tuple works fine though.
        #if unit.is_quantity(box_size):
        #    box_size = box_size.value_in_unit(unit.nanometer)
        #elif box_size is not None:
        #    box_size = box_size*unit.nanometers
        #self.box_size = openmm.vec3.Vec3(box_size[0],box_size[1],box_size[2])
        #self.box_size = box_size
        
        if box_size is not None:
            self.box_size = box_size
        elif solvent_padding is not None:
            if unit.is_quantity(solvent_padding):
                solvent_padding = solvent_padding.value_in_unit(unit.nanometer)
            else:
                solvent_padding = solvent_padding*unit.nanometer
            self.solvent_padding = solvent_padding
        else:
            print(self.box_size)
            print('Warning: you must select either a box size or padding')

        if unit.is_quantity(ionic_strength):
            ionic_strength = ionic_strength.value_in_unit(unit.millimolar)
        else:
            ionic_strength = ionic_strength * unit.millimolar
        self.ionic_strength = ionic_strength
        if unit.is_quantity(pressure):
            pressure.value_in_unit(unit.atmospheres)
        else:
            pressure = pressure * unit.atmospheres
        self.pressure = pressure
        ## TO DO: allow unit conversion of collision rate
        ## right now collosion rate must be given in picoseconds^-1
        if not unit.is_quantity(collision_rate):
            collision_rate = collision_rate / unit.picoseconds
        self.collision_rate = collision_rate
        if unit.is_quantity(temperature):
            temperature.value_in_unit(u.kelvin)
        else:
            temperature = temperature * unit.kelvin
        self.temperature = temperature
        if unit.is_quantity(timestep):
            timestep.value_in_unit(femtoseconds)
        else:
            timestep = timestep * unit.femtoseconds
        self.timestep = timestep
        
        self.nsteps_equil = nsteps_equil
        self.protein_forcefield = protein_forcefield # 'amber14/protein.ff14SB.xml'
        self.small_molecule_forcefield = small_molecule_forcefield # 'openff-1.3.0' or 'gaff-2.11'
        self.solvation_forcefield = solvation_forcefield # 'amber14/tip3p.xml'
        self.write_gromacs_files = write_gromacs_files
        self.write_openmm_files = write_openmm_files
        if unit.is_quantity(hydrogen_mass):
            hydrogen_mass.value_in_unit(unit.amu)
        else:
            hydrogen_mass = hydrogen_mass * unit.amu
        self.hydrogen_mass = hydrogen_mass
        self.removeCMMotion = False
        self.rigidWater = False
        self.project_directory = './frag_box/'
        self.minimize = minimize
        self.equilibrate = equilibrate

    def _convert_mol2_to_sdf(self, structure_prefix):
        for output in ['pdb','sdf']:
            obConversion = openbabel.OBConversion()
            obConversion.SetInAndOutFormats("mol2", output)
            mol = openbabel.OBMol()
            obConversion.ReadFile(mol, f'{structure_prefix}.mol2')
            outMDL = obConversion.WriteFile(mol,f'{self.project_directory}/LIG_temp.{output}')
        ligand_structure = parmed.load_file(f'{self.project_directory}/LIG_temp.pdb')
        ligand = Molecule.from_file(f'{self.project_directory}/LIG_temp.sdf')
        os.remove(f'{self.project_directory}/LIG_temp.sdf')
        os.remove(f'{self.project_directory}/LIG_temp.pdb')
        return ligand_structure, ligand

    def _make_sampl_index_files(self, conf):
        ''' Make gourps to pipe to gmx make_ndx for the
        Sampl9 host-guest challenge
        inputs:
            conf: string
            path to gro, pdb, or tpr file
        outputs:
            list of index names
            writes groups to file
                index_options.txt
                index_strings.txt
        Warning: these options will require a tpr with make_ndx
        example:
        cat index_options.txt | gmx make_ndx -f host_lig7.tpr
        '''

        import mdtraj as md

        deprotonated_edges = [0,5,
                            25, 30,
                            50, 55,
                            75, 80,
                            100, 105,
                            125, 130]

        top = md.load(conf).top
        true_edges = top.select('resname == WP6')[deprotonated_edges]
        new_list = np.array(true_edges).reshape(6,2)
        s_out = [f'a {i[0] +1}-{i[1] +1}' for i in new_list]
        f_out = ''
        for i in range(len(s_out)):
            f_out += s_out[i]
            if i <= 4:
                f_out += '|'
            else:
                break

        with open(f'{self.project_directory}/index_options.txt', 'w') as f:
            f.write('2 | 3\n')
            f.write(f_out + '\n')
            f.write('2 & ! t H*\n')
            f.write('3 & ! t H*\n')
            f.write('quit\n')
            f.close()

        f_out.replace(' ', '_')
        f_out.replace('|', '-')

        with open(f'{self.project_directory}/index_string.txt', 'w') as f1:
            f1.write('LIG_WP6\n')
            f1.write(f_out + '\n')
            f1.write('LIG_&_!H*\n')
            f1.write('WP6_&_!H*\n')
            f1.close()
        return ['LIG_WP6', f_out, 'LIG_&_!H*', 'WP6_&_!H*']

    def _make_sampl_index_files_8_10(self, conf):
        ''' Make gourps to pipe to gmx make_ndx for the
        Sampl9 host-guest challenge
        inputs:
            conf: string
            path to gro, pdb, or tpr file
        outputs:
            list of index names
            writes groups to file
                index_options.txt
                index_strings.txt
        Warning: these options will require a tpr with make_ndx
        example:
        cat index_options.txt | gmx make_ndx -f host_lig7.tpr
        '''

        import mdtraj as md

        deprotonated_edges = [0,35]

        top = md.load(conf).top
        true_edges = top.select('resname == WP6')[deprotonated_edges]
        f_out = f'a {true_edges[0] +1}-{true_edges[1] +1}'

        with open(f'{self.project_directory}/index_options.txt', 'w') as f:
            f.write('2 | 3\n')
            f.write(f_out + '\n')
            f.write('2 & ! t H*\n')
            f.write('3 & ! t H*\n')
            f.write('quit\n')
            f.close()

        f_out.replace(' ', '_')
        f_out.replace('|', '-')

        with open(f'{self.project_directory}/index_string.txt', 'w') as f1:
            f1.write('LIG_WP6\n')
            f1.write(f_out + '\n')
            f1.write('LIG_&_!H*\n')
            f1.write('WP6_&_!H*\n')
            f1.close()
        return ['LIG_WP6', f_out, 'LIG_&_!H*', 'WP6_&_!H*']


    def insert_posre_to_top(self, top):
        ''' inserts an if statment in the middle of a Gromacs topology file
            to include position restraint itp file'''

        molecule_line_list = []
        with open(top, 'r') as f3:
            for i, line in enumerate(f3):
                if '[ moleculetype ]' in line:
                    molecule_line_list.append(i)

        with open(top, 'r') as f4: 
            contents = f4.readlines() 

        contents.insert(molecule_line_list[2], '#ifdef POSRE\n#include "WP6_posre.itp"\n#endif\n\n')
        contents.insert(molecule_line_list[1], '#ifdef POSRE\n#include "LIG_posre.itp"\n#endif\n\n')

        with open(top, 'w') as f5: 
            contents = "".join(contents) 
            f5.write(contents)

    def write_gromacs_min_mdp(self, filename, steps=100000):
        ''' writes a mdp file for minimization in gromacs'''
        with open(filename, 'w') as f2:
            f2.write(f'''integrator          = steep
nsteps                 = {steps}
emtol                  = 100
; Bond parameters
continuation           = no                          ; Initial simulation
constraint_algorithm   = lincs                       ; holonomic constraints
constraints            = h-bonds                     ; all bonds (even heavy atom-H bonds) constrained
lincs_iter             = 1                           ; accuracy of LINCS
lincs_order            = 4                           ; also related to accuracy
; Neighborsearching
ns_type                = grid                        ; search neighboring grid cels
nstlist                = 10                          ; 20 fs
rlist                  = 1.4                         ; short-range neighborlist cutoff (in nm)
rcoulomb               = 1.4                         ; short-range electrostatic cutoff (in nm)
rvdw                   = 1.4                         ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype            = PME                         ; Particle Mesh Ewald for long-range electrostatics
pme_order              = 4                           ; cubic interpolation
fourierspacing         = 0.16    ''')
            f2.close()

    def write_gromacs_equil_mdp(self, filename):
        ''' writes a gromacs mdp file for equilibration'''
        with open(filename, 'w') as f:
            f.write(f'''define                 = -DPOSRE          ; position restrain the protein
; Run parameters
integrator             = md                          ; leap-frog integrator
nsteps                 = 100000                      ; 2 * 100000 = 200 ps
dt                     = 0.002                       ; 2 fs
; Output control
nstxout                = 1000                        ; save coordinates every 2 ps
nstvout                = 1000                        ; save velocities every 2 ps
nstenergy              = 1000                        ; save energies every 2 ps
nstlog                 = 1000                        ; update log file every 2 ps
; Bond parameters
continuation           = no                          ; Initial simulation
constraint_algorithm   = lincs                       ; holonomic constraints
constraints            = h-bonds                     ; all bonds (even heavy atom-H bonds) constrained
lincs_iter             = 1                           ; accuracy of LINCS
lincs_order            = 4                           ; also related to accuracy
; Neighborsearching
ns_type                = grid                        ; search neighboring grid cels
nstlist                = 10                          ; 20 fs
rlist                  = 1.4                         ; short-range neighborlist cutoff (in nm)
rcoulomb               = 1.4                         ; short-range electrostatic cutoff (in nm)
rvdw                   = 1.4                         ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype            = PME                         ; Particle Mesh Ewald for long-range electrostatics
pme_order              = 4                           ; cubic interpolation
fourierspacing         = 0.16                        ; grid spacing for FFT
; Temperature coupling is on
tcoupl                 = V-rescale                   ; Weak coupling for initial equilibration
tc-grps                = non-Water   Water           ; two coupling groups - more accurate
tau_t                  = 0.1       0.1               ; time constant, in ps
ref_t                  = 298       298               ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                 = Berendsen                   ; Pressure coupling on in NPT, also weak coupling
pcoupltype             = isotropic                   ; uniform scaling of x-y-z box vectors
tau_p                  = 2.0                         ; time constant, in ps
ref_p                  = 1.0                         ; reference pressure (in bar)
compressibility        = 4.5e-5                      ; isothermal compressibility, bar^-1
refcoord_scaling       = com
; Periodic boundary conditions
pbc                    = xyz                         ; 3-D PBC
; Dispersion correction
DispCorr               = EnerPres                    ; account for cut-off vdW scheme
; Velocity generationd
gen_vel                = yes                         ; Velocity generation is on
gen_temp               = 310                         ; temperature for velocity generation
gen_seed               = -1                          ; random seed
; These options remove COM motion of the system
nstcomm                = 10
comm-mode              = Linear
comm-grps              = System''')
        f.close()

    def build_fragment(self, receptor_prefix=None, ligand_prefix=None, ligand_pre=None, convert_mol2_to_sdf_ligand=False, small_mol_rec=False, convert_mol2_to_sdf_receptor=False):
        """ Prepare system for simulation
        will prepare receptor only, ligand only, or receptor-ligand complex
        Inputs
        receptor_prefix: string
            prefix of receptor pdb
            or if it's a small molecule (small_mol_rec=Ture) pdb and sdf
        ligand_prefix: string
            prefix of ligand pdb and sdf or mol2
        convett_mol2_to_sdf_ligand: Boolean
            supply sdf or used openbabel to convert mol2 to pdb and sdf
        small_mol_rec: Boolean
            receptor is a non-protein and needs to be paramaterized as a small molecule
        convert_mol2_to_sdf_receptor: Boolean
            receptor is a small molecule and needs to be converted from mol2 to pdb and sdf
        Outputs
        Files to run a Gromacs and/or OpenMM simulaton
        """
        if receptor_prefix:
            if convert_mol2_to_sdf_receptor==True:
                receptor_structure, receptor = self._convert_mol2_to_sdf(receptor_prefix)
            elif small_mol_rec==True:
                try: 
                    receptor_structure = parmed.load_file(f'{receptor_prefix}.pdb')
                except:
                    print ('There must be a properly named fragment.pdb and fragment.sdf file within the frag_box folder or other fragment box recepticle folder')
                    sys.exit(1)
                try:                
                    receptor = Molecule.from_file(f'{receptor_prefix}.sdf')
                except:
                    print ('There must be a properly named fragment.sdf and fragment.pdb file within the frag_box folder or other fragment box recepticle folder')
                    sys.exit(1)
            else:
                receptor_structure = parmed.load_file(f'{receptor_prefix}.pdb')

        if ligand_prefix:
            if convert_mol2_to_sdf_ligand == True:
                print(ligand_prefix)
                ligand_structure, ligand = self._convert_mol2_to_sdf(ligand_prefix)
            else:
                try:
                   ligand_structure = parmed.load_file(f'{ligand_prefix}.pdb')
                except:
                    print('There must be a properly named fragment.sdf and fragment.pdb file within the frag_box folder or other fragment box recepticle folder')
                    sys.exit(1)
                try:
                    ligand = Molecule.from_file(f'{ligand_prefix}.sdf')
                except:
                    print('There must be a properly named fragment.sdf and fragment.pdb file within the frag_box folder or other fragment box recepticle folder')
                    sys.exit(1)
        if receptor_prefix and ligand_prefix:
            structure = ligand_structure + receptor_structure
        elif receptor_prefix:
            structure = receptor_structure
        elif ligand_prefix:
            structure = ligand_structure
        else:
            raise(f'prepare_system() requires receptor_pdb and/or ligand_sdf')

           
        parmed_forcefield_kwargs = {'removeCMMotion': self.removeCMMotion, 'ewaldErrorTolerance': 5e-04,
            'nonbondedMethod': app.PME, 'constraints': False, 'rigidWater': self.rigidWater, 'hydrogenMass': self.hydrogen_mass}
        openmm_forcefield_kwargs = {'removeCMMotion': self.removeCMMotion, 'ewaldErrorTolerance': 5e-04,
            'nonbondedMethod': app.PME, 'constraints': True, 'rigidWater': self.rigidWater, 'hydrogenMass': self.hydrogen_mass}

        if ligand_structure:
            if small_mol_rec==True:
                parmed_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                     periodic_forcefield_kwargs=parmed_forcefield_kwargs, molecules=[receptor, ligand],
                    small_molecule_forcefield=self.small_molecule_forcefield, cache=f'{self.project_directory}/LIG.json')
                openmm_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                    periodic_forcefield_kwargs=openmm_forcefield_kwargs, molecules=[receptor, ligand],
                    small_molecule_forcefield=self.small_molecule_forcefield, cache=f'{self.project_directory}/LIG.json')
            else:
                parmed_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                    periodic_forcefield_kwargs=parmed_forcefield_kwargs, molecules=[ligand],
                    small_molecule_forcefield=self.small_molecule_forcefield, cache=f'{self.project_directory}/LIG.json')
                openmm_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                     periodic_forcefield_kwargs=openmm_forcefield_kwargs, molecules=[ligand],
                    small_molecule_forcefield=self.small_molecule_forcefield, cache=f'{self.project_directory}/LIG.json')

        else:
            parmed_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                periodic_forcefield_kwargs=parmed_forcefield_kwargs)
            openmm_system_generator = SystemGenerator(forcefields=[self.protein_forcefield, self.solvation_forcefield],
                 periodic_forcefield_kwargs=openmm_forcefield_kwargs)

        modeller = app.Modeller(structure.topology, structure.positions)
        
        parmed_system = parmed_system_generator.create_system(modeller.topology)
        openmm_system = openmm_system_generator.create_system(modeller.topology)
        
        integrator = openmm.LangevinIntegrator(self.temperature, self.collision_rate, self.timestep)
        
        try:
            platform = openmm.Platform.getPlatformByName('CPU')
            context = openmm.Context(openmm_system, integrator, platform)
        except Exception as e:
            platform = openmm.Platform.getPlatformByName('OpenCL')
            context = openmm.Context(openmm_system, integrator, platform)
        
        context.setPositions(modeller.positions)
        
        state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
        parmed_system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())
        parmed_system = parmed.openmm.load_topology(modeller.topology,
            parmed_system, xyz=state.getPositions(asNumpy=True))
        if self.write_gromacs_files:
            print('Saving Gromacs files...')
            parmed_system.save(f'{self.project_directory}/{ligand_pre}.gro', overwrite=True)
            parmed_system.save(f'{self.project_directory}/{ligand_pre}.top', overwrite=True)
            
          
        if self.write_openmm_files:
            with open(f'{self.project_directory}/integrator.xml', 'w') as f:
                f.write(openmm.XmlSerializer.serialize(integrator))
            with open(f'{self.project_directory}/state.xml','w') as f:
                f.write(openmm.XmlSerializer.serialize(state))
            with open(f'{self.project_directory}/system.xml','w') as f:
                f.write(openmm.XmlSerializer.serialize(parmed_system))
 