from flow import FlowProject
from flow.environment import DefaultSlurmEnvironment
import sys
import mosdef_cassandra as mc
import generate_mc as gen   # Chris' script
import mbuild as mb
import foyer
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
import random
import unyt as u    # For Cassandra's physical unit specification 
                    # https://mosdef-cassandra.readthedocs.io/en/latest/guides/unyts.html


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"
    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition", default="batch", help="Specify the partition to submit to."
        )
    
def gen_system(dims=[2, 2, 2]):
    defaults_dict = {'stoichiometry': {'Mo': 1, 'V': 0.15, 'Nb': 0.13, 'Te': 0.12},
                     'dimensions': dims,
                     'template': '/home/erjank_project/nealeellyson/CassandraSim_signac/M1UnitCell_mod.pdb',
                     'crystal_x': 2.148490,
                     'crystal_y': 2.664721,
                     'crystal_z': 0.400321,
                     'z_reactor_size': 20.0,
                     'forcefield': None
                     }

    system = gen.Surface(defaults_dict['dimensions'],
                        defaults_dict['template'],
                         defaults_dict['stoichiometry'],
                         True,
                         defaults_dict['crystal_x'],
                         defaults_dict['crystal_y'],
                         defaults_dict['crystal_z']
                        )
    system.translate_to(np.zeros(3)) # Added so that periodicity doesn't error out
    system.periodicity = system.boundingbox.lengths # necessary
    #system.to_parmed() # Triggers an error
    return system

def apply_ff(system, systemFF, organicFF=None):
    typed_surface = systemFF.apply(system.children[0],
                            assert_angle_params=False,
                            assert_bond_params=False,
                            assert_dihedral_params=False,
                            assert_improper_params=False)
    if organicFF:
        for idx, child in enumerate(system.children):
            if idx == 0:
                pass
            else:
                typed_comp = organicFF.apply(child)
                typed_surface += typed_comp
    return typed_surface

def run_cassandra_surf(chem_pot, temp,steps_run, steps_restart):
    # Create sim box
    L = 10
    box = mb.Box(mins=[-L/2]*3, maxs=[L/2]*3)

    # Create materials
    surface = gen_system(dims=[1, 1, 2]) 
    surface.periodicity = surface.periodicity + np.array([0,0,1.4-0.85032099]) # Make (14 in A, convert from nm), only room for an ethane
    surface.translate_to([0,0,0]) # Centering
    ethane = mb.load("CC", smiles=True)

    # Load forcefields
    opls_uff = foyer.forcefields.Forcefield(forcefield_files = '/home/erjank_project/nealeellyson/CassandraSim_signac/forcefields/FF_opls_uff.xml')
    oplsaa = foyer.forcefields.load_OPLSAA()

    # Use foyer to apply forcefields
    typed_surface = apply_ff(surface, opls_uff)
    typed_ethane = oplsaa.apply(ethane, assert_bond_params=False, assert_angle_params=False, 
        assert_dihedral_params=False)
    
    # Create box and species list
    box_list = [surface]
    species_list = [typed_surface,typed_ethane]

    # Since we have an occupied box we need to specify
    # the number of each species present in the intial config
    mols_in_boxes = [[2,0]]

    system = mc.System(box_list, species_list, mols_in_boxes=mols_in_boxes)
    moveset = mc.MoveSet("gcmc", species_list)

    custom_args = {
        "run_name": f"surfequil_{chem_pot:.0f}_{temp:.0f}",
        "chemical_potentials": ["none",chem_pot*u.Unit('kJ/mol')],
        "rcut_min": 0.3980 * 2.5* u.angstrom, #(or 3.0)
        "vdw_cutoff": min(box.lengths)/2.1* u.angstrom,
        "charge_style": "none",
        #"charge_cutoff": 14.0,
        "coord_freq": 100,
        "prop_freq": 10,
    }

    mc.run(
        system=system, 
        moveset=moveset, 
        run_type="equilibration", 
        run_length= steps_run, #10000, # To reach ~1.33 hours
        temperature=temp*u.K, 
        **custom_args
    )
    
    # Set max translate and volume for production
    moveset.max_translate = [[0* u.angstrom,10.0* u.angstrom]] # angstroms

    # Update run_name and restart_name
    custom_args["run_name"] = f"surfprod_{chem_pot:.0f}_{temp:.0f}"
    custom_args["restart_name"] = f"surfequil_{chem_pot:.0f}_{temp:.0f}"

    mc.restart(
        system=system,
        moveset=moveset,
        run_type="production",
        run_length= steps_restart, #50000, # To reach ~6.67 hours, total 8 hrs
        temperature=temp*u.K,
        **custom_args,
    )
@FlowProject.operation
@FlowProject.post.isfile("box1.in.xyz")
def run_sim(job):
    with job:
        

        # Run simulation for spread of chemical potentials at set temperature
        run_cassandra_surf(job.sp.chem_pot,
                           job.sp.T, 
                           job.doc.steps_restart,
                           job.doc.steps_run)

if __name__ == '__main__':
    FlowProject().main()