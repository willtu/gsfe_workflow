from images_generation import StructureGenerator
import numpy as np
import pandas as pd
import json
import os
import math
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.matproj.rest import MPRester
from pymatgen import Lattice
from numpy import pi, dot, transpose, radians
from pymatgen.io.vasp import Poscar, Kpoints
from pymatgen import Structure, Lattice, Element
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.util.coord_utils import lattice_points_in_supercell
from fireworks import LaunchPad
from fireworks.core.launchpad import LAUNCHPAD_LOC
from atomate.vasp.workflows.base.neb import get_wf_neb_from_endpoints
import numpy as np
from numpy import linalg as LA
import os
import copy

# Design some parameters for building a slab,
# including the layer, element, slip system
layer = 12
mp_id = 'mp-1'
kpoint = 2100
yaml_file = "slip_system_1.yaml"
# instantiate the class StructureGenerator and get the initial bulk structure
#stru_gen = StructureGenerator(mp_id, make_supercell=True, supercell_matrix=[2,2,1])
#stru_gen = StructureGenerator(mp_id, make_supercell=True, supercell_matrix=[2,2,1], isdope=True, dopant="Ca")
stru_gen = StructureGenerator(mp_id)
formula = str(stru_gen.init_struc.composition.reduced_formula)
store_path = "/lustre/home/acct-umjzhh/umjzhh-1/tupeng/practice/Workflow/GSFE/example/Output/%s_%s/"%(mp_id,formula)
print(store_path)
if not os.path.isdir(store_path):
    os.mkdir(store_path)
endpoints_structures = stru_gen.get_endpoints_structure_from_yaml(layer, yaml_file, store_path)
print("-"*20+"First Part Ending"+"-"*20)
print(endpoints_structures)
print("-"*20+"Second Part Starting"+"-"*20)
lpad = LaunchPad.auto_load()
for endpoints_structures_key, endpoints_structures_value in endpoints_structures.items():
    print("-"*20+"Welcome to the magic world of Atomate"+"-"*20)
    print("endpoints_structures_key:\n",endpoints_structures_key)
    print('delta_coords:\n', endpoints_structures_value[-1])
    wf = get_wf_neb_from_endpoints(endpoints_structures_value[0],
    endpoints_structures_value[0:2], user_incar_settings=[
    {'ICHARG':2, 'EDIFF':5E-5, 'EDIFFG':-0.01, 'IBRION':2},
    {'ICHARG':2, 'EDIFF':5E-5, 'EDIFFG':-0.01, 'IBRION':2},
    {'IMAGES':5, 'EDIFF':1E-4, 'EDIFFG':-0.03, 'IBRION':3,'IOPT':2,'ENCUT':520,'ISYM':1, 'ICHAIN':0, 'SPRING':-5, 'LCLIMB':'TRUE', 'NSW':200, 'PREC':'Accurate', 'ALGO':'Normal', 'LREAL':'Auto', 'ISIF':2, 'ISMEAR':1, 'SIGMA':0.1}],
    user_kpoints_settings=[{"grid_density": kpoint},{"grid_density": kpoint},{"grid_density": kpoint}],
    additional_spec={'wf_name':'%s_%s_%d_%d_GSFE_NEB_Calcualtion_FCC_20180725'%(formula, endpoints_structures_key.replace(" ",""), kpoint,layer), 'delta_coords':endpoints_structures_value[-1]})
    lpad.add_wf(wf)
~