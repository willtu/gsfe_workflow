# coding: utf-8
"""
This module defines a workflow for calculating generalized stacking fault
energy (GSFE)
"""

import os
import datetime
import logging
from images_generation import StructureGenerator
from fireworks import LaunchPad
from atomate.vasp.workflows.base.neb import get_wf_neb_from_endpoints

# Design some parameters for building a slab,
# including the layer, element, slip system
layer = 12
mp_id = 'mp-75'
kpoint = 2100
yaml_file = "slip_system.yaml"
# instantiate the class StructureGenerator and get the initial bulk structure
# stru_gen = StructureGenerator(mp_id, make_supercell=True, supercell_matrix=[2,2,1],
#                               isdope=True, dopant="Ca")
stru_gen = StructureGenerator(mp_id)
formula = str(stru_gen.init_struc.composition.reduced_formula)
store_path = "Output/%s_%s/" % (mp_id, formula)
if not os.path.isdir(store_path):
    os.mkdir(store_path)
endpoints_structures = stru_gen.get_endpoints_structure_from_yaml(
    layer, yaml_file, store_path)
logging.info("-"*20+"First Part Ending"+"-"*20)
logging.info("endpoints_structures: {}".format(endpoints_structures))
logging.info("-"*20+"Second Part Starting"+"-"*20)
lpad = LaunchPad.auto_load()
for key, value in endpoints_structures.items():
    logging.info("-"*20+"Welcome to the magic world of Atomate"+"-"*20)
    logging.info('delta_coords: {} \n'.format(value[-1]))
    wf = get_wf_neb_from_endpoints(value[0], value[0:2], user_incar_settings=[
        {'ICHARG': 2, 'EDIFF': 5E-5, 'EDIFFG': -0.01, 'IBRION': 2},
        {'ICHARG': 2, 'EDIFF': 5E-5, 'EDIFFG': -0.01, 'IBRION': 2},
        {'IMAGES': 5, 'EDIFF': 1E-4, 'EDIFFG': -0.03, 'IBRION': 3, 'IOPT': 2,
         'ENCUT': 520, 'ISYM': 1, 'ICHAIN': 0, 'SPRING': -5, 'LCLIMB': 'TRUE',
         'NSW': 200, 'PREC': 'Accurate', 'IALGO': 48, 'ALGO': 'Fast',
         'LREAL': 'Auto', 'ISIF': 2, 'ISMEAR': 1, 'SIGMA': 0.1}],
        user_kpoints_settings=[{"grid_density": kpoint},
        {"grid_density": kpoint}, {"grid_density": kpoint}],
        additional_spec={'wf_name': '%s_%s_%d_%d_GSFE_NEB_Calcualtion_%s' %
        (formula, key.replace(" ", ""), kpoint, layer,
         datetime.datetime.today().strftime("%Y%m%d")),
         'delta_coords': value[-1]})
    lpad.add_wf(wf)
