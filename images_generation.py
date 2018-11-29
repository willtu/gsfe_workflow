# coding: utf-8
# Get the input files for Ag3PO4 with different layers(***)
from pymatgen import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import *
from pymatgen.io.vasp.sets import MITRelaxSet
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.site_transformations import \
    ReplaceSiteSpeciesTransformation as rsst
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import numpy as np
import pandas as pd
import os
import copy
import yaml
import logging
pi = np.pi
mpr = MPRester()
# rsst = ReplaceSiteSpeciesTransformation


class StructureGenerator():

    """
    A class to obtain struture from initial poscar and get the new poscar for 
    establishing slab model
    """

    def __init__(self, user_input, make_supercell=False, supercell_matrix=None, 
                 isdope=False, dopant=None):
        """
        Args:
            input: the POSCAR file of initial structure, or the material_id
        """
        self.make_supercell = make_supercell
        self.supercell_matrix = supercell_matrix
        self.isdope = isdope
        self.dopant = dopant
        # if user_input is a path where store the POSCAR file
        if(os.path.isfile(user_input)):
            pos = Poscar.from_file(user_input)
            struc = SpacegroupAnalyzer(pos.structure)
            self.init_struc = struc.get_conventional_standard_structure()
        # user_input is a mp_id
        else:
            prim_struc = mpr.get_structure_by_material_id(user_input)
            struc = SpacegroupAnalyzer(prim_struc)
            self.init_struc = struc.get_conventional_standard_structure()
        self.formula = str(self.init_struc.composition.reduced_formula)
        logging.info("self.formula: {}".format(self.formula))

    def get_interplanar_spacing(self, structure, surface):
        xc = XRDCalculator()
        xrd_data = xc.get_xrd_data(structure)
        logging.info("xrd_data: {}".format(xrd_data))
        sort_data = pd.DataFrame(xrd_data, columns=["two_theta", "intensity", 
                                 "miller_index_dict", "d_hkl"])
        matcher = (1, 1, 1)
        logging.info("sort_data: {}".format(sort_data))
        flag = 0
        for mid in sort_data.miller_index_dict:
            logging.info("mid: {}".format(mid))
            for key, value in mid.items():
                key = reduce_vector(key)
                if (surface == key):
                    matcher = mid
                    flag = 1
                    break
                else:
                    logging.info("*"*50)
                    logging.info("key: {}".format(key))
                    logging.info("value {}".format(value))
                    logging.info("*"*50)
            if flag == 1:
                break
        spacing = sort_data[sort_data.miller_index_dict == matcher].d_hkl.values[0]
        return spacing

    def get_endpoints_structure(self, path, layer=5, surface='1 1 1', 
                                burger_vector='0 0 0', crystal_lattice='FCC'):
        """
        Based on the parameters introduced, constructing a slab structure 
        oriented on the surface.
        Args:
            layer: the layer of slab.
            surface: the slip plane.
            burger_vector: the burger vector of slip system, is the multiple
            of a constant and a
            slip direction which is in miller index
        Returns:
            a list contains:
                Structure: includs the input parameters for function
                    get_wf_neb_from_endpoints(), parent_slab and slip_slab,
                    and delta_coords is the input parameter I add for
                    get_wf_neb_from_endpoints to reserve the slip
                    direction from parent_slab to slip_slab.
                POSCAR_path: the path of the new POSCAR file
        """
        # step1 get the initial structure from POSCAR file

        title = "files4layer%d_%s_%s" % (layer, surface, burger_vector)
        new_path = os.path.join(path, title)
        if not os.path.isdir(new_path):
            os.mkdir(new_path)
        struc = self.init_struc
        Poscar(struc).write_file(new_path + '/initial_POSCAR')
        if (crystal_lattice == 'HCP'):
            logging.info("This is a HCP Structure")
            h, k, l = surface.split(' ')
            surface = (int(h), int(k), int(l))
            hcp_surface_dict = {(0, 0, 1): (0, 0, 0, 1), 
                                (1, 0, 0): (1, 0, -1, 0),
                                (1, 0, 1): (1, 0, -1, 1), 
                                (2, -1, 2): (2, -1, -1, 2)}
            four_coords_surface = hcp_surface_dict[surface]
            interplanar_spacing = self.get_interplanar_spacing(
                struc, four_coords_surface)
        else:
            h, k, l = surface.split(' ')
            surface = (int(h), int(k), int(l))
            interplanar_spacing = self.get_interplanar_spacing(struc, surface)
        slab_size = np.ceil(interplanar_spacing * (layer-1))
        lattice = struc.lattice.matrix
        slabs = generate_all_slabs(struc, max_index=3, min_slab_size=slab_size, 
                                   min_vacuum_size=15.0, center_slab=True)
        parent_slab = [slab for slab in slabs if slab.miller_index == surface][0]
        if self.make_supercell is True:
            parent_slab.make_supercell(self.supercell_matrix)
        ps_lattice_matrix = parent_slab.lattice.matrix
        cuc_lattice_matrix = np.matrix(struc.lattice.matrix)
        ab_across = np.cross(ps_lattice_matrix[0], ps_lattice_matrix[1])
        area = np.linalg.norm(ab_across)
        logging.info("slip area: {}".format(area))
        transform_factor = ps_lattice_matrix * cuc_lattice_matrix.I
        dburger_vector = np.matrix(burger_vector) * transform_factor.I
        dburger_vector = np.ma.round(dburger_vector, decimals=6)
        order = zip(range(len(parent_slab.frac_coords)), 
                    parent_slab.frac_coords, parent_slab.species)
        c_order = sorted(order, key=lambda x: x[1][2])
        coord_seq = []
        for j in c_order:
            coord_seq.append(j[1])
        new_lattice = parent_slab.lattice
        new_species = parent_slab.species
        parent_slab = Structure(lattice=new_lattice, coords=coord_seq, 
                                species=new_species)
        atom_number = len(parent_slab.frac_coords)
        # To judge whether to dope
        if self.isdope is True:
            logging.info("We are substitute the %s by %s" % (
                self.formula, self.dopant))
            if isinstance(self.dopant, dict):
                ism = self.dopant
            else:
                dope_num = atom_number/2
                ism = {dope_num: self.dopant}
            Dope = rsst(ism)
            parent_slab = Dope.apply_transformation(parent_slab)
        Poscar(parent_slab).write_file(new_path + '/parent_slab_POSCAR')
        a = float(dburger_vector[0][0])
        b = float(dburger_vector[0][1])
        c = float(dburger_vector[0][2])
        coord_seq_1 = []
        d_order = copy.deepcopy(c_order)
        critical_atom = int(atom_number/2)
        logging.info("critical_atom: {}".format(critical_atom))
        for i in d_order[0:critical_atom]:
            i[1][0] -= a/2
            i[1][1] -= b/2
            i[1][2] -= c/2
            coord_seq_1.append(i[1])
        logging.info("d_order: {} \n".format(d_order))
        logging.info("c_order after: {} \n".format(c_order))
        logging.info("parent_slab_1: {} \n".format(parent_slab))
        for i in d_order[critical_atom:int(atom_number)]:
            i[1][0] += a/2
            i[1][1] += b/2
            i[1][2] += c/2
            coord_seq_1.append(i[1])
        new_lattice = parent_slab.lattice
        new_species = parent_slab.species
        frac_coords = coord_seq_1
        slip_slab = Structure(lattice=new_lattice, coords=frac_coords, 
                              species=new_species)
        Poscar(slip_slab).write_file(new_path + '/slip_slab_POSCAR')
        parent_set = MITRelaxSet(parent_slab)
        parent_struc = parent_set.poscar.structure
        slip_set = MITRelaxSet(slip_slab)
        slip_struc = slip_set.poscar.structure
        delta_coords = slip_struc.frac_coords - parent_struc.frac_coords
        return({"Structure": [parent_slab, slip_slab, delta_coords], 
                "POSACR_path": new_path})

    def get_endpoints_structure_from_yaml(self, layer, yaml_file, store_path):
        finder = SpacegroupAnalyzer(self.init_struc)
        stream = open(yaml_file, 'r')
        doc = yaml.load(stream)
        tag = 0
        endpoints_structures = {}
        for key, value in doc.items():
            crystal_lattice = key
            content = value
            eg_struc = mpr.get_structure_by_material_id(content['example'])
            finder_eg = SpacegroupAnalyzer(eg_struc)
            if(finder.get_space_group_symbol() == finder_eg.get_space_group_symbol()):
                tag = 1
                break
        if(tag == 0):
            logging.error("Sorry, we only support the structure of FCC, BCC \
                           and HCP right now.")
            return None
        elif(tag == 1):
            slip_system = content['slip system']
            for slip_system_key, slip_system_value in slip_system.items():
                surface = slip_system_key
                if isinstance(slip_system_value, list):
                    logging.info("One surface corresponding to more than one \
                                  Burgers vector.")
                    for direction in slip_system_value:
                        burger_vector = direction
                        logging.info("surface: {}, burger_vector: {}".format(
                            surface, burger_vector))
                        data_dict = self.get_endpoints_structure(store_path, layer,
                                surface, burger_vector, crystal_lattice)
                        endpoints_structures["%s_%s" %(
                            surface, burger_vector)] = data_dict["Structure"]
                        POSACR_path = data_dict["POSACR_path"]
                        logging.info("The new poscar file is: {}".format(POSACR_path))
                else:
                    burger_vector = slip_system_value
                    logging.info("surface: {}, burger_vector: {}".format(
                            surface, burger_vector))
                    data_dict = self.get_endpoints_structure(store_path, 
                        layer, surface, burger_vector, crystal_lattice)
                    endpoints_structures["%s_%s" %(
                        surface, burger_vector)] = data_dict["Structure"]
                    POSACR_path = data_dict["POSACR_path"]
                    logging.info("The new poscar file is: {}".format(POSACR_path))
        return(endpoints_structures)