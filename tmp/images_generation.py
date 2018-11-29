# Get the input files for Ag3PO4 with different layers(***)
from pymatgen import Structure, Lattice, Element, Molecule
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.adsorption import *
from pymatgen.core.surface import *
from pymatgen.io.vasp import Vasprun
from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.io.vasp.sets import MPRelaxSet, MITRelaxSet
from pymatgen.matproj.rest import MPRester
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.coord_utils import lattice_points_in_supercell
from pymatgen.transformations.site_transformations import ReplaceSiteSpeciesTransformation
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import numpy as np
import pandas as pd
from numpy import linalg as LA
import os
import copy
import yaml
pi = np.pi
mpr = MPRester()
rsst = ReplaceSiteSpeciesTransformation

class StructureGenerator():


    """
    A class obtains the struture from initial poscar and get the new poscar for establishing slab model
    """

    def __init__(self, user_input, make_supercell=False,supercell_matrix=None, isdope=False, dopant=None):
        """
        Args:
            input: the POSCAR file of initial structure, or the material_id
        """
        self.make_supercell = make_supercell
        self.supercell_matrix = supercell_matrix
        self.isdope = isdope
        self.dopant = dopant
        #user_input is a path where store the POSCAR file
        if(os.path.isfile(user_input)):
            pos = Poscar.from_file(user_input)
            struc = SpacegroupAnalyzer(pos.structure)
            self.init_struc = struc.get_conventional_standard_structure()
        #user_input is a mp_id
        else:
            prim_struc = mpr.get_structure_by_material_id(user_input)
            struc = SpacegroupAnalyzer(prim_struc)
            self.init_struc = struc.get_conventional_standard_structure()
        self.formula = str(self.init_struc.composition.reduced_formula)
        print("self.formula:", self.formula)

    def get_interplanar_spacing(self, structure, surface):
        xc = XRDCalculator()
        xrd_data = xc.get_xrd_data(structure)
        print(xrd_data)
        sort_data = pd.DataFrame(xrd_data, columns=["two_theta", "intensity", "miller_index_dict", "d_hkl"])
        matcher = (1, 1, 1)
        print(sort_data)
        flag = 0
        for mid in sort_data.miller_index_dict:
            print("mid:", mid)
            for key, value in mid.items():
                key = reduce_vector(key)
                if (surface == key):
                    print("Bravo!!!!")
                    matcher = mid
                    flag = 1
                    break;
                else:
                    print("*"*50)
                    print("key:",key)
                    print("value:",value)
                    print("*"*50)
            if flag == 1:
                break;
        spacing = sort_data[sort_data.miller_index_dict == matcher].d_hkl.values[0]
        print(spacing)
        return spacing


    def get_endpoints_structure(self, path, layer=5, surface = '1 1 1', burger_vector = '0 0 0', crystal_lattice = 'FCC'):
        """
        Based on the parameters introduced, constructing a new poscar file.
        Return new poscar file.

        Args:
            layer: the layer of slab.
            surface: the slip plane.
            burger_vector: the burger vector of slip system, is the multiple of a constant and a
            slip direction which is in miller index
        """
        ##step1 get the initial structure from POSCAR file

        title = "files4layer%d_%s_%s" %(layer,surface,burger_vector)
        new_path = os.path.join(path, title)
        if not os.path.isdir(new_path):
            os.mkdir(new_path)
        struc = self.init_struc
        print("initial struc:",struc)
        Poscar(struc).write_file(new_path + '/initial_POSCAR')
        if (crystal_lattice == 'HCP'):
            print("HCP Structure")
            h, k, l = surface.split(' ')
            surface = (int(h), int(k), int(l))
            hcp_surface_dict = {(0, 0, 1):(0, 0, 0, 1), (1, 0, 0):(1, 0, -1, 0), (1, 0, 1):(1, 0, -1, 1), (2, -1, 2):(2, -1, -1, 2)}
            four_coords_surface = hcp_surface_dict[surface]
            print("surface: ", four_coords_surface)
            interplanar_spacing = self.get_interplanar_spacing(struc, four_coords_surface)
        else:
            h, k, l = surface.split(' ')
            surface = (int(h), int(k), int(l))
            print("surface: ", surface)
            interplanar_spacing = self.get_interplanar_spacing(struc, surface)
        slab_size = np.ceil(interplanar_spacing * (layer-1))
        lattice = struc.lattice.matrix
        slabs = generate_all_slabs(struc, max_index=3, min_slab_size= slab_size, min_vacuum_size=15.0, center_slab=True)
        print(slabs)
        parent_slab = [slab for slab in slabs if slab.miller_index== surface][0]
        if self.make_supercell == True:
            parent_slab.make_supercell(self.supercell_matrix)
        ps_lattice_matrix = parent_slab.lattice.matrix
        cuc_lattice_matrix = np.matrix(struc.lattice.matrix)
        print("ps_lattice_matrix:", ps_lattice_matrix)
        print("ps_lattice_matrix.a:", ps_lattice_matrix[0])
        ab_across = np.cross(ps_lattice_matrix[0], ps_lattice_matrix[1])
        area = np.linalg.norm(ab_across)
        print("area:", area)
        print("cuc_lattice_matrix:", cuc_lattice_matrix)
        transform_factor = ps_lattice_matrix * cuc_lattice_matrix.I
        print("transform_factor:\n", transform_factor)
        dburger_vector = np.matrix(burger_vector) * transform_factor.I
        print("dburger_vector:\n", dburger_vector)
        dburger_vector = np.ma.round(dburger_vector, decimals=6)
        print("Transformed Burger_vector:\n",dburger_vector)
        order = zip(range(len(parent_slab.frac_coords)), parent_slab.frac_coords,  parent_slab.species)
        c_order = sorted(order, key = lambda x:x[1][2])
        coord_seq = []
        for j in c_order:
            coord_seq.append(j[1])
        print("test","*"*50)
        new_lattice = parent_slab.lattice
        new_species = parent_slab.species
        parent_slab = Structure(lattice=new_lattice, coords=coord_seq, species=new_species)
        atom_number = len(parent_slab.frac_coords)
        print("atom_number:",atom_number)
        # To judge whether to dope
        if self.isdope == True:
            print("We are substitute the %s by %s"%(self.formula, self.dopant))
            if isinstance(self.dopant, dict):
                ism = self.dopant
            else:
                dope_num = atom_number/2
                ism = {dope_num:self.dopant}
            Dope = rsst(ism)
            parent_slab = Dope.apply_transformation(parent_slab)
        Poscar(parent_slab).write_file(new_path + '/parent_slab_POSCAR')
        print("parent_slab: \n",parent_slab)
        a = float(dburger_vector[0][0])
        b = float(dburger_vector[0][1])
        c = float(dburger_vector[0][2])
        print(a,b,c)
        coord_seq_1 = []

        d_order = copy.deepcopy(c_order)
        print("id:",id(c_order),id(d_order))
        print("c_order before",c_order)
        #critical_atom = int(atom_number/2 - atom_number/layer)
        critical_atom = int(atom_number/2)
        print("critical_atom:",critical_atom)
        for i in d_order[0:critical_atom]:
            i[1][0] -= a/2
            i[1][1] -= b/2
            i[1][2] -= c/2
            coord_seq_1.append(i[1])
        print("d_order:\n", d_order)
        print("c_order after",c_order)
        print("parent_slab_1: \n",parent_slab)
        for i in d_order[critical_atom:int(atom_number)]:
            i[1][0] += a/2
            i[1][1] += b/2
            i[1][2] += c/2
            coord_seq_1.append(i[1])
        print("coord_seq_1:\n", coord_seq_1)
        new_lattice = parent_slab.lattice
        new_species = parent_slab.species
        frac_coords = coord_seq_1
        slip_slab = Structure(lattice=new_lattice, coords=frac_coords, species=new_species)
        print("c_order after again",c_order)
        print("parent_slab_2: \n",parent_slab)
        print("distance_matrix:",parent_slab.distance_matrix)
        print('slip slab:/n',slip_slab)
        Poscar(slip_slab).write_file(new_path + '/slip_slab_POSCAR')
        parent_set = MITRelaxSet(parent_slab)
        parent_struc = parent_set.poscar.structure
        slip_set = MITRelaxSet(slip_slab)
        slip_struc = slip_set.poscar.structure
        delta_coords = slip_struc.frac_coords-parent_struc.frac_coords
        return({"Structure":[parent_slab, slip_slab, delta_coords], "POSACR_path": new_path})

    def get_endpoints_structure_from_yaml(self, layer, yaml_file, store_path):
        finder = SpacegroupAnalyzer(self.init_struc)
        stream = open(yaml_file,'r')
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
        if(tag==0):
            print("Sorry, we only support the structure of FCC, BCC and HCP right now.")
        elif(tag==1):
            slip_system = content['slip system']
            for slip_system_key, slip_system_value in slip_system.items():
                surface = slip_system_key
                if isinstance(slip_system_value, list):
                    print("One surface to more than one burger vector!!!")
                    for direction in slip_system_value:
                        burger_vector = direction
                        print("surface:", surface, "burger_vector:", burger_vector)
                        data_dict = self.get_endpoints_structure(store_path, layer, surface, burger_vector, crystal_lattice)
                        print("data_dict:",data_dict)
                        endpoints_structures["%s_%s"%(surface, burger_vector)] = data_dict["Structure"]
                        POSACR_path = data_dict["POSACR_path"]
                        print("Hey, get the poscar file:", POSACR_path)
                else:
                    burger_vector = slip_system_value
                    print("surface:", surface, "burger_vector:", burger_vector)
                    data_dict = self.get_endpoints_structure(store_path, layer, surface, burger_vector, crystal_lattice)
                    print("data_dict:",data_dict)
                    endpoints_structures["%s_%s"%(surface, burger_vector)] = data_dict["Structure"]
                    POSACR_path = data_dict["POSACR_path"]
                    print("Hey, get the poscar file:", POSACR_path)
        return(endpoints_structures)