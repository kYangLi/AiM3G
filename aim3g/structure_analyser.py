import json
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from prettytable import PrettyTable
from tqdm import tqdm

from ._misc_tools import find
from ._misc_tools import min_with_index
from ._misc_tools import max_with_index
from ._misc_tools import calc_rms
from ._misc_tools import count
from ._misc_tools import sort_lists_by_list
from ._misc_tools import printdict
from ._poscar_process import get_struc_tag_from_m3g_poscar
from ._poscar_process import get_poscar_energy

from .CONSTANT import M3GNET_RELAXED_POSCAR_SUBFIX
from .CONSTANT import STRUC_GROUPS_INFO_JSON
from .CONSTANT import UNDETERMINED_SPACE_GROUP



class StructureAnalyzer:
    def __init__(self, symprec=0.01) -> None:
        self.symprec = symprec
        self.structures = []
        self.grouped_structures = []


    def _find_symmtery(self, pymg_struc, struc_tag):
        analyzer = SpacegroupAnalyzer(pymg_struc, symprec=self.symprec)
        try:
            space_group = analyzer.get_space_group_symbol()
        except BaseException:
            print("")
            print(f"[error] Fail to get the space group of structre '{struc_tag}', set to '{UNDETERMINED_SPACE_GROUP}'.")
            space_group = UNDETERMINED_SPACE_GROUP
        return space_group


    def _readin_str_and_energy_from_poscar(self, poscar_filename):
        poscar_obj = Poscar.from_file(poscar_filename)
        pymg_struc = poscar_obj.structure
        struc_tag = get_struc_tag_from_m3g_poscar(poscar_filename)
        energy = get_poscar_energy(poscar_obj.comment)
        symmetry = self._find_symmtery(pymg_struc, struc_tag)
        return struc_tag, energy, symmetry, pymg_struc
        

    def _push_poscar_data_to_memory(self, poscar_filename):
        """All of the data is pushed into the `self.structures_dict`"""
        struc_tag, energy, symmetry, pymg_struc = \
            self._readin_str_and_energy_from_poscar(poscar_filename)
        pymg_struc.info = {
            'tag': struc_tag, 
            'E': energy, 
            'symm': symmetry} # Be Careful Here!
        self.structures.append(pymg_struc)


    def _push_all_poscars_data_to_memory(self, poscar_list):
        """Push all of the poscar data in the list to memory dict. The `poscar_list` here indicate the poscar filename list."""
        for poscar_filename in tqdm(poscar_list, leave=True):
            self._push_poscar_data_to_memory(poscar_filename)


    def analysis_all_m3g_poscar_in_folder(self, folder='.'):
        print("[do] Reading in all relaxed structures into memory...")
        poscar_list = find(folder, M3GNET_RELAXED_POSCAR_SUBFIX)
        self._push_all_poscars_data_to_memory(poscar_list)


    def get_grouped_structures(self, ltol=0.2, stol=0.3, angle_tol=5):
        """
        Args:
            ltol (float): Fractional length tolerance. Default is 0.2.
            stol (float): Site tolerance. Defined as the fraction of the
                average free length per atom := ( V / Nsites ) ** (1/3)
                Default is 0.3.
            angle_tol (float): Angle tolerance in degrees. Default is 5 degrees.
        """
        print("[do] Grouping the structures...")
        matcher = StructureMatcher(ltol=ltol, stol=stol, angle_tol=angle_tol)
        self.grouped_structures = matcher.group_structures(self.structures)


    def _get_all_info_in_group(self, group, key):
        """Return all of the info to a list with the specific key."""
        key_val_list = []
        for ele in group:
            key_val_list.append(ele.info[key])
        return key_val_list


    def _get_group_infos(self, group):
        """Get the necessary group informations for printting to user."""
        # Get all of the extra info stored in the group
        tags = self._get_all_info_in_group(group, 'tag')
        energys = self._get_all_info_in_group(group, 'E')
        symmetrys = self._get_all_info_in_group(group, 'symm')
        # Sort the data using energys
        energys, tags, symmetrys = sort_lists_by_list(energys, tags, symmetrys)
        # Analysis those informations
        # - Reapeat times
        struc_repeat_time = len(group)
        # - Everage energy
        energy_average = sum(energys) / len(group)
        # - Energy RMS
        energy_rms = calc_rms(energys)
        # - Min energy structure infos
        min_energy, i_min = min_with_index(energys)
        min_energy_struc = {
            'tag' : tags[i_min],
            'energy' : min_energy,
            'symm' : symmetrys[i_min]
        }
        # - Max energy structure infos
        max_energy, i_max = max_with_index(energys)
        max_energy_struc = {
            'tag' : tags[i_max],
            'energy' : max_energy,
            'symm' : symmetrys[i_max]
        }
        # - Symmetry counts
        symmetrys_count = count(symmetrys)
        # Packup the group info
        group_info = {
            'repeat_time' : struc_repeat_time,
            'energy_rms' : energy_rms,
            'energy_aver' : energy_average,
            'min_energy_struc' : min_energy_struc,
            'max_energy_struc' : max_energy_struc,
            'symmetrys_count' : symmetrys_count,
            'metadata' : {
                'tags' : tags,
                'energys' : energys,
                'symmetrys' : symmetrys,
            },
        }
        return group_info


    def collect_groups_struc_infos(self):
        groups_info_list = []
        for group in self.grouped_structures:
            group_info = self._get_group_infos(group)            
            groups_info_list.append(group_info)
        # Sort it!
        self.groups_info_list = sorted(
            groups_info_list, 
            key=lambda x: x['min_energy_struc']['energy'])


    def print_grouped_structures(self):
        table = PrettyTable(['Group', 'Min Struc. Symm', 'Max Struc. Symm', 'Min Energy (eV)', 'Max Energy (eV)', 'Energy RMS (eV)', 'Repeat Time'])
        for i, group in enumerate(self.groups_info_list, start=1):
            min_energy_struc = group['min_energy_struc']
            max_energy_struc = group['max_energy_struc']
            energy_rms = group['energy_rms']
            repeat_time = group['repeat_time']
            min_energy = "{:.6f}".format(min_energy_struc['energy'])
            max_energy = "{:.6f}".format(max_energy_struc['energy'])
            row = [f"#{i}", f"{min_energy_struc['tag']}. {min_energy_struc['symm']}", f"{max_energy_struc['tag']}. {max_energy_struc['symm']}", min_energy, max_energy, "{:.6f}".format(energy_rms), repeat_time]
            table.add_row(row)
        print(table)

    
    def dump_groups_struc_infos_to_json_file(self):
        print(f"[do] Dumping the data dict to '{STRUC_GROUPS_INFO_JSON}' ...")
        with open(STRUC_GROUPS_INFO_JSON, 'w') as jfwp:
            json.dump(self.groups_info_list, jfwp, indent=4)


    def read_groups_struc_infos_from_json_file(self):
        with open(STRUC_GROUPS_INFO_JSON, 'r') as jfrp:
            self.groups_info_list = json.load(jfrp)


    def publish_analysis_results(self):
        print("[do] Printing the results...")
        self.collect_groups_struc_infos()
        self.print_grouped_structures()
        self.dump_groups_struc_infos_to_json_file()
        print("[done]")


    def glimpse_the_group(self, n_group):
        i_group = n_group - 1 # The NO. of the group and index of the group
        selected_group = self.groups_info_list[i_group]
        # Data reshow
        repeat_time = selected_group['repeat_time']
        energy_rms = "{:.6f}".format(selected_group['energy_rms'])
        energy_aver = "{:.6f}".format(selected_group['energy_aver'])
        min_energy_structure = selected_group['min_energy_struc']
        max_energy_structure = selected_group['max_energy_struc']
        symmetrys_count = selected_group['symmetrys_count']
        #:: The metadata has already be well sorted via energy in the self._get_group_infos() method. ::
        tags = selected_group['metadata']['tags']
        energys = selected_group['metadata']['energys']
        symmetrys = selected_group['metadata']['symmetrys']
        # Print the result
        log_filename = f'group_{n_group}.list'
        print(f"[info] Group Order         : {n_group}     ")
        print(f"[info] Structure Number    : {repeat_time} ")
        print(f"[info] Energy RMS (eV)     : {energy_rms}  ")
        print(f"[info] Average Energy (eV) : {energy_aver} ")
        print(f"[info] Minimal Energy Strucutre info.s :")
        printdict(min_energy_structure)
        print(f"[info] Maximal Energy Strucutre info.s :")
        printdict(max_energy_structure)
        print(f"[info] Symmetry Counts:")
        printdict(symmetrys_count)
        print(f'[do] Saving detailed info to {log_filename} ...')
        table = PrettyTable(['Structure', 'Energy (eV)', 'Space Group'])
        for tag, energy, spg in zip(tags, energys, symmetrys):
            table.add_row([tag, "{:.6f}".format(energy), spg])
        with open(log_filename, 'w') as fwp:
            fwp.write(str(table))
        print("[done]")


