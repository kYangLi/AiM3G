import os
import warnings
from pymatgen.io.vasp import Poscar

from ._misc_tools import get_uuid
from ._misc_tools import find
from ._misc_tools import get_filename_from_path
from ._poscar_process import gen_struc_tag
from ._poscar_process import gen_m3g_poscar_filename
from ._poscar_process import gen_poscar_comment_line
from ._poscar_process import gen_init_poscar_filename
from ._poscar_process import get_struc_tag_from_init_poscar
from ._poscar_process import get_struc_tag_from_m3g_poscar
from .CONSTANT import M3GNET_MAX_RELAX_STEPS
from .CONSTANT import AIRSS_INIT_POSCAR_SUBFIX
from .CONSTANT import M3GNET_RELAXED_POSCAR_SUBFIX



class M3GNetRelaxer:
    def __init__(self, output_path) -> None:
        from m3gnet.models import Relaxer
        self.relaxer = Relaxer()
        self.output_path = output_path # The relax result file output path
        self.relaxed_struc_num = 0
        self.tot_init_struc_num = -1
        self._ignore_warnings()


    def _ignore_warnings(self):
        for category in (UserWarning, DeprecationWarning):
            warnings.filterwarnings(
                "ignore", category=category, module="tensorflow")


    def relax_pymg_struc(self, pymg_struct):
        """Relax the structure in pymatgen format"""
        # This loads the default pre-trained model
        relax_results = self.relaxer.relax(
                            atoms=pymg_struct,
                            fmax=0.1,
                            interval=1,
                            steps=M3GNET_MAX_RELAX_STEPS, 
                            verbose=False)
        return relax_results


    def refresh_pymg_result_list(self):
        self.relaxed_pymg_strucs_list = []


    def _get_struc_tag(self, pymg_struc):
        if hasattr(pymg_struc, 'info') and 'tag' in pymg_struc.info:
            return pymg_struc.info['tag']
        else:
            return False


    def _sync_struc_tag(self, pymg_struc, relax_results):
        struc_tag = self._get_struc_tag(pymg_struc)
        if not struc_tag:
            uuid  = get_uuid()
            struc_tag = gen_struc_tag(uuid)
        relax_results['final_structure'].info = {'tag': struc_tag}


    def _print_pre_info(self, i_struc):
        print(f"[do] Relaxing Structure : [ {i_struc} of {self.tot_init_struc_num} ]")


    def _print_after_info(self, relax_results):
        final_energy  = relax_results['trajectory'].energies[-1]
        relax_steps = len(relax_results['trajectory'].energies)
        print(f" +--[result] Relax step: {relax_steps}/{M3GNET_MAX_RELAX_STEPS}")
        print(f" +--[result] Final energy: {final_energy} eV")


    def _get_pymg_struc_from_poscar(self, poscar_filename):
        pymg_struc = Poscar.from_file(poscar_filename).structure
        struc_tag = get_struc_tag_from_init_poscar(poscar_filename)
        pymg_struc.info = {'tag' : struc_tag}
        return pymg_struc

    
    def _poscar_in_exclude_files(self, exclude_files, poscar_path):
        poscar_filename = get_filename_from_path(poscar_path)
        return poscar_filename in exclude_files


    def _get_pymg_struc_list(self, poscar_file_list, exclude_files):
        pymg_struc_list = []
        for poscar_file in poscar_file_list:
            if self._poscar_in_exclude_files(exclude_files, poscar_file):
                print(f"[skip] Skipping file: {poscar_file}")
                continue
            pymg_struc = self._get_pymg_struc_from_poscar(poscar_file)
            pymg_struc_list.append(pymg_struc)
        return pymg_struc_list


    def relax_all_pymg_strucs_in_list(self, pymg_struc_list):
        self.refresh_pymg_result_list()
        self.tot_init_struc_num = len(pymg_struc_list) + self.relaxed_struc_num
        for i_struc, pymg_struc in enumerate(pymg_struc_list):
            self._print_pre_info(self.relaxed_struc_num + i_struc)
            relax_results = self.relax_pymg_struc(pymg_struc)
            self._sync_struc_tag(pymg_struc, relax_results)
            self.relaxed_pymg_strucs_list.append(relax_results)
            self.save_relaxed_res_to_poscar_file(relax_results)
            self. _print_after_info(relax_results)            


    def relax_all_list_poscar_files(self, poscar_file_list, exclude_files=[]):
        pymg_struc_list = \
            self._get_pymg_struc_list(poscar_file_list, exclude_files)
        self.relax_all_pymg_strucs_in_list(pymg_struc_list)


    def get_relaxed_init_poscars(self, output_folder):
        relaxed_init_poscars = [] # The list of filename of relaxed poscars
        relaxed_poscars_list = find(output_folder, M3GNET_RELAXED_POSCAR_SUBFIX)
        for relaxed_poscar in relaxed_poscars_list:
            relaxed_poscar_tag = get_struc_tag_from_m3g_poscar(relaxed_poscar)
            relaxed_init_poscar = gen_init_poscar_filename(relaxed_poscar_tag)
            relaxed_init_poscars.append(relaxed_init_poscar)
        return relaxed_init_poscars


    def get_poscar_file_list(self, input_folder):
        return find(input_folder, AIRSS_INIT_POSCAR_SUBFIX)


    def continue_relax(self, input_folder):
        init_poscars = self.get_poscar_file_list(input_folder)
        relaxed_init_poscars = self.get_relaxed_init_poscars(self.output_path)
        self.relaxed_struc_num = len(relaxed_init_poscars)
        self.relax_all_list_poscar_files(init_poscars, relaxed_init_poscars)


    def save_relaxed_res_to_poscar_file(self, relax_results):
        struc_tag = relax_results['final_structure'].info['tag']
        filename = gen_m3g_poscar_filename(struc_tag)
        with open(os.path.join(self.output_path, filename), 'w') as fwp:
            final_energy = relax_results['trajectory'].energies[-1]
            comment = gen_poscar_comment_line(final_energy)
            poscar = Poscar(relax_results['final_structure'], comment)
            fwp.write(str(poscar))
