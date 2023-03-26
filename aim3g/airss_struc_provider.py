import os
import io
import warnings
import subprocess
from tqdm import tqdm

from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.io.vasp import Poscar

from ._misc_tools import get_uuid
from ._poscar_process import gen_struc_tag
from ._poscar_process import gen_init_poscar_filename

class AirssStrucProvider:
    def __init__(self, 
                 buildcell_exec, 
                 seedcell_file, 
                 cache_struc_in_memory=True) -> None:
        self.buildcell_exec = buildcell_exec
        self.seedcell_file = seedcell_file
        self.cache_struc_in_memory = cache_struc_in_memory
        self.cached_mymg_struc_list = []
        self.struc_tag = None
        self._ignore_warnings()


    def enable_cache_struc_in_memory(self, cache_struc_in_memory=True):
        self.cache_struc_in_memory = cache_struc_in_memory


    def clean_cached_struc(self):
        self.cached_mymg_struc_list = []


    def _ignore_warning(self, warning_head):
        warnings.filterwarnings("ignore", warning_head)
    def _ignore_warnings(self):
        self._ignore_warning("Generating CASTEP keywords JSON file... hang on.")
        self._ignore_warning("Could not determine the version of your CASTEP")
        self._ignore_warning("read_cell: Warning")


    def buildcell(self):
        """Run the buildcell process with the given seed.cell file"""
        with open(self.seedcell_file, "r") as seedinfo:
            buildcell_raw_out = subprocess.run(
                self.buildcell_exec, 
                stdin=seedinfo, 
                check=True, 
                text=True, 
                capture_output=True
            )
        return buildcell_raw_out
    
    def _castep_cell_to_pymatgen(self, castap_res_struc):
        """Tranform the castep cell file to ase structure data first, and then transfer the ase to pymatgen"""
        with io.StringIO(castap_res_struc) as res_file:
            ase_atoms = read(res_file, format="castep-cell")
            return AseAtomsAdaptor.get_structure(ase_atoms)


    def build_struc_pymg(self):
        buildcell_raw_out = self.buildcell()
        castap_res_struc = buildcell_raw_out.stdout
        mymg_struc = self._castep_cell_to_pymatgen(castap_res_struc)
        mymg_struc.info = {'tag': self.struc_tag}
        # If set to cache the struc, then cache it
        if self.cache_struc_in_memory:
            self.cached_mymg_struc_list.append(mymg_struc)
        return mymg_struc


    def build_struc_poscar(self):
        mymg_struc = self.build_struc_pymg()
        return Poscar(mymg_struc)


    def gen_batch_init_mymg_struc_in_memory(self, struc_quantity=1):
        self.enable_cache_struc_in_memory()
        job_uuid = get_uuid()
        print(f"[do] Generating init pymatgen structures into memory with UUID:{job_uuid}.")
        for i_struc in tqdm(range(struc_quantity)):
            self.struc_tag = gen_struc_tag(job_uuid, i_struc)
            self.build_struc_pymg()


    def gen_batch_init_poscar_files(self, output_path='.', struc_quantity=1):
        job_uuid = get_uuid()
        print(f"[do] Generating init POSCAR files with UUID:{job_uuid}.")
        for i_struc in tqdm(range(struc_quantity)):
            self.struc_tag = gen_struc_tag(job_uuid, i_struc)
            filename = gen_init_poscar_filename(self.struc_tag)
            filepath = os.path.join(output_path, filename)
            poscar = self.build_struc_poscar()
            with open(filepath, 'w') as fwp:
                fwp.write(str(poscar))




if __name__ == '__main__':
    buildcell_exec = "/home/cmtmd/PROJECT/M3GNet/airm3g/airss/bin/buildcell"
    seedcell_file = "/home/cmtmd/PROJECT/M3GNet/airm3g/example/LiFeSe/LiFeSe.cell"
    struc_provider = AirssStrucProvider(buildcell_exec, seedcell_file)

    struc_provider.gen_batch_init_poscar_files('.', 20)


