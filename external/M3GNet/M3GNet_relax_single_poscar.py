import os
import warnings

from m3gnet.models import Relaxer
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.structure_matcher import StructureMatcher

INIT_POSCAR_FILE = 'POSCAR.INIT'
DFT_RELAXED_POSCAR_FILE = 'POSCAR.DFT.RELAXED'
M3G_RELAXED_POSCAR_FILE = 'POSCAR.M3G.RELAXED'
RELAXED_STEPS = 10000
EXTRA_RELAX_TIMES = 5


# Parepare before the relax
for category in (UserWarning, DeprecationWarning):
    warnings.filterwarnings("ignore", category=category, module="tensorflow")
has_compare_dft_strcutre = os.path.isfile(DFT_RELAXED_POSCAR_FILE)


# Read in the POSCAR structures
init_structure = Poscar.from_file(INIT_POSCAR_FILE).structure


# Relax the structure and get the relaxed result
print("[do] Relaxing the structure...")
relaxer = Relaxer()  # This loads the default pre-trained model
relax_results = relaxer.relax(
    atoms=init_structure, steps=RELAXED_STEPS, verbose=True)
m3g_relaxed_structure = relax_results['final_structure']
for _ in range(EXTRA_RELAX_TIMES):
    relax_results = relaxer.relax(
        atoms=m3g_relaxed_structure, steps=RELAXED_STEPS, verbose=True)
    m3g_relaxed_structure = relax_results['final_structure']


# Save the relaxed POSCAR
print(f"[do] Printing the POSCAR file: '{M3G_RELAXED_POSCAR_FILE}' ...")
m3g_relaxed_poscar = Poscar(m3g_relaxed_structure)
with open(M3G_RELAXED_POSCAR_FILE, 'w') as fwp:
    fwp.write(str(m3g_relaxed_poscar))


# Print the result out
print(f"[do] Comparing the structures...")
matcher = StructureMatcher(ltol=1, stol=1, angle_tol=20,)
im_rms_dist = matcher.get_rms_dist(init_structure, m3g_relaxed_structure)
print(f"[info] init. vs. m3g : {im_rms_dist} aver/max")
if has_compare_dft_strcutre:
    dft_related_structure = Poscar.from_file(DFT_RELAXED_POSCAR_FILE).structure
    id_rms_dist = matcher.get_rms_dist(
        init_structure, dft_related_structure)
    dm_rms_dist = matcher.get_rms_dist(
        dft_related_structure, m3g_relaxed_structure)
    print(f"[info] init. vs. dft : {id_rms_dist} aver/max")
    print(f"[info] m3g   vs. dft : {dm_rms_dist} aver/max")







