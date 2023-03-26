import os
import time
from aim3g import AirssStrucProvider
from aim3g import M3GNetRelaxer


# Useful constant
STRUC_QUANTITY = 50000
SEEDCELL_FILE = "./C2.cell"
# - May not need to be changed on the same PC
BUILDCELL_EXEC = "../../external/AIRSS/bin/buildcell"
STRUC_INPUT_FOLDER = './input'
STRUC_OUTPUT_FOLDER = './result'
CACHE_STRUC_IN_MEMORY = True

# The file perperation
if not os.path.exists(STRUC_INPUT_FOLDER):
    os.mkdir(STRUC_INPUT_FOLDER)
if not os.path.exists(STRUC_OUTPUT_FOLDER):
    os.mkdir(STRUC_OUTPUT_FOLDER)

# Time stamp record
start_time = time.time()

# Structure provider wrok
struc_provider = AirssStrucProvider(
    BUILDCELL_EXEC, 
    SEEDCELL_FILE,
    cache_struc_in_memory=CACHE_STRUC_IN_MEMORY)
struc_provider.gen_batch_init_poscar_files(STRUC_INPUT_FOLDER, STRUC_QUANTITY)
all_strucs = struc_provider.cached_mymg_struc_list

# Time stamp record
struc_gen_time = time.time() - start_time

# Structure relaxer wrok
m3g_relaxer = M3GNetRelaxer(output_path=STRUC_OUTPUT_FOLDER)
m3g_relaxer.relax_all_pymg_strucs_in_list(all_strucs)

# Time stamp record and output
end_time = time.time()
run_time = end_time - start_time
print(f"[info] Generate structures time: {struc_gen_time:.2f} secondes.")
print(f"[info] Total time: {run_time:.2f} secondes.")

