import time
from aim3g import M3GNetRelaxer

# Useful constant
STRUC_INPUT_FOLDER = './input'
STRUC_OUTPUT_FOLDER = './result'

# Time stamp record
start_time = time.time()

# Structure relaxer work
m3g_relaxer = M3GNetRelaxer(output_path=STRUC_OUTPUT_FOLDER)
m3g_relaxer.continue_relax(STRUC_INPUT_FOLDER)

# Time stamp record and output
end_time = time.time()
run_time = end_time - start_time
print(f"[info] Total time: {run_time:.2f} secondes.")
