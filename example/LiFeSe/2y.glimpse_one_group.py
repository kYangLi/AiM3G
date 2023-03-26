import argparse
from aim3g.structure_analyser import StructureAnalyzer

# Get the Group ip user want glimpse
parser = argparse.ArgumentParser(description='Process Observation Group ID.')
parser.add_argument('group_id', type=int, default=1,
                    help='Observation Group ID')
args = parser.parse_args()
group_id = args.group_id

# Glimpse the selected group
struc_analyzer = StructureAnalyzer()
struc_analyzer.read_groups_struc_infos_from_json_file()
struc_analyzer.glimpse_the_group(group_id)
