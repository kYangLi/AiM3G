from aim3g.structure_analyser import StructureAnalyzer

struc_analyzer = StructureAnalyzer()
struc_analyzer.read_groups_struc_infos_from_json_file()
struc_analyzer.print_grouped_structures()
