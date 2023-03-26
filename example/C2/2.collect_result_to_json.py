from aim3g.structure_analyser import StructureAnalyzer

RES_FOLDER = './result'

struc_analyzer = StructureAnalyzer(symprec=0.05)
struc_analyzer.analysis_all_m3g_poscar_in_folder(RES_FOLDER)
struc_analyzer.get_grouped_structures(ltol=0.1, stol=0.1, angle_tol=5)
struc_analyzer.publish_analysis_results()
