from .CONSTANT import AIRSS_INIT_POSCAR_SUBFIX
from .CONSTANT import M3GNET_RELAXED_POSCAR_SUBFIX


def gen_struc_tag(job_uuid, index=0):
    return f"{job_uuid}-{index}"


def gen_init_poscar_filename(struc_tag):
    return f"{struc_tag}{AIRSS_INIT_POSCAR_SUBFIX}"


def gen_m3g_poscar_filename(struc_tag):
    return f"{struc_tag}{M3GNET_RELAXED_POSCAR_SUBFIX}"


def get_struc_tag_from_init_poscar(poscar_filename):
    struc_tag = poscar_filename.replace(AIRSS_INIT_POSCAR_SUBFIX, '').split('/')[-1]
    return struc_tag


def get_struc_tag_from_m3g_poscar(poscar_filename):
    """Get the structure-tag from the poscar filename."""
    struc_tag = poscar_filename.replace(M3GNET_RELAXED_POSCAR_SUBFIX, '').split('/')[-1]
    return struc_tag


def gen_poscar_comment_line(tot_energy):
    comment = f"Energy(eV):{tot_energy}"
    return comment


def get_poscar_energy(poscar_comment):
    """Get the energy form the poscar comment."""
    tot_energy = float(poscar_comment.split(':')[1])
    return tot_energy



