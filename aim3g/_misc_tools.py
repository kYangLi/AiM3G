import os
import uuid
import numpy as np
import json

from .CONSTANT import JOB_UUID_LENGTH



def get_uuid(length=JOB_UUID_LENGTH):
    return str(uuid.uuid4()).upper().replace("-", "")[:length]


def find(folder, string):
    # sort filenames based on numeric value
    def sort_by_numeric_value(filename):
        numeric_part = ''.join(filter(str.isdigit, filename))
        return int(numeric_part)
    
    filepath_list = []
    files = os.listdir(folder)
    filtered_files = [f for f in files if string in f]
    filtered_files.sort(key=sort_by_numeric_value)
    for filename in filtered_files:
        filepath_list.append(os.path.join(folder, filename))
    return filepath_list


def get_filename_from_path(path):
    return os.path.basename(path)


def max_with_index(alist):
    """Return the max value and its index in the list"""
    max_value = max(alist)
    its_index = alist.index(max_value)
    return max_value, its_index


def min_with_index(alist):
    """Return the max value and its index in the list"""
    min_value = min(alist)
    its_index = alist.index(min_value)
    return min_value, its_index


def calc_rms(alist):
    return np.sqrt(np.var(alist))


def count(alist):
    """Count the number of each element inside the list, return there count by a dict."""
    # return dict(zip(*np.unique(alist, return_counts=True))) # Not JSON-able
    keys, counts = np.unique(alist, return_counts=True)
    counts = list(map(int, counts))
    return dict(zip(keys, counts))


def sort_lists_by_list(ref_list, *args):
    sorted_idx = np.argsort(ref_list)
    sorted_lists = []
    for alist in args:
        sorted_lists.append(list(np.array(alist)[sorted_idx]))
    sorted_ref_list = list(np.array(ref_list)[sorted_idx])
    sorted_lists.insert(0, sorted_ref_list)
    return tuple(sorted_lists)


def printdict(adict, tab_space=6):
    formatted_dict = json.dumps(adict, indent=4, sort_keys=True)
    for line in formatted_dict.splitlines():
        print(f"{' '*tab_space}{line}")
