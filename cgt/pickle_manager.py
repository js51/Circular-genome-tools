import pickle
import os

def retrieve_irrep_of_z(n=None):
    """
    Retrieve the irreps of z
    """
    dir_path = os.path.dirname(os.path.realpath(__file__)) + "/pickles/"
    full_path = dir_path + f"irreps_of_z_n{n}.pickle"
    try:
        with open(full_path, 'rb') as handle:
            irreps_of_z = pickle.load(handle)
    except FileNotFoundError:
        print(f"File not found: irreps_of_z_n{n}.pickle")
        return None
    
    return irreps_of_z