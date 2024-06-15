#%%
import pickle
import os

for n in [ 6, 7, 8, 9 ]:
    """
    Retrieve the irreps of z
    """
    dir_path = os.path.dirname(os.path.realpath("pickling.py")) + "/cgt/pickles_old/"
    dir_path_new = os.path.dirname(os.path.realpath("pickling.py")) + "/cgt/pickles/"
    full_path = dir_path + f"irreps_of_z_n{n}.pickle"
    full_path_new = dir_path_new + f"irreps_of_z_n{n}.pickle"
    print(full_path)
    try:
        with open(full_path, 'rb') as handle:
            irreps_of_z = pickle.load(handle)
    except FileNotFoundError:
        print(f"File not found: irreps_of_z_n{n}.pickle")
    
    irreps_of_z_new = []
    for irrep in irreps_of_z:
        irreps_of_z_new.append(irrep.transpose())

    try:
        with open(full_path_new, 'wb') as handle:
            pickle.dump(irreps_of_z_new, handle)
    except:
        print("A BADNESS HAPPENED")

# %%
