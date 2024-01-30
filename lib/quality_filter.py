import numpy as np

def filter_extension(region_file, threshold=0.6):
    ext_list = region_file['perc_res_list']

    func = lambda x: len(x[x > threshold]) >= 1

    bool_list = [func(np.asarray(i)) for i in ext_list]

    return bool_list