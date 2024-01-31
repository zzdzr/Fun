import cooler
import numpy as np
import pandas as pd
import pylab
import math


def _transform(theta):

    import numpy as np
    return theta * np.pi / 180

def cos(theta):
    return np.cos(_transform(theta))

def sin(theta):
    return np.sin(_transform(theta))

def rotate(
    mat, theta, center, x_cord, y_cord,
    fill_mat = False, fill_value = 10
):
    '''
    :param theta: theta > 0, counter-clockwise, else clockwise
    '''

    # judge the orientation of rotation
    if theta >= 0:
        rotation_mat = np.array([[cos(theta), -sin(theta)],
                                 [sin(theta), cos(theta)]])
    else:
        theta = abs(theta)
        rotation_mat = np.array([[cos(theta), sin(theta)],
                                [-sin(theta), cos(theta)]])

    # Let's rotate the cord!
    cord_list = []
    for i, j in zip(x_cord, y_cord):
        x, y = i - center[0], j - center[1]
        tmp_cord = np.dot(rotation_mat, (x, y))
        tmp_cord[0] = tmp_cord[0] + center[0]
        tmp_cord[1] = tmp_cord[1] + center[1]
        cord_list.append(tmp_cord)

    # Extract the values
    value_list = []
    for i in cord_list:
        value_list.append(mat[int(i[0]), int(i[1])])
    value_list = np.asarray(value_list)
    value_list[np.isinf(value_list)] = np.nan

    # Fill the matrix
    if fill_mat:
        for i in cord_list:
            mat[int(i[0]), int(i[1])] = fill_value
        return mat

    return value_list



