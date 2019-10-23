import numpy as np

# Reduce time array size, taking points before AND after a big step

def red_time1d(min_step, arr):

    nt = len(arr)
    ind_red = np.zeros(nt, dtype=np.bool)
    for jt in range(1, nt):
        if ( abs(arr[jt] - arr[jt-1]) > min_step ):
            ind_red[jt-1] = True
            ind_red[jt] = True

    ind_red[0]  = True
    ind_red[-1] = True

    return ind_red
 
    
def red_time(min_step, arr_2d):

# Several channels
    nt, nx = arr_2d.shape
    ind_red = [0]
    for jt in range(1, nt):
        for jx in range(nx):
            if ( abs(arr_2d[jt, jx] - arr_2d[jt-1, jx]) > min_step ):
                ind_red.append(jt-1)
                ind_red.append(jt)
                break

    ind_red.append(nt - 1)

    return np.unique(ind_red)
