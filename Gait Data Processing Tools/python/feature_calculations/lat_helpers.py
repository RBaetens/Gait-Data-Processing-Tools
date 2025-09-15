##########################
# Very general functions #
##########################

import numpy as np
from scipy.signal import butter, sosfiltfilt

import copy

# Resample 2d grid along axis 0 using linear interpolation
def interp1d(data_arr, n_new):
    n_old, n_col = np.shape(data_arr)
    x_old = np.linspace(0, 1, num=n_old)
    x_new = np.linspace(0, 1, num=n_new)
    data_new = np.zeros((n_new, n_col))

    for i in np.arange(n_col):
        data_new[:, i] = np.interp(x_new, x_old, data_arr[:, i])

    return data_new

# Get elements from one dimensional array (column) as masked by a boolean array of the same size
def mask_squeeze(vector, bool_ind):
    masked_squeezed = np.zeros((np.sum(bool_ind), 1))
    j = 0
    for i, e in enumerate(bool_ind):
        if e:
            masked_squeezed[j] = vector[i]
            j += 1
    return masked_squeezed

# Take difference of elements in 2d array along axis 0, pad with order rows of nan at the start
def diff2d(arr, f_sample, order=1):
    if np.shape(arr)[0] <= order:
        raise ValueError("'arr' should have at least 'order'+1 rows")
    d_arr = np.diff(arr, n=order, axis=0)
    d_arr_pad = np.empty((order, np.shape(arr)[1]))
    d_arr_pad[:] = np.nan
    d_arr = np.concatenate((d_arr_pad, d_arr), axis=0)
    d_arr = np.multiply(d_arr, f_sample**order)
    return d_arr

# Shift 2d array along axis 0 for n_shift rows and replace n_shift rows with nan
def shiftnan2d(arr, n_shift):
    # Sanity checks
    n_pts = np.shape(arr)[0]
    if n_shift > n_pts:
        raise ValueError("The size of arr along axis 0 should be larger than or equal to n_shift.")

    # Actual shift
    res_arr = np.roll(arr, n_shift, axis=0)

    # Substitute part with nans
    if n_shift > 0:
        res_arr[:n_shift, :] = np.nan
    elif n_shift < 0:
        res_arr[n_shift:, :] = np.nan
    
    return res_arr

# Expand a 2d array of data to include shifted versions of the same data
def expand_data(arr_dat, n_shift, n_frames, crop=True):
    # Sanity check
    if n_frames <= 1:
        raise ValueError("n_frames should be greater than 1.")

    # Initialize
    arr_dat_exp = copy.copy(arr_dat)

    # Append shifted frames
    for i in np.arange(n_frames-1):
        arr_dat_exp = np.concatenate((arr_dat_exp, shiftnan2d(arr_dat, n_shift*(i+1))), axis=1)

    # Crop to not include rows with nans
    arr_dat_exp = arr_dat_exp[(n_frames-1)*n_shift:, :]
        
    return arr_dat_exp

# Impute NaN values and filter 2d grid along axis 0
def fill_filter(arr, f_sample, f_filter):
    n_row, n_col = np.shape(arr)
    nan_mask = np.isnan(arr)

    # Find blocks of NaN values
    i = 0
    while i < n_col:
        j = 0
        while j < n_row:
            if nan_mask[j, i]:
                start_block = j
                k = 0
                while nan_mask[j+k, i] and ((j+k) < n_row):
                    k += 1
                    if ((j+k) == n_row):
                        break
                end_block = j + k - 1
                n_block = end_block - start_block + 1
                
                # Impute (~linear interpolation)
                if n_block == 1:
                    arr[start_block, i] = np.nanmean(arr[:, i]) # value doesn't really matter as a singular point will be filtered out strongly
                elif ((start_block-n_block) >= 0) and ((end_block+n_block) < n_row):
                    if n_block <= 10:
                        arr[start_block:end_block+1, i] = np.ones(n_block) * np.nanmean(arr[start_block-n_block:end_block+1+n_block, i])
                    elif n_block <= 40:
                        mean1 = np.nanmean(arr[start_block-n_block:start_block, i])
                        rico1 = np.nanmean(np.diff(arr[start_block-n_block:start_block, i]))
                        mean2 = np.nanmean(arr[end_block+1:end_block+1+n_block, i])
                        rico2 = np.nanmean(np.diff(arr[end_block+1:end_block+1+n_block, i]))
                        arr[start_block:end_block+1, i] = np.ones(n_block) * (mean1 + n_block*rico1 + mean2 - n_block*rico2)/2
                    else:
                        mean1 = np.nanmean(arr[start_block-40:start_block, i])
                        rico1 = np.nanmean(np.diff(arr[start_block-40:start_block, i]))
                        mean2 = np.nanmean(arr[end_block+1:end_block+1+40, i])
                        rico2 = np.nanmean(np.diff(arr[end_block+1:end_block+1+40, i]))
                        arr[start_block:end_block+1, i] = np.ones(n_block) * (mean1 + (20+n_block/2)*rico1 + mean2 - (20+n_block/2)*rico2)/2
                elif ((start_block-n_block) >= 0):
                    mean1 = np.nanmean(arr[start_block-n_block:start_block, i])
                    rico1 = np.nanmean(np.diff(arr[start_block-n_block:start_block, i]))
                    arr[start_block:end_block+1, i] = np.ones(n_block) * (mean1 + n_block*rico1)
                elif ((end_block+n_block) < n_row):
                    mean2 = np.nanmean(arr[end_block+1:end_block+1+n_block, i])
                    rico2 = np.nanmean(np.diff(arr[end_block+1:end_block+1+n_block, i]))
                    arr[start_block:end_block+1, i] = np.ones(n_block) * (mean2 - n_block*rico2)
                else:
                    raise ValueError("Too many NaNs to handle meaningfully")
                j = end_block     
            j += 1
        i += 1
    
    filter = butter(2, f_filter, btype='low', analog=False, output='sos', fs=f_sample)
    arr = sosfiltfilt(filter, arr, axis=0)
    
    return arr

# Impute NaN values and filter 2d grid along axis 0
def fill_filter_robust(input_arr, f_sample, f_filter, n_imp=10):
    arr = copy.copy(input_arr)
    n_row, n_col = np.shape(arr)
    nan_mask = np.isnan(arr)

    # Find blocks of NaN values
    i = 0
    while i < n_col:
        j = 0
        while j < n_row:
            if nan_mask[j, i]:

                # Define block
                start_block = j
                k = 0
                while nan_mask[j+k, i] and ((j+k) < n_row):
                    k += 1
                    if ((j+k) == n_row):
                        break
                end_block = j + k - 1
                n_block = end_block - start_block + 1
                
                # Impute
                if ((start_block-n_imp) >= 0) and ((end_block+n_imp) < n_row):
                    mean1 = np.nanmean(arr[start_block-n_imp:start_block, i])
                    mean2 = np.nanmean(arr[end_block+1:end_block+n_imp+1, i])
                    arr[start_block:end_block+1, i] = np.ones(n_block) * (mean1 + mean2)/2
                elif (start_block-n_imp) >= 0:
                    mean1 = np.nanmean(arr[start_block-n_imp:start_block, i])
                    arr[start_block:end_block+1, i] = np.ones(n_block) * mean1
                elif (end_block+n_imp) < n_row:
                    mean2 = np.nanmean(arr[end_block+1:end_block+n_imp+1, i])
                    arr[start_block:end_block+1, i] = np.ones(n_block) * mean2
                else:
                    arr[start_block:end_block+1, i] = np.ones(n_block) * np.nanmean(arr[:, i])
                
                j = end_block     
            j += 1
        i += 1

    # Filter
    filter = butter(2, f_filter, btype='low', analog=False, output='sos', fs=f_sample)
    arr = sosfiltfilt(filter, arr, axis=0)
    
    return arr