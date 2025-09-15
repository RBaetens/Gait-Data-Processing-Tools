##################################################
# Getter functions to get stuff from marker data #
##################################################

import numpy as np
from scipy.signal import butter, sosfiltfilt
from lat_helpers import *

# Get angle of orientation of foot with respect to global x-axis
def get_foot_prog_angle(arrHEE, arrTOE, arrGRF):
    contact_mask = np.sum(arrGRF, axis=1)
    
    HEE = np.array([np.mean(mask_squeeze(arrHEE[:, 0], contact_mask)), np.mean(mask_squeeze(arrHEE[:, 1], contact_mask))])
    TOE = np.array([np.mean(mask_squeeze(arrTOE[:, 0], contact_mask)), np.mean(mask_squeeze(arrTOE[:, 1], contact_mask))])
    foot_dir = TOE - HEE

    # The following also works, despite no longer representing marker locations
    #HEE = np.array([np.mean(np.multiply(contact_mask, arrHEE[:, 0])), np.mean(np.multiply(contact_mask, arrHEE[:, 1]))])
    #TOE = np.array([np.mean(np.multiply(contact_mask, arrTOE[:, 0])), np.mean(np.multiply(contact_mask, arrTOE[:, 1]))])
    #foot_dir = TOE - HEE

    if np.sum(contact_mask) != 0:
        angle = np.arccos(foot_dir[0] / np.linalg.norm(foot_dir))
    else:
        angle = 0
        
    if foot_dir[1] < 0:
        angle = -angle
    
    return angle

# Get speed of a marker
def get_marker_speed(marker_arr, f_sample=100, filt_freq=None):
    if np.shape(marker_arr)[1] != 3:
        raise ValueError("Marker array must have second dimension of size 3")
    if np.shape(marker_arr)[0] <= 1:
        raise ValueError("Marker array must have first dimension of size larger than 1")
        
    diff_arr = np.diff(marker_arr, axis=0)
    speed_arr = np.sqrt(np.sum(np.power(diff_arr, 2), axis=1)) * (f_sample/1000)

    if filt_freq and not(np.sum(np.isnan(speed_arr))):
        filter = butter(2, filt_freq, btype='low', analog=False, output='sos', fs=f_sample)
        speed_arr = sosfiltfilt(filter, speed_arr)
    
    return speed_arr

# Get acceleration of a marker
def get_marker_accel(marker_arr, f_sample=100, filt_freq=None):
    if np.shape(marker_arr)[1] != 3:
        raise ValueError("Marker array must have second dimension of size 3")
    if np.shape(marker_arr)[0] <= 1:
        raise ValueError("Marker array must have first dimension of size larger than 1")
        
    ddiff_arr = np.diff(np.diff(marker_arr, axis=0), axis=0)
    accel_arr = np.sqrt(np.sum(np.power(ddiff_arr, 2), axis=1)) * ((f_sample**2)/1000)

    if filt_freq and not(np.sum(np.isnan(accel_arr))):
        filter = butter(2, filt_freq, btype='low', analog=False, output='sos', fs=f_sample)
        accel_arr = sosfiltfilt(filter, accel_arr)
    
    return accel_arr

# Get indices of "marker strikes"
def get_marker_strikes(arrMARK, arrGRF, f_sample=100, threshold_speed=0.2, filt_freq=None):
    # Find when foot first strikes plate, bit earlier
    start_ind = np.argmax(np.sum(arrGRF, axis=1))-int(0.1*f_sample)

    # Get marker speed
    speed_arr = get_marker_speed(arrMARK, f_sample=f_sample, filt_freq=filt_freq)

    # Get estimate of strikes
    marker_strikes = []
    n_samples = np.shape(speed_arr)[0]
    for i, speed in enumerate(speed_arr):
        
        # Conditions on index
        if (i >= start_ind) & (i > 5) & (i < n_samples-3):
            
            # Conditions on value and neighbouring values
            cond1 = (speed_arr[i] < threshold_speed)
            cond2 = (speed_arr[i-1] >= threshold_speed)
            cond3 = (np.sum(speed_arr[i+1:i+4] < threshold_speed) >= 2)
            cond4 = (np.sum(speed_arr[i-3:i] >= threshold_speed) >= 2)
            cond5 = (speed_arr[i-int(0.2*f_sample)] >= threshold_speed)
            if cond1 & cond2 & cond3 & cond4 & cond5:
                
                # Conditions on neighbouring strikes
                if len(marker_strikes) == 0:
                    marker_strikes.append(i)
                else:
                    if (i - marker_strikes[-1]) > 0.5*f_sample:
                        marker_strikes.append(i)

    return marker_strikes

# Get indices of "marker strikes", for training data calculations
def get_marker_strikes_train(arrMARK, f_sample=100, threshold_speed=0.2, filt_freq=None):
    # Get marker speed
    speed_arr = get_marker_speed(arrMARK, f_sample=f_sample, filt_freq=filt_freq)

    # Get estimate of strikes
    marker_strikes = []
    n_samples = np.shape(speed_arr)[0]
    
    for i, speed in enumerate(speed_arr):
        # Conditions on index
        if (i > 5) & (i < n_samples-3):
            # Conditions on value and neighbouring values
            cond1 = (speed_arr[i] < threshold_speed)
            cond2 = (speed_arr[i-1] >= threshold_speed)
            cond3 = (np.sum(speed_arr[i+1:i+4] < threshold_speed) >= 2)
            cond4 = (np.sum(speed_arr[i-3:i] >= threshold_speed) >= 2)
            cond5 = (speed_arr[i-int(0.2*f_sample)] >= threshold_speed)
            
            if cond1 & cond2 & cond3 & cond4 & cond5:
                # Conditions on neighbouring strikes
                if len(marker_strikes) == 0:
                    marker_strikes.append(i)
                else:
                    if (i - marker_strikes[-1]) > 0.5*f_sample:
                        marker_strikes.append(i)

    return marker_strikes

# Get "foot strikes" = heel or toe strikes, depending on what comes first
def get_foot_strikes(arrHEE, arrTOE, arrGRF, f_sample=100, threshold_speed=0.2, filt_freq=None):
    heel_strikes = get_marker_strikes(arrHEE, arrGRF,f_sample=f_sample, threshold_speed=threshold_speed, filt_freq=filt_freq)
    toe_strikes = get_marker_strikes(arrTOE, arrGRF,f_sample=f_sample, threshold_speed=threshold_speed, filt_freq=filt_freq)

    n_heel_strikes = len(heel_strikes)
    n_toe_strikes = len(toe_strikes)

    # If both strikes are found the same amount: select the ones coming first
    if n_heel_strikes == n_toe_strikes:
        foot_strikes = []
        for i in np.arange(n_heel_strikes):
            if heel_strikes[i] <= toe_strikes[i]:
                foot_strikes.append(heel_strikes[i])
            else:
                foot_strikes.append(toe_strikes[i])

    # If some toe strikes are missing: select heel strikes
    elif n_heel_strikes > n_toe_strikes:
        foot_strikes = heel_strikes

    # If some heel strikes are missing: select toe strikes
    else:
        foot_strikes = toe_strikes
        
    return foot_strikes

# Get "foot strikes" = heel or toe strikes, depending on what comes first, for training data calculations
def get_foot_strikes_train(arrHEE, arrTOE, f_sample=100, threshold_speed=0.2, filt_freq=None):
    heel_strikes = get_marker_strikes_train(arrHEE, f_sample=f_sample, threshold_speed=threshold_speed, filt_freq=filt_freq)
    toe_strikes = get_marker_strikes_train(arrTOE, f_sample=f_sample, threshold_speed=threshold_speed, filt_freq=filt_freq)

    n_heel_strikes = len(heel_strikes)
    n_toe_strikes = len(toe_strikes)

    # If both strikes are found the same amount: select the ones coming first
    if n_heel_strikes == n_toe_strikes:
        foot_strikes = []
        for i in np.arange(n_heel_strikes):
            if heel_strikes[i] <= toe_strikes[i]:
                foot_strikes.append(heel_strikes[i])
            else:
                foot_strikes.append(toe_strikes[i])

    # If some toe strikes are missing: select heel strikes
    elif n_heel_strikes > n_toe_strikes:
        foot_strikes = heel_strikes

    # If some heel strikes are missing: select toe strikes
    else:
        foot_strikes = toe_strikes
        
    return foot_strikes

# Get indices of "toe off"
def get_toe_offs(arrMARK, f_sample=100, threshold_speed=0.2, filt_freq=None):
    start_ind = 0

    # Get marker speed
    speed_arr = get_marker_speed(arrMARK, f_sample=f_sample, filt_freq=filt_freq)

    # Get estimate of toe offs
    toe_offs = []
    n_samples = np.shape(speed_arr)[0]
    for i, speed in enumerate(speed_arr):
        
        # Conditions on index
        if (i >= start_ind) & (i > 2) & (i < n_samples-int(0.2*f_sample)):
            
            # Conditions on value and neighbouring values
            cond1 = (speed_arr[i] > threshold_speed)
            cond2 = (speed_arr[i-1] <= threshold_speed)
            cond3 = (np.sum(speed_arr[i+1:i+4] > threshold_speed) >= 2)
            cond4 = (np.sum(speed_arr[i-3:i] <= threshold_speed) >= 2)
            cond5 = (speed_arr[i+int(0.2*f_sample)] >= threshold_speed)
            if cond1 & cond2 & cond3 & cond4 & cond5:
                
                # Conditions on neighbouring toe offs
                if len(toe_offs) == 0:
                    toe_offs.append(i)
                else:
                    if (i - toe_offs[-1]) > 0.5*f_sample:
                        toe_offs.append(i)

    return toe_offs