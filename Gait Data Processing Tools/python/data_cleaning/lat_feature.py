##############################################################
# Functions to help calculate features from measurement data #
##############################################################

import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
from scipy.spatial.transform import Rotation

import os

from lat_globals import *
from lat_helpers import *


#################
### EMG stuff ###
#################

# Root mean square
def rms_feat(_, signal):
    rms = np.mean(np.power(signal, 2))
    return np.sqrt(rms)

# Returns the powerspectrum and accompanying frequency axis in a form that is nice for feature extraction
def freq_feat_spec(t_axis, signal):
    Nsamp = np.shape(t_axis)[0]
    Tmeasure = t_axis[-1]-t_axis[0]
    dt = Tmeasure/(Nsamp-1)

    spectrum = np.fft.fft(signal)
    spectrum = np.power(np.abs(spectrum), 2)
    f_axis = np.fft.fftfreq(Nsamp, d=dt)
    
    n_freq = np.shape(f_axis)[0]
    n_freq = int(np.ceil(n_freq/2))
    f_axis = f_axis[0:n_freq]
    spectrum = spectrum[0:n_freq]
    return f_axis, spectrum

# Mean frequency
def mean_freq_feat(t_axis, signal):
    f_axis, spectrum = freq_feat_spec(t_axis, signal)
    mean_freq_out = np.sum(np.multiply(f_axis, spectrum))/np.sum(spectrum)
    return mean_freq_out

# Get norms for EMG signals based on the 'Normal level walking' trials
def get_emg_norms(foldersToRead, f_emg=1000):
    # Define number of muscles and filter
    n_muscles = len(EMG_HEADER)
    filter = butter(2, 6, btype='low', analog=False, output='sos', fs=f_emg)
    
    # List the files with raw EMG data for the 'Normal level walking' trials
    emg_files = []
    for folder in foldersToRead:
        files = os.listdir(folder)
        for file in files:
            if (("EMG." in file) or ("EMG - " in file)) and ("Normal level walking" in file) and ("Rec1" in file):
                emg_files.append(folder + "/" + file)

    # To store results
    n_files = len(emg_files)
    emg_rms = np.zeros((n_files, n_muscles))
    emg_mnf = np.zeros((n_files, n_muscles))
    emg_act = np.zeros((n_files, n_muscles))
    
    # Loop over all these files and process
    for i, emg_file in enumerate(emg_files):
        # Read raw EMG data
        emg_arr = pd.read_csv(emg_file).to_numpy()
        
        # Calculate RMS
        emg_rms[i, :] = np.sqrt(np.mean(np.power(emg_arr, 2), axis=0))

        # Calculate mean frequency
        n_pts = np.shape(emg_arr)[0]
        t_axis = np.linspace(0, (n_pts-1)/f_emg, n_pts)
        for j in np.arange(n_muscles):
            emg_mnf[i, j] = mean_freq_feat(t_axis, emg_arr[:, j])

        # Calculate activation
        emg_arr = np.abs(emg_arr)
        emg_arr = sosfiltfilt(filter, emg_arr, axis=0)
        emg_act[i, :] = np.mean(emg_arr, axis=0)

    # Mean
    emg_rms = pd.DataFrame(np.mean(emg_rms, axis=0).reshape((1, -1)))
    emg_mnf = pd.DataFrame(np.mean(emg_mnf, axis=0).reshape((1, -1)))
    emg_act = pd.DataFrame(np.mean(emg_act, axis=0).reshape((1, -1)))

    # Save
    sync_folder_path = "/".join(emg_file.split("/")[:-1])
    if "Proc_Ref" in sync_folder_path:
        sync_folder_path = sync_folder_path.replace("Proc_Ref", "Sync")
    elif "Proc_Trn" in sync_folder_path:
        raise ValueError("You should use the Proc_Ref data to calculate the norms.")
    else:
        raise ValueError("Cannot find Sync folder.")
    emg_rms.to_csv(sync_folder_path + "/RMS_norms_walking.csv", index=False, float_format="%.6g", header=False)
    emg_mnf.to_csv(sync_folder_path + "/MNF_norms_walking.csv", index=False, float_format="%.6g", header=False)
    emg_act.to_csv(sync_folder_path + "/ACT_norms_walking.csv", index=False, float_format="%.6g", header=False)
    return

# Read the raw EMG data, calculate desired features and store in new files
def create_emg_feat_files(foldersToRead, f_emg=1000, f_markers=100, window=0.3, norm="cycle"):
    # Define what features you calculate
    n_muscles = len(EMG_HEADER)
    td_feat_func_dict = {"RMS": rms_feat}
    fd_feat_func_dict = {"MNF": mean_freq_feat}
    n_td_feat_func = len(td_feat_func_dict)
    n_td_feat = n_td_feat_func * n_muscles
    n_fd_feat_func = len(fd_feat_func_dict)
    n_fd_feat = n_fd_feat_func * n_muscles

    n_pts_window_emg = int(np.round(window*f_emg))
    n_pts_window_feat = int(np.round(window*f_markers))

    t_axis_window_emg = np.linspace(0, window, n_pts_window_emg)
    
    # List all the files with raw EMG data
    emg_files = []
    for folder in foldersToRead:
        if "Proc_Ref" in folder:
            sync_folder = folder.replace("Proc_Ref", "Sync")
        elif "Proc_Trn" in folder:
            sync_folder = folder.replace("Proc_Trn", "Sync")
        else:
            raise ValueError("Cannot find Sync folder.")
        
        if norm == "cycle":
            RMS_file = sync_folder + "/RMS_norms_walking.csv"
            RMS_vals = np.loadtxt(RMS_file, delimiter=",")
            MNF_file = sync_folder + "/MNF_norms_walking.csv"
            MNF_vals = np.loadtxt(MNF_file, delimiter=",")
            ACT_file = sync_folder + "/ACT_norms_walking.csv"
            ACT_vals = np.loadtxt(ACT_file, delimiter=",")
        else:
            RMS_file = sync_folder + "/MVC_norms.csv"
            RMS_vals = np.loadtxt(RMS_file, delimiter=",")
            MNF_file = sync_folder + "/MNF_norms_walking.csv" # Should still be changed!!!
            MNF_vals = np.loadtxt(MNF_file, delimiter=",")
            ACT_file = sync_folder + "/ACT_norms.csv"
            ACT_vals = np.loadtxt(ACT_file, delimiter=",")
        files = os.listdir(folder)
        for file in files:
            if ("EMG." in file) or ("EMG - " in file):
                emg_files.append((folder + "/" + file, RMS_vals, MNF_vals, ACT_vals))

    # Define filter for activation calculation
    filter = butter(2, 6, btype='low', analog=False, output='sos', fs=f_emg)

    # Loop over all these files and process
    for emg_file in emg_files:
        # Read raw EMG data
        emg_df = pd.read_csv(emg_file[0])
        RMS_vals = emg_file[1]
        MNF_vals = emg_file[2]
        ACT_vals = emg_file[3]
        
        n_pts_emg = emg_df.shape[0]
        n_pts_feat = int(np.round(n_pts_emg*f_markers/f_emg))

        # Calculate time domain features and store in DataFrame
        feat_arr = np.empty((n_pts_feat, n_td_feat))
        feat_arr[:] = np.nan
        emg_feat_header = []
        for i, feature_name in enumerate(td_feat_func_dict.keys()):
            emg_feat_header = emg_feat_header + [" ".join((feature_name, muscle_name)) for muscle_name in EMG_HEADER_SHORT]
            for j, muscle_name in enumerate(EMG_HEADER):
                emg_arr = emg_df[muscle_name].to_numpy()
                for k in np.arange(n_pts_window_feat-1, n_pts_feat):
                    k_emg = int(np.round((k+1)*f_emg/f_markers))
                    feat_arr[k, i*n_muscles+j] = td_feat_func_dict[feature_name](t_axis_window_emg, emg_arr[k_emg-n_pts_window_emg:k_emg])
                feat_arr[:, i*n_muscles+j] = feat_arr[:, i*n_muscles+j] / RMS_vals[j]
        emg_feat_df = pd.DataFrame(feat_arr, columns=emg_feat_header)

        # Calculate frequency domain features and append to DataFrame
        feat_arr = np.empty((n_pts_feat, n_fd_feat))
        feat_arr[:] = np.nan
        emg_feat_header = []
        for i, feature_name in enumerate(fd_feat_func_dict.keys()):
            emg_feat_header = emg_feat_header + [" ".join((feature_name, muscle_name)) for muscle_name in EMG_HEADER_SHORT]
            for j, muscle_name in enumerate(EMG_HEADER):
                emg_arr = emg_df[muscle_name].to_numpy()
                for k in np.arange(n_pts_window_feat-1, n_pts_feat):
                    k_emg = int(np.round((k+1)*f_emg/f_markers))
                    feat_arr[k, i*n_muscles+j] = fd_feat_func_dict[feature_name](t_axis_window_emg, emg_arr[k_emg-n_pts_window_emg:k_emg])
                feat_arr[:, i*n_muscles+j] = feat_arr[:, i*n_muscles+j] / MNF_vals[j]
        emg_feat_df[emg_feat_header] = feat_arr

        # Calculate activations and append to DataFrame
        act_arr = np.empty((n_pts_feat, 8))
        act_arr[:] = np.nan
        for i, muscle_name in enumerate(EMG_HEADER):
            emg_arr = emg_df[muscle_name].to_numpy()
            emg_arr = np.abs(emg_arr)
            emg_arr = sosfiltfilt(filter, emg_arr, axis=0)
            x_old = np.linspace(0, 1, n_pts_emg)
            x_new = np.linspace(0, 1, n_pts_feat)
            emg_arr = np.interp(x_new, x_old, emg_arr)
            act_arr[:, i] = emg_arr / ACT_vals[i]
        emg_feat_df[EMG_HEADER_ACT] = act_arr

        # Save in a new file
        emg_feat_file = emg_file[0].replace("EMG", "EMG_features")
        emg_feat_df.to_csv(emg_feat_file, index=False, float_format="%.6g")
    return


#################
### IMU stuff ###
#################

def get_simulated_IMU_measurements(segment_orientation_position, f_sample=100, seq="XYZ"):    
    # Restructure
    rot_vecs = segment_orientation_position[:, 0:3] # in radians
    translations = segment_orientation_position[:, 3:6] # in meters

    # Change orientation representation, rotation vector to rotation matrices
    n_pts = np.shape(rot_vecs)[0]
    rot_matrices = np.zeros((n_pts, 3, 3))
    inv_rot_matrices = np.zeros((n_pts, 3, 3))

    for i in np.arange(n_pts):
        rotation = Rotation.from_rotvec(rot_vecs[i, :])
        rot_matrices[i, :, :] = rotation.as_matrix()
        inv_rot_matrices[i, :, :] = np.linalg.inv(rot_matrices[i, :, :])
        #inv_rot_matrices[i, :, :] = rotation.inv().as_matrix()

    # Calculate rotational velocity and linear acceleration a measured from the segment frame of reference
    delta_rot_matrices = np.zeros((n_pts-1, 3, 3))
    delta_euler_angles = np.zeros((n_pts-1, 3))

    accelerations = np.diff(translations, n=2, axis=0)
    accelerations *= f_sample**2
    
    gravity = np.concatenate((np.zeros((n_pts-2, 2)), np.ones((n_pts-2, 1))*9.81), axis=1)
    accelerations += gravity

    for i in np.arange(n_pts-2):
        delta_rot_matrices[i, :, :] = np.matmul(rot_matrices[i+1, :, :], inv_rot_matrices[i, :, :])
        try:
            rotation = Rotation.from_matrix(delta_rot_matrices[i, :, :])
            delta_euler_angles[i, :] = np.matmul(inv_rot_matrices[i, :, :], rotation.as_euler(seq).reshape((-1, 1))).reshape((1, -1))
            accelerations[i, :] = np.matmul(inv_rot_matrices[i, :, :], accelerations[i, :].reshape((-1, 1))).reshape((1, -1))
        except (np.linalg.LinAlgError):
            delta_euler_angles[i, :] = np.nan
            accelerations[i, :] = np.nan

    i = n_pts-2
    delta_rot_matrices[i, :, :] = np.matmul(rot_matrices[i+1, :, :], inv_rot_matrices[i, :, :])
    try:
        rotation = Rotation.from_matrix(delta_rot_matrices[i, :, :])
        delta_euler_angles[i, :] = np.matmul(inv_rot_matrices[i, :, :], rotation.as_euler(seq).reshape((-1, 1))).reshape((1, -1))
    except (np.linalg.LinAlgError):
        delta_euler_angles[i, :] = np.nan

    delta_euler_angles *= f_sample
    
    return delta_euler_angles, accelerations

def create_imu_feat_files(foldersToRead, f_markers=100, seq="XYZ"):
    deg_to_rad = np.pi/180
    mm_to_m = 1/1000
    
    # List all the files with the global angles
    glo_ang_files = []
    for folder in foldersToRead:
        files = os.listdir(folder)
        for file in files:
            if ("GloAng" in file) or ("GlobalAngles" in file):
                glo_ang_files.append(folder + "/" + file)
    
    # Loop over all these files and process
    for glo_ang_file in glo_ang_files:
        # Read global angle data
        glo_ang_df = pd.read_csv(glo_ang_file)
        glo_ang_df = glo_ang_df[COLS_SEGMENTS_FLAT]
        n_pts = glo_ang_df.shape[0]
        
        # Initialize IMU feature DataFrame
        imu_feat_df = pd.DataFrame(np.zeros((n_pts, N_IMU_FEAT)), columns=COLS_IMU_FEAT_FLAT)
    
        # Calculate features and store in DataFrame, per body segment
        for i, COLS in enumerate(COLS_SEGMENTS):
            # Read one segment
            segment_orientation_position = glo_ang_df[COLS].to_numpy()
            segment_orientation_position[:, 0:3] = segment_orientation_position[:, 0:3] * deg_to_rad
            segment_orientation_position[:, 3:6] = segment_orientation_position[:, 3:6] * mm_to_m

            # Calculate measurements and filtered measurements
            delta_euler_angles, accelerations = get_simulated_IMU_measurements(segment_orientation_position, f_sample=f_markers, seq=seq)

            # Append rows of nan to match n_pts
            nan_row = np.zeros((1, 3))
            nan_row[:] = np.nan
            nan_rows = np.zeros((2, 3))
            nan_rows[:] = np.nan
            delta_euler_angles = np.concatenate((delta_euler_angles, nan_row), axis=0)
            accelerations = np.concatenate((accelerations, nan_rows), axis=0)

            # Update output dataFrame
            COLS_FEAT = COLS_IMU_FEAT[i]
            imu_feat_df[COLS_FEAT] = np.concatenate((delta_euler_angles, accelerations), axis=1)
            
        # Save in a new file
        if "GloAng" in glo_ang_file:
            imu_feat_file = glo_ang_file.replace("GloAng", "IMU_features")
        elif "GlobalAngles" in glo_ang_file:
            imu_feat_file = glo_ang_file.replace("GlobalAngles", "IMU_features")
        else:
            raise ValueError("Something went wrong with the name of the global angle file.")
        imu_feat_df.to_csv(imu_feat_file, index=False, float_format="%.6g")
    return


################################
### Encoder/goniometer stuff ###
################################

def create_gon_feat_files(foldersToRead, f_markers=100):    
    # List all the files with the model outputs
    mod_out_files = []
    for folder in foldersToRead:
        files = os.listdir(folder)
        for file in files:
            if ("ModOut" in file) or ("ModelOutputs" in file):
                mod_out_files.append(folder + "/" + file)
    
    # Loop over all these files and process
    for mod_out_file in mod_out_files:
        # Read model output data
        mod_out_df = pd.read_csv(mod_out_file)
        joint_angle_arr = mod_out_df[COLS_JOINT_ANGLES].to_numpy()

        # Calculate joint velocities
        joint_vel_arr = np.diff(joint_angle_arr, axis=0) * f_markers
        nan_row = np.zeros((1, np.shape(joint_angle_arr)[1]))
        nan_row[:] = np.nan
        joint_vel_arr = np.concatenate((joint_vel_arr, nan_row), axis=0)
        
        # Create GON feature DataFrame
        gon_feat_df = pd.DataFrame(joint_angle_arr, columns=COLS_JOINT_ANGLE)
        gon_feat_df[COLS_JOINT_VELOCITY] = joint_vel_arr
            
        # Save in a new file
        if "ModOut" in mod_out_file:
            gon_feat_file = mod_out_file.replace("ModOut", "GON_features")
        elif "ModelOutputs" in mod_out_file:
            gon_feat_file = mod_out_file.replace("ModelOutputs", "GON_features")
        else:
            raise ValueError("Something went wrong with the name of the mode outputs file.")
        gon_feat_df.to_csv(gon_feat_file, index=False, float_format="%.6g")
    return