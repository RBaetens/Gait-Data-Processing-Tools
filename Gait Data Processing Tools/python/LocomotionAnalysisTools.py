# Dependencies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

import os
import copy

from lat_globals import *
from lat_helpers import *
from lat_visuals import *
from lat_markers import *
from lat_feature import *


####################
# Class definition #
####################

class ProcessedData:
    def __init__(self, foldersToRead):
        initDictMean, initDictCounter = initialize_activity_dictionaries(foldersToRead)
        initDictSds, _ = initialize_activity_dictionaries(foldersToRead)
        initDictSkw, _ = initialize_activity_dictionaries(foldersToRead)
        self.meanDict = initDictMean
        self.sdsDict = initDictSds
        self.skwDict = initDictSkw
        self.counterDict = initDictCounter
        return

    def update_means_sds(self, filesToRead, instrumented_leg):
        self.meanDict, self.sdsDict, self.skwDict, self.counterDict = calculate_means_sds(filesToRead, self.meanDict, 
                                                                                          self.sdsDict, self.skwDict, 
                                                                                          self.counterDict, instrumented_leg)
        return

    def store_results(self, folder, style="shaded", activities_to_save=None):
        # Check which activities to save data for: default is all
        if activities_to_save:
            pass
        else:
            activities_to_save = self.meanDict.keys()

        # Save data
        for activity in activities_to_save:
            # Classic gait plot
            filename = folder + "/" + activity + " - Classical gait plot.png"
            classic_gait_plot(self, activity, style=style, filename=filename)

            filename = folder + "/" + activity + " - Classical gait plot filt.png"
            classic_gait_plot(self, activity, style=style, filtered=True, filename=filename)

            filename = folder + "/" + activity + " - Classical gait plot filt offset.png"
            classic_gait_plot(self, activity, style=style, filtered=True, offset=True, filename=filename)

            # EMG plot
            filename = folder + "/" + activity + " - EMG plot.png"
            emg_plot(self, activity, style=style, signal="ACT", filename=filename)

            filename = folder + "/" + activity + " - EMG plot line.png"
            emg_plot(self, activity, style="line", signal="ACT", filename=filename)

            filename = folder + "/" + activity + " - EMG plot RMS.png"
            emg_plot(self, activity, style=style, signal="RMS", filename=filename)

            filename = folder + "/" + activity + " - EMG plot RMS line.png"
            emg_plot(self, activity, style="line", signal="RMS", filename=filename)

            # GRF plot
            filename = folder + "/" + activity + " - GRF plot.png"
            grf_plot(self, activity, style=style, filename=filename)

            # CoP plot
            filename = folder + "/" + activity + " - CoP plot.png"
            cop_traj_plot(self, activity, filename=filename)

            # Numbers
            self.meanDict[activity].to_csv((folder + "/" + activity + " - means.csv"), index=False)
            self.sdsDict[activity].to_csv((folder + "/" + activity + " - sds.csv"), index=False)
            self.skwDict[activity].to_csv((folder + "/" + activity + " - skw.csv"), index=False)
            self.counterDict[activity].to_csv((folder + "/" + activity + " - cnt.csv"), index=False)
            
        return


##########################
# Manipulation functions #
##########################

# Rotate forces to align force plate x-axis with foot orientation
def rotate_grf(forcesXY, angle):
    rotation_matrix = np.array([[np.cos(-angle), -np.sin(-angle)],
                                [np.sin(-angle), np.cos(-angle)]])
    
    forcesXY = np.matmul(rotation_matrix, forcesXY.transpose())
    forcesXY = forcesXY.transpose()
    
    return forcesXY

# Calculate pelvic tilt based on Root segment orientation, input in radians, output in degrees
def get_pelvis_angles(pel_rot_vecs, f_sample=100): 
    # Constants
    up = np.array([0, 0, 1])
    to_deg = 180/np.pi
    
    # Change orientation representation, rotation vector to rotation matrices
    n_pts = np.shape(pel_rot_vecs)[0]
    rot_matrices = np.zeros((n_pts, 3, 3))
    z_axes = np.zeros((n_pts, 3))
    pelvis_angles = np.zeros(n_pts)

    for i in np.arange(n_pts):
        rotation = Rotation.from_rotvec(pel_rot_vecs[i, :])
        rot_matrices[i, :, :] = rotation.as_matrix()
        z_axes[i, :] = np.transpose(rot_matrices[i, 2, :])
        pelvis_angles[i] = np.arccos(np.dot(z_axes[i, :], up)) * to_deg - 90

    mean_angle = np.mean(pelvis_angles)
    pelvis_angles = -pelvis_angles + 2*mean_angle
    
    return pelvis_angles


###################################
# Reference trajectory generation #
###################################

# The segmentation process for cyclic locomotion
def segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, flag_mask, left_right, filt_freq=None):
    # Find foot strikes
    if left_right == "left":
        arrHEE = dfFootMarkers[["LHEEX", "LHEEY", "LHEEZ"]].to_numpy()
        arrTOE = dfFootMarkers[["LTOEX", "LTOEY", "LTOEZ"]].to_numpy()
        arrGRF = arrGRF[:, [0, 2, 4, 6]]
        arrCOP = dfCOP[["LCOPX_FP", "LCOPY_FP"]].to_numpy()
    else:
        arrHEE = dfFootMarkers[["RHEEX", "RHEEY", "RHEEZ"]].to_numpy()
        arrTOE = dfFootMarkers[["RTOEX", "RTOEY", "RTOEZ"]].to_numpy()
        arrGRF = arrGRF[:, [1, 3, 5, 7]]
        arrCOP = dfCOP[["RCOPX_FP", "RCOPY_FP"]].to_numpy()
    
    if activity == "Ramp ascent (10deg)":
        foot_strikes = get_toe_offs(arrTOE, filt_freq=filt_freq) # not actually foot_strikes in this case, but handy to use same variable
    else:
        foot_strikes = get_foot_strikes(arrHEE, arrTOE, arrGRF, filt_freq=filt_freq)

    # Store data between the first 2 consecutive foot strikes
    if len(foot_strikes) >= 2:
        
        # Correct foot_strikes
        if activity == "Ramp ascent (10deg)":
            over_corr = 20
            if left_right == "left":
                GRFZ = dfModelOutputs["LGroundReactionForceZ"].to_numpy()[foot_strikes[0]+over_corr:foot_strikes[1]+over_corr]
            else:
                GRFZ = dfModelOutputs["RGroundReactionForceZ"].to_numpy()[foot_strikes[0]+over_corr:foot_strikes[1]+over_corr]
            under_corr = np.argmax((np.flip(GRFZ) > 0.001))
            corr = over_corr - under_corr
            foot_strikes = np.array(foot_strikes) + corr
        elif ("(SBS)" in activity):
            pass # no correction because the when two feet are on the same force plate, the GRF suddenly switches from one foot to the other
        else:
            over_corr = 20
            if left_right == "left":
                GRFZ = dfModelOutputs["LGroundReactionForceZ"].to_numpy()[foot_strikes[0]-over_corr:foot_strikes[1]-over_corr]
            else:
                GRFZ = dfModelOutputs["RGroundReactionForceZ"].to_numpy()[foot_strikes[0]-over_corr:foot_strikes[1]-over_corr]
            under_corr = np.argmax((GRFZ > 0.001))
            corr = over_corr - under_corr
            foot_strikes = np.array(foot_strikes) - corr
        
        # Segment
        arrHEE = arrHEE[foot_strikes[0]:foot_strikes[1]]
        arrTOE = arrTOE[foot_strikes[0]:foot_strikes[1]]
        arrGRF = arrGRF[foot_strikes[0]:foot_strikes[1]]
        arrCOP = arrCOP[foot_strikes[0]:foot_strikes[1]]
        
        if left_right == "left":
            arrModelOutputs = dfModelOutputs[COLS_REFERENCE_DATA_LEFT].to_numpy()[foot_strikes[0]:foot_strikes[1]]
        else:
            arrModelOutputs = dfModelOutputs[COLS_REFERENCE_DATA_RIGHT].to_numpy()[foot_strikes[0]:foot_strikes[1]]
            
        # Align ground reaction forces correctly
        angle = get_foot_prog_angle(arrHEE, arrTOE, arrGRF)
        forcesXY = arrModelOutputs[:, FORCES_SLICE] # X and Y component of GRF
        if left_right == "left":
            forcesXY = rotate_grf(forcesXY, angle)
            forcesXY[:, 1] = -forcesXY[:, 1] # sign reversal is result of different frame of reference definition
            arrModelOutputs[:, FORCES_SLICE] = forcesXY
            arrModelOutputs[:, FRICTION_IND] = -arrModelOutputs[:, FRICTION_IND]
        else:
            arrModelOutputs[:, FORCES_SLICE] = rotate_grf(forcesXY, angle)

        # Set COP = np.nan in regions where there is no contact with the ground
        arrCOP = np.where(np.logical_and((arrCOP < 10e-5), (arrCOP > -10e-5)), np.nan, arrCOP)

        # Append EMG, COP and pelvis data to output array
        arrEMGFeat = dfEMGFeat.to_numpy()[foot_strikes[0]:foot_strikes[1]]
        if left_right == "left": # still hard coded that instrumented_leg is the right leg!
            arrEMGFeat[:] = np.nan
        arrPelvis = arrPelvis[foot_strikes[0]:foot_strikes[1]]
        arrModelOutputs = np.concatenate((arrModelOutputs, arrEMGFeat), axis=1)
        arrModelOutputs = np.concatenate((arrModelOutputs, arrCOP), axis=1)
        arrModelOutputs = np.concatenate((arrModelOutputs, arrPelvis), axis=1)

        # Apply flag mask
        arrModelOutputs = np.multiply(arrModelOutputs, flag_mask[foot_strikes[0]:foot_strikes[1]])
        
        # Resample
        arrModelOutputs = interp1d(arrModelOutputs, N_PTS_REPR)

    else:
        arrModelOutputs = np.empty((N_PTS_REPR, N_COLS_REFERENCE_DATA))
        arrModelOutputs[:] = np.nan
        
    return arrModelOutputs

# Code for the actual update shared by all activity types
def update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs):
    # Add data to dictionary where it's not None and keep track with counter
    nan_mask = np.logical_not(np.isnan(arrModelOutputs))
    meanDict[activity] = np.where(nan_mask, np.add(meanDict[activity], arrModelOutputs), meanDict[activity])
    sdsDict[activity] = np.where(nan_mask, np.add(sdsDict[activity], np.power(arrModelOutputs, 2)), sdsDict[activity])
    skwDict[activity] = np.where(nan_mask, np.add(skwDict[activity], np.power(arrModelOutputs, 3)), skwDict[activity])
    counterDict[activity] = np.add(counterDict[activity], nan_mask)
    
    return meanDict, sdsDict, skwDict, counterDict

# Add data from one file to the activity dictionaries for a static activity
def update_activity_dictionaries_SA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                    arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                    flag_mask, instrumented_leg, filt_freq=None):
    
    # Take average and resample to the right size to store in the dictionary
    
    if instrumented_leg == "left":
        # Get left foot data
        arrHEE = dfFootMarkers[["LHEEX", "LHEEY", "LHEEZ"]].to_numpy()
        arrTOE = dfFootMarkers[["LTOEX", "LTOEY", "LTOEZ"]].to_numpy()
        arrModelOutputs = dfModelOutputs[COLS_REFERENCE_DATA_LEFT].to_numpy()
        arrGRF = arrGRF[:, [0, 2, 4, 6]]
        arrCOP = dfCOP[["LCOPX_FP", "LCOPY_FP"]].to_numpy()

        # Align ground reaction forces correctly
        angle = get_foot_prog_angle(arrHEE, arrTOE, arrGRF)
        forcesXY = arrModelOutputs[:, FORCES_SLICE] # X and Y component of GRF
        forcesXY = rotate_grf(forcesXY, angle)
        forcesXY[:, 1] = -forcesXY[:, 1] # sign reversal is result of different frame of reference definition
        arrModelOutputs[:, FORCES_SLICE] = forcesXY
        arrModelOutputs[:, FRICTION_IND] = -arrModelOutputs[:, FRICTION_IND]
            
    elif instrumented_leg == "right":
        # Get right foot data
        arrHEE = dfFootMarkers[["RHEEX", "RHEEY", "RHEEZ"]].to_numpy()
        arrTOE = dfFootMarkers[["RTOEX", "RTOEY", "RTOEZ"]].to_numpy()
        arrModelOutputs = dfModelOutputs[COLS_REFERENCE_DATA_RIGHT].to_numpy()
        arrGRF = arrGRF[:, [1, 3, 5, 7]]
        arrCOP = dfCOP[["RCOPX_FP", "RCOPY_FP"]].to_numpy()

        # Align ground reaction forces correctly
        angle = get_foot_prog_angle(arrHEE, arrTOE, arrGRF)
        forcesXY = arrModelOutputs[:, FORCES_SLICE] # X and Y component of GRF
        arrModelOutputs[:, FORCES_SLICE] = rotate_grf(forcesXY, angle)
            
    else:
        raise ValueError("Please define 'instrumented_leg' as 'left' or 'right'")

    # Append EMG, COP and pelvis data to output array
    arrEMGFeat = dfEMGFeat.to_numpy()
    arrModelOutputs = np.concatenate((arrModelOutputs, arrEMGFeat), axis=1)
    arrModelOutputs = np.concatenate((arrModelOutputs, arrCOP), axis=1)
    arrModelOutputs = np.concatenate((arrModelOutputs, arrPelvis), axis=1)

    # Apply flag mask
    arrModelOutputs = np.multiply(arrModelOutputs, flag_mask)
    
    # Get mean and 'resample'
    arrModelOutputs = np.nanmean(arrModelOutputs, axis=0)
    arrModelOutputs = np.tile(arrModelOutputs, (N_PTS_REPR, 1))

    # Add data to dictionary where it's not None and keep track with counter
    meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
    
    return meanDict, sdsDict, skwDict, counterDict

# Add data from one file to the activity dictionaries for a transitional activity
def update_activity_dictionaries_TA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                    arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                    flag_mask, instrumented_leg, filt_freq=None):
    
    # Define the included activities
    if activity == "Calf raises":
        first_act = "Calf raise up"
        last_act = "Calf raise down"
        
    elif activity == "Sit-to-stand and stand-to-sit":
        first_act = "Sit-to-stand"
        last_act = "Stand-to-sit"
        
    elif activity == "Squat":
        first_act = "Squat down"
        last_act = "Squat up"
        
    else:
        raise ValueError("Something went wrong defining the transitional activities")

    # Sort out left vs. right
    if instrumented_leg == "left":
        # Get left foot data
        arrHEE = dfFootMarkers[["LHEEX", "LHEEY", "LHEEZ"]].to_numpy()
        arrTOE = dfFootMarkers[["LTOEX", "LTOEY", "LTOEZ"]].to_numpy()
        arrModelOutputs = dfModelOutputs[COLS_REFERENCE_DATA_LEFT].to_numpy()
        arrPowers = dfModelOutputs[["LAnklePowerZ", "LKneePowerZ", "LHipPowerZ"]].to_numpy()
        arrGRF = arrGRF[:, [0, 2, 4, 6]]
        arrCOP = dfCOP[["LCOPX_FP", "LCOPY_FP"]].to_numpy()
            
    elif instrumented_leg == "right":
        # Get right foot data
        arrHEE = dfFootMarkers[["RHEEX", "RHEEY", "RHEEZ"]].to_numpy()
        arrTOE = dfFootMarkers[["RTOEX", "RTOEY", "RTOEZ"]].to_numpy()
        arrModelOutputs = dfModelOutputs[COLS_REFERENCE_DATA_RIGHT].to_numpy()
        arrPowers = dfModelOutputs[["RAnklePowerZ", "RKneePowerZ", "RHipPowerZ"]].to_numpy()
        arrGRF = arrGRF[:, [1, 3, 5, 7]]
        arrCOP = dfCOP[["RCOPX_FP", "RCOPY_FP"]].to_numpy()
            
    else:
        raise ValueError("Please define 'instrumented_leg' as 'left' or 'right'")
    
    # Find stationary points to base segmentation on
    # Start by looking at the velocity of the modeled bones in the lower body
    bones = dfModelOutputs[COLS_BONES].to_numpy()
    bones = np.power(np.diff(bones, axis=0)*100, 2)
    n_bones = int(np.round(np.shape(bones)[1]/3))
    n_samples = np.shape(arrGRF)[0]
    bone_speeds = np.zeros((n_samples-1, n_bones))
    for j in np.arange(n_bones):
        bone_speeds[:, j] = np.sqrt(np.sum(bones[:, 3*j:3*(j+1)], axis=1))
    bone_speeds = fill_filter_robust(bone_speeds, 100, 3)
    avg_bone_speed = np.nanmean(bone_speeds, axis=1)

    # Find where the velocity is the lowest
    bone_speed_peaks, _ = find_peaks(-avg_bone_speed, height=-150, distance=50, prominence=40, width=20, rel_height=0.5)
    max_inter_peak_dist = np.max(np.diff(bone_speed_peaks))
    first = max((bone_speed_peaks[0] - max_inter_peak_dist), 0)
    last = min((bone_speed_peaks[-1] + max_inter_peak_dist), n_samples-2)
    segment_points = np.concatenate(([first], bone_speed_peaks, [last])).tolist()

    # Manual correction for some errors
    if (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [45, 280, 429, 580, 722, 867, 1022, 1160, 1239, 1311, 1452, 1687, 1807]):
        segment_points = [45, 280, 429, 580, 722, 867, 1022, 1160, 1311, 1452, 1687]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [59, 362, 539, 710, 908, 1095, 1283, 1490, 1592, 1739, 2042, 2345]):
        segment_points = [59, 362, 539, 710, 908, 1095, 1283, 1490, 1739, 2042, 2345]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [93, 329, 565, 722, 916, 1000, 1098, 1294, 1459, 1641, 1802, 2038]):
        segment_points = [93, 329, 565, 722, 916, 1098, 1294, 1459, 1641, 1802, 2038]
    elif (activity == "Squat") and (segment_points == [119, 338, 395, 500, 719, 844, 1055, 1183, 1384, 1531, 1722, 1941]):
        segment_points = [119, 338, 500, 719, 844, 1055, 1183, 1384, 1531, 1722, 1941]
    elif (activity == "Calf raises") and (segment_points == [102, 205, 304, 383, 481, 569, 665, 768, 859, 962]):
        segment_points = [102, 205, 304, 383, 481, 569, 665, 768, 859, 962, 1100]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [73, 300, 505, 674, 901, 1045, 1250, 1427, 1619, 1774, 1972, 2168, 2395]):
        segment_points = [73, 300, 505, 674, 901, 1045, 1250, 1427, 1619, 1774, 1972, 2168, 2600]

    # Refine this first guess by looking at joint powers
    segment_points_fine = []
    joint_power_thresh = 0.05
    for i in np.arange(len(segment_points)-1):
        arrPowers_seg = arrPowers[segment_points[i]:segment_points[i+1]]
        start_seg = segment_points[i]
        end_seg = segment_points[i+1]
        
        # Forward sweep
        j = 0
        while (j < (end_seg-start_seg)) and ((np.max(np.abs(arrPowers_seg[j, :])) < joint_power_thresh) or (np.sum(np.isnan(arrPowers_seg[j, :])) == 3)):
            j += 1
        if j != (end_seg-start_seg):
            start_seg += j
        
        # Backward sweep
        j = 1
        while (j < (end_seg-start_seg)) and ((np.max(np.abs(arrPowers_seg[-j, :])) < joint_power_thresh) or (np.sum(np.isnan(arrPowers_seg[-j, :])) == 3)):
            j += 1
        if j != (end_seg-start_seg):
            end_seg = end_seg - (j-1)

        # Update
        segment_points_fine.append(start_seg)
        segment_points_fine.append(end_seg)

    # Manual correction for some errors
    if (activity == "Squat") and (segment_points == [50, 264, 396, 524, 670, 796, 968, 1105, 1246, 1378, 1555, 1669, 1832, 1973, 2112, 2277, 2491, 2604, 2818]):
        segment_points_fine = [157, 264, 264, 377, 412, 524, 524, 629, 693, 796, 796, 910, 1004, 1105, 1105, 1217, 1270, 1378, 1379, 1493, 1568, 1669, 1669, 1784, 1868, 1973, 1973, 2088, 2161, 2277, 2279, 2390, 2510, 2604, 2604, 2726]
    elif (activity == "Squat") and (2005 in segment_points):
        segment_points.remove(2005)
        segment_points_fine.remove(2005)
        segment_points_fine.remove(2005)
    elif (activity == "Sit-to-stand and stand-to-sit") and (2695 in segment_points):
        segment_points.remove(2695)
        segment_points_fine.remove(2695)
        segment_points_fine.remove(2640)
    elif (activity == "Calf raises") and (1503 in segment_points):
        segment_points.remove(1503)
        segment_points_fine.remove(1503)
        segment_points_fine.remove(1347)
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [0, 321, 794, 1268, 1756, 2285, 2667]):
        segment_points = [84, 321, 558, 794, 1031, 1268, 1512, 1756, 2021, 2285, 2549]
        segment_points_fine = [151, 307, 336, 530, 622, 766, 849, 1003, 1098, 1255, 1297, 1484, 1580, 1741, 1799, 1993, 2089, 2249, 2309, 2521]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [67, 286, 505, 650, 842, 997, 1184, 1361, 1550, 1728, 1900, 2068, 2263, 2421, 2640]):
        segment_points_fine = [138, 274, 300, 505, 505, 633, 674, 789, 848, 993, 1025, 1158, 1213, 1343, 1388, 1524, 1574, 1710, 1741, 1868, 1923, 2050, 2085, 2210, 2288, 2395, 2436, 2549]
    elif (activity == "Squat") and (segment_points == [137, 251, 365, 455, 564, 650, 762, 846, 960, 1048, 1161, 1246, 1360]):
        segment_points_fine = [137, 264, 280, 360, 366, 468, 479, 559, 564, 665, 679, 757, 762, 861, 866, 955, 965, 1063, 1068, 1156, 1162, 1259, 1280, 1358]
    elif (activity == "Calf raises") and (segment_points == [101, 242, 383, 516, 638, 769, 899, 1035, 1157, 1292, 1433]):
        segment_points_fine = [165, 215, 249, 336, 436, 487, 524, 597, 696, 744, 769, 843, 958, 1001, 1061, 1120, 1222, 1259, 1292, 1369]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [136, 318, 499, 636, 796, 926, 1108, 1252, 1425, 1570, 1752]):
        segment_points_fine = [193, 308, 321, 475, 510, 627, 660, 770, 807, 926, 926, 1091, 1115, 1232, 1259, 1407, 1436, 1553, 1583, 1706]
    elif (activity == "Squat") and (segment_points == [266, 390, 467, 591, 683, 764, 863, 978, 1076, 1168, 1273, 1371, 1495]):
        segment_points_fine = [266, 367, 387, 454, 494, 574, 587, 683, 683, 770, 782, 853, 894, 978, 978, 1056, 1092, 1174, 1178, 1257, 1292, 1377, 1388, 1457]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [46, 328, 610, 790, 1056, 1260, 1518, 1719, 1992, 2206, 2482, 2687, 2947, 3196, 3478]):
        segment_points = [610, 790, 1056, 1260, 1518, 1719, 1992, 2206, 2482, 2687, 2947]
        segment_points_fine = [610, 789, 794,   1016, 1056, 1248, 1279,   1478, 1518, 1708, 1736, 1957, 1996, 2197, 2212, 2398, 2484, 2677, 2699,   2914]
    elif (activity == "Squat") and (segment_points == [97, 254, 404, 514, 665, 775, 931, 1041, 1198, 1317, 1463, 1576, 1733]):
        segment_points_fine = [162, 284, 284, 391, 421, 544, 544, 655, 682, 805, 805, 920, 945, 1071, 1071, 1191, 1207, 1347, 1347, 1459, 1475, 1606, 1606, 1729]
    elif (activity == "Calf raises") and (segment_points == [109, 252, 395, 506, 632, 767, 887, 999, 1136, 1269, 1412]):
        segment_points_fine = [148, 252, 252, 365, 429, 506, 506, 610, 681, 767, 767, 856, 922, 992, 999, 1102, 1184, 1269, 1269, 1392]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [77, 224, 344, 445, 592, 690, 809, 905, 1016, 1126, 1273]):
        segment_points_fine = [116, 214, 229, 326, 348, 432, 455, 582, 592, 685, 690, 799, 809, 900, 905, 1000, 1021, 1121, 1126, 1240]
    elif (activity == "Calf raises") and (segment_points == [127, 209, 291, 356, 437, 515, 595, 667, 741, 815, 897]):
        segment_points_fine = [157, 204, 209, 277, 308, 356, 355, 413, 467, 512, 519, 571, 621, 667, 667, 722, 769, 815, 815, 873]
    elif (activity == "Squat") and (segment_points == [110, 213, 303, 406, 492, 591, 680, 774, 860, 953, 1056]):
        segment_points_fine = [120, 216, 216, 300, 310, 409, 409, 488, 497, 594, 594, 679, 687, 777, 777, 855, 865, 956, 956, 1033]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [172, 480, 788, 933, 1114, 1266, 1413, 1601, 1786, 1949, 2146, 2319, 2627]):
        segment_points_fine = [252, 376, 567, 714, 810, 918, 955, 1064, 1128, 1257, 1286, 1396, 1465, 1584, 1628, 1750, 1809, 1937, 1963, 2089, 2167, 2314, 2329, 2453]
    elif (activity == "Squat") and (segment_points == [171, 281, 382, 492, 592, 698, 795, 901, 1006, 1110, 1220]):
        segment_points_fine = [193, 283, 283, 370, 387, 494, 494, 592, 592, 700, 700, 788, 800, 903, 903, 1006, 1006, 1112, 1112, 1220]
    elif (activity == "Calf raises") and (segment_points == [64, 234, 404, 554, 691, 833, 980, 1125, 1282, 1440, 1610]):
        segment_points_fine = [172, 234, 234, 349, 461, 554, 554, 645, 731, 827, 846, 974, 1026, 1125, 1125, 1234, 1348, 1432, 1446, 1542]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [45, 280, 429, 580, 722, 867, 1022, 1160, 1311, 1452, 1687]):
        segment_points_fine = [137, 280, 284, 403, 440, 579, 585, 716, 728, 867, 867, 1007, 1027, 1160, 1160, 1306, 1318, 1452, 1452, 1580]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [59, 362, 539, 710, 908, 1095, 1283, 1490, 1739, 2042, 2345]):
        segment_points_fine = [123, 272, 369, 498, 562, 703, 714, 865, 924, 1087, 1110, 1237, 1310, 1485, 1490, 1709, 1739, 1911, 2087, 2237]
    elif (activity == "Calf raises") and (segment_points == [135, 220, 302, 371, 456, 531, 602, 684, 750, 824, 909]):
        segment_points_fine = [168, 215, 220, 281, 320, 366, 371, 436, 476, 521, 531, 580, 619, 674, 684, 725, 767, 819, 824, 874]
    elif (activity == "Calf raises") and (segment_points == [0, 207, 344, 467, 615, 832, 965, 1099, 1250, 1411, 1682, 1880, 1995]):
        segment_points_fine = [123, 197, 217, 285, 373, 453, 520, 593, 663, 725, 858, 904, 1002, 1062, 1128, 1186, 1300, 1362, 1433, 1493, 1789, 1840, 1890, 1932]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [112, 328, 544, 706, 922, 1082, 1247, 1448, 1632, 1819, 2021, 2196, 2412]):
        segment_points_fine = [166, 315, 345, 482, 552, 696, 739, 881, 926, 1064, 1086, 1247, 1297, 1439, 1465, 1617, 1663, 1796, 1850, 1981, 2021, 2183, 2240, 2357]
    elif (activity == "Calf raises") and (segment_points == [121, 267, 400, 517, 663, 778, 876, 1020, 1134, 1246, 1392]):
        segment_points_fine = [196, 251, 299, 337, 451, 510, 554, 600, 718, 768, 796, 851, 948, 999, 1021, 1075, 1185, 1236, 1257, 1311]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [73, 300, 505, 674, 901, 1045, 1250, 1427, 1619, 1774, 1972, 2168, 2600]):
        segment_points_fine = [154, 285, 317, 464, 544, 662, 702, 848, 910, 1027, 1076, 1207, 1258, 1386, 1453, 1604, 1627, 1765, 1796, 1939, 2004, 2141, 2361, 2537]
    elif (activity == "Calf raises") and (segment_points == [102, 205, 304, 383, 481, 569, 665, 768, 859, 962, 1100]):
        segment_points_fine = [155, 202, 202, 262, 342, 380, 380, 447, 525, 564, 566, 623, 716, 765, 765, 822, 904, 958, 960, 1015]
    elif (activity == "Squat") and (segment_points == [152, 275, 386, 490, 597, 690, 813, 923, 1032, 1138, 1261]):
        segment_points_fine = [170, 278, 278, 382, 392, 493, 493, 596, 605, 693, 693, 811, 822, 926, 926, 1032, 1032, 1151, 1151, 1253]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [30, 312, 594, 768, 1011, 1214, 1436, 1682, 1939, 2186, 2468]):
        segment_points_fine = [153, 288, 373, 476, 616, 740, 804, 933, 1041, 1195, 1248, 1400, 1500, 1661, 1706, 1862, 1987, 2152, 2239, 2375]
    elif (activity == "Calf raises") and (segment_points == [8, 147, 251, 366, 495, 632, 771, 880, 1005, 1128, 1267]):
        segment_points_fine = [96, 145, 147, 195, 301, 366, 366, 454, 552, 629, 636, 694, 818, 880, 880, 950, 1065, 1120, 1134, 1200]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [81, 281, 461, 620, 820, 965, 1157, 1295, 1495]):
        segment_points_fine = [145, 277, 287, 404, 467, 612, 633, 765, 820, 962, 974, 1097, 1157, 1288, 1305, 1434]
    elif (activity == "Squat") and (segment_points == [113, 229, 345, 436, 547, 648, 745, 848, 954, 1056, 1172]):
        segment_points_fine = [132, 230, 230, 338, 350, 437, 437, 539, 555, 649, 649, 745, 745, 849, 849, 951, 958, 1057, 1057, 1165]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [42, 226, 410, 531, 686, 801, 965, 1094, 1215, 1358, 1542]):
        segment_points_fine = [102, 222, 228, 387, 426, 528, 536, 673, 706, 796, 810, 934, 973, 1078, 1098, 1211, 1241, 1336, 1363, 1490]
    elif (activity == "Calf raises") and (segment_points == [97, 221, 337, 416, 525, 622, 740, 846, 970, 1075, 1199]):
        segment_points_fine = [158, 221, 221, 291, 368, 416, 416, 481, 559, 618, 631, 697, 777, 833, 848, 915, 1006, 1064, 1080, 1150]
    
    elif (activity == "Calf raises") and (segment_points == [0, 268, 471, 620, 841, 1072, 1270, 1526, 1787, 2031, 2315, 2556, 2812, 3081, 3319, 3489, 3757, 3957, 4203]):
        segment_points_fine = [160, 255, 301, 400, 519, 612, 642, 762, 918, 1017, 1078, 1163, 1358, 1437, 1606, 1708, 1876, 1964, 2081, 2190, 2399, 2483, 2621, 2727, 2915, 3031, 3098, 3187, 3386, 3470, 3563, 3682, 3831, 3926, 3994, 4104]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [0, 432, 693, 926, 1166, 1388, 1859, 2213, 2522, 2771, 3187, 3438, 3664, 4116, 4591, 4801, 5208, 5604, 5870, 6345]):
        segment_points = [0, 432, 693, 926, 1166, 1388, 1859, 2213, 2522, 2771, 3187, 3664, 4116, 4591, 4801, 5208, 5604, 5870, 6345]
        segment_points_fine = [253, 404, 480, 658, 727, 865, 981, 1142, 1197, 1359, 1439, 1601, 1929, 2083, 2270, 2469, 2560, 2719, 2887, 3117, 3232, 3410, 3709, 3887, 4212, 4390, 4608, 4768, 4857, 5031, 5383, 5542, 5658, 5839, 5976, 6151]
        

    # # If something goes wrong, uncomment this plot to troubleshoot
    # print(segment_points)
    # print(segment_points_fine)
    # plt.figure()
    # plt.title(activity)
    # plt.plot(avg_bone_speed, label="avg_bone_speed")
    # plt.scatter(segment_points, avg_bone_speed[segment_points], label="segment_points")
    # plt.scatter(segment_points_fine, avg_bone_speed[segment_points_fine], marker="+", label="segment_points_fine")
    # plt.legend()
    # plt.show()

    arrEMGFeat = dfEMGFeat.to_numpy()
    
    # Segment and update
    for i in np.arange(int(np.round(len(segment_points_fine)/2))):
        # Label correctly
        if (i%2) == 0:
            activity = first_act
        else:
            activity = last_act

        # Segment
        arrHEE_seg = arrHEE[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        arrTOE_seg = arrTOE[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        arrModelOutputs_seg = arrModelOutputs[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        arrGRF_seg = arrGRF[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        flag_mask_seg = flag_mask[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        arrEMGFeat_seg = arrEMGFeat[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        arrCOP_seg = arrCOP[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        arrPelvis_seg = arrPelvis[segment_points_fine[2*i]:segment_points_fine[2*i+1]]
        
        if instrumented_leg == "left":
            # Align ground reaction forces correctly
            angle = get_foot_prog_angle(arrHEE_seg, arrTOE_seg, arrGRF_seg)
            forcesXY = arrModelOutputs_seg[:, FORCES_SLICE] # X and Y component of GRF
            forcesXY = rotate_grf(forcesXY, angle)
            forcesXY[:, 1] = -forcesXY[:, 1] # sign reversal is result of different frame of reference definition
            arrModelOutputs_seg[:, FORCES_SLICE] = forcesXY
            arrModelOutputs_seg[:, FRICTION_IND] = -arrModelOutputs_seg[:, FRICTION_IND]
                
        else:    
            # Align ground reaction forces correctly
            angle = get_foot_prog_angle(arrHEE_seg, arrTOE_seg, arrGRF_seg)
            forcesXY = arrModelOutputs_seg[:, FORCES_SLICE] # X and Y component of GRF
            arrModelOutputs_seg[:, FORCES_SLICE] = rotate_grf(forcesXY, angle)

        # Append EMG, COP and pelvis data to output array
        arrModelOutputs_seg = np.concatenate((arrModelOutputs_seg, arrEMGFeat_seg), axis=1)
        arrModelOutputs_seg = np.concatenate((arrModelOutputs_seg, arrCOP_seg), axis=1)
        arrModelOutputs_seg = np.concatenate((arrModelOutputs_seg, arrPelvis_seg), axis=1)
        
        # Apply flag mask
        arrModelOutputs_seg = np.multiply(arrModelOutputs_seg, flag_mask_seg)
        
        # Resample
        arrModelOutputs_seg = interp1d(arrModelOutputs_seg, N_PTS_REPR)

        # Add data to dictionary where it's not None and keep track with counter
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs_seg)
    
    return meanDict, sdsDict, skwDict, counterDict

# Add data from one file to the activity dictionaries for a symmetric rhythmic activity
def update_activity_dictionaries_SRA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                     arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                     flag_mask, left_right, filt_freq=None):

    # Different names for stair ascent/descent
    if (activity == "Stair ascent (SOS)(L first)") or (activity == "Stair ascent (SOS)(R first)"):
        activity = "Stair ascent (SOS)"
    elif (activity == "Stair descent (SOS)(L first)") or (activity == "Stair descent (SOS)(R first)"):
        activity = "Stair descent (SOS)"
    
    # Perform segmentation = to determine for what range of samples to save the data
    # Then select the right data and resample it to the right size to store in the dictionary
    arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                     flag_mask, left_right, filt_freq=filt_freq)

    # Add data to dictionary where it's not None and keep track with counter
    meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
    
    return meanDict, sdsDict, skwDict, counterDict

# Add data from one file to the activity dictionaries for an asymmetric rhythmic activity
def update_activity_dictionaries_ARA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                     arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                     flag_mask, flag_mask_sbs, left_right, filt_freq=None):

    # Redefine activities and assign corresponding leg, then update as in SRA
    # Save data of both feet for side stepping and SBS stair ascent/descent, as both feet will land on the force plate(s)
    # Save data for the one foot on the force plate for circular walking
    
    if activity == "Side stepping (L)":
        activity = "Side stepping (leading leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask, "left", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
        activity = "Side stepping (following leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask, "right", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
    
    elif activity == "Side stepping (R)":
        activity = "Side stepping (leading leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask, "right", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
        activity = "Side stepping (following leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask, "left", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
    elif activity == "Stair ascent (SBS)(L first)":
        activity = "Stair ascent (SBS)(leading leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "left", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
        activity = "Stair ascent (SBS)(following leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "right", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
    elif activity == "Stair ascent (SBS)(R first)":
        activity = "Stair ascent (SBS)(leading leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "right", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
        activity = "Stair ascent (SBS)(following leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "left", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
    elif activity == "Stair descent (SBS)(L first)":
        activity = "Stair descent (SBS)(leading leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "left", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
        activity = "Stair descent (SBS)(following leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "right", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
    elif activity == "Stair descent (SBS)(R first)":
        activity = "Stair descent (SBS)(leading leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "right", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
        activity = "Stair descent (SBS)(following leg)"
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask_sbs, "left", filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
    elif activity == "Circular walking (diam 2m)(CW)":
        # Redefine the activity using which foot was on the force plate
        if left_right == "left":
            activity = "Circular walking (outer leg)"
        else:
            activity = "Circular walking (inner leg)"
            
        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask, left_right, filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
    elif activity == "Circular walking (diam 2m)(CCW)":
        # Redefine the activity using which foot was on the force plate
        if left_right == "left":
            activity = "Circular walking (inner leg)"
        else:
            activity = "Circular walking (outer leg)"

        arrModelOutputs = segment_cyclic(activity, arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                         flag_mask, left_right, filt_freq=filt_freq)
        meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries(meanDict, sdsDict, skwDict, counterDict, activity, arrModelOutputs)
        
    #else:
    #    raise ValueError("Something went wrong defining the asymmetric rhythmic activities")
    
    return meanDict, sdsDict, skwDict, counterDict

# Calculate the means of all the model outputs for each locomotion activity
def calculate_means_sds(filesToRead, meanDict, sdsDict, skwDict, counterDict, instrumented_leg):
    # Define filter frequency for marker stuff
    filt_freq = 25
    
    # Loop over all files to read, add model outputs together, keep track of number of non-NaN values (per output per gait phase)
    for i, file in enumerate(filesToRead):
        # Get all five files for one step and the corresponding activity
        activity, fileGRF, fileFootMarkers, fileModelOutputs, fileEMGFeat, fileGlobalAngles = file
        dfGRF = pd.read_csv(fileGRF)
        dfFootMarkers = pd.read_csv(fileFootMarkers)
        dfModelOutputs = pd.read_csv(fileModelOutputs)
        dfEMGFeat = pd.read_csv(fileEMGFeat)
        dfCOP = dfGRF[["LCOPX_FP", "LCOPY_FP", "RCOPX_FP", "RCOPY_FP"]]
        dfPelvis = pd.read_csv(fileGlobalAngles)[["Root_RX", "Root_RY", "Root_RZ"]]

        # Calculate pelvic tilt
        pel_rot_vecs = dfPelvis.to_numpy()
        pel_rot_vecs = pel_rot_vecs * (np.pi/180) # to radians
        arrPelvis = get_pelvis_angles(pel_rot_vecs, f_sample=100).reshape((-1, 1))
    
        # Find out which foot stepped on the force plate
        arrGRF = dfGRF[["LEFTFOOT_FP1", "RIGHTFOOT_FP1", "LEFTFOOT_FP2", "RIGHTFOOT_FP2", 
                        "LEFTFOOT_FP3", "RIGHTFOOT_FP3", "LEFTFOOT_FP4", "RIGHTFOOT_FP4"]].to_numpy()
        n_samples = np.shape(arrGRF)[0]
        n_left = np.sum(arrGRF[:, [0, 2, 4, 6]])
        n_right = np.sum(arrGRF[:, [1, 3, 5, 7]])
        if n_left >= n_right:
            left_right = "left"
        else:
            left_right = "right"

        # Mask for step-by-step stair ascent and descent, don't look at "BAD_FP_MASK" and "BAD_STEP_MASK" because "bad steps" are necessary
        # Will make rows 'None' if "BAD_FEET_MASK" flag went off
        flag_mask_sbs = dfGRF[["BAD_FEET_MASK"]].to_numpy()
        flag_mask_sbs = np.tile(flag_mask_sbs.reshape((-1, 1)), (1, N_COLS_REFERENCE_DATA))
        flag_mask_sbs = np.where(flag_mask_sbs, flag_mask_sbs, None)
        flag_mask_sbs = np.array(flag_mask_sbs, dtype=np.float64)

        # Mask for other activities, makes force plate related columns 'None' at rows where "BAD_FP_MASK" and "BAD_STEP_MASK" went off
        flag_mask_upd1 = dfGRF[["BAD_FP_MASK", "BAD_STEP_MASK"]].to_numpy()
        flag_mask_upd1 = np.logical_not(np.all(flag_mask_upd1, axis=1))
        flag_mask_upd1 = np.tile(flag_mask_upd1.reshape((-1, 1)), (1, N_COLS_REFERENCE_DATA))
        
        flag_mask_upd2 = np.zeros(np.shape(flag_mask_sbs))
        for col in IDX_FP_RELATED_HEADERS:
            flag_mask_upd2[:, col] = np.ones(np.shape(flag_mask_upd2[:, col]))
            
        flag_mask_upd = np.logical_not(np.logical_and(flag_mask_upd1, flag_mask_upd2))
        flag_mask_upd = np.where(flag_mask_upd, flag_mask_upd, None)
        flag_mask_upd = np.array(flag_mask_upd, dtype=np.float64)

        flag_mask = np.multiply(flag_mask_sbs, flag_mask_upd)
    
        # Different segmentation strategy depending on locomotion category
        if activity in SA:
            meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries_SA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                                                             arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                                                             flag_mask, instrumented_leg, filt_freq=filt_freq)
                
        elif activity in TA:
            meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries_TA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                                                             arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                                                             flag_mask, instrumented_leg, filt_freq=filt_freq)
                
        elif activity in SRA:
            meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries_SRA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                                                              arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                                                              flag_mask, left_right, filt_freq=filt_freq)
                
        elif activity in ARA:
            meanDict, sdsDict, skwDict, counterDict = update_activity_dictionaries_ARA(meanDict, sdsDict, skwDict, counterDict, activity, 
                                                                              arrGRF, dfFootMarkers, dfModelOutputs, dfEMGFeat, dfCOP, arrPelvis, 
                                                                              flag_mask, flag_mask_sbs, left_right, filt_freq=filt_freq)
        
        else:
            raise ValueError("Unrecognized activity")

    # Actually calculate mean and sds
    for activity in meanDict.keys():    
        # Force plate measures from ModelOutputs are None where they should be zero --> correct this by making the counter at least 1
        counterDict[activity][:, ALL_FORCES_SLICE] = np.where(counterDict[activity][:, ALL_FORCES_SLICE] == 0, 1, counterDict[activity][:, ALL_FORCES_SLICE])
        counterDict[activity][:, COP_SLICE] = np.where(counterDict[activity][:, COP_SLICE] == 0, 1, counterDict[activity][:, COP_SLICE])
        
        # Return None if no samples were added for a given activity and gait phase
        zero_mask = np.logical_not(counterDict[activity] == 0)

        # Calculate mean, mean of squares and mean of cubes
        meanDict[activity] = np.where(zero_mask, np.divide(meanDict[activity], counterDict[activity]), np.nan)
        sdsDict[activity] = np.where(zero_mask, np.divide(sdsDict[activity], counterDict[activity]), np.nan)
        skwDict[activity] = np.where(zero_mask, np.divide(skwDict[activity], counterDict[activity]), np.nan)

        # Calculate sds and skew
        sdsDict[activity] = np.sqrt(np.subtract(sdsDict[activity], np.power(meanDict[activity], 2)))
        skwDict[activity] = (skwDict[activity] - 3*meanDict[activity]*(sdsDict[activity]**2) - meanDict[activity]**3)/(sdsDict[activity]**3)

        # Change dtype to pd.DataFrame
        meanDict[activity] = pd.DataFrame(meanDict[activity], columns=COLS_REFERENCE_DATA)
        sdsDict[activity] = pd.DataFrame(sdsDict[activity], columns=COLS_REFERENCE_DATA)
        skwDict[activity] = pd.DataFrame(skwDict[activity], columns=COLS_REFERENCE_DATA)
        counterDict[activity] = pd.DataFrame(counterDict[activity], columns=COLS_REFERENCE_DATA)
            
    return meanDict, sdsDict, skwDict, counterDict


############################
# Training data generation #
############################

# Append data to the training data file for static activities
def append_static_data(fileTrainingData, activity, arr_mot, arr_imu, arr_gon, arr_emg):
    # Construct "Label", "GaitPhase", "GaitSpeed" and "NewFile" column
    n_pts_data = np.shape(arr_mot)[0]
    labels_dummy = np.zeros((n_pts_data, 1))
    labels = [activity for i in np.arange(n_pts_data)]
    
    gait_phase = np.empty((n_pts_data, 1))
    gait_phase[:] = np.nan

    gait_speed = np.empty((n_pts_data, 1))
    gait_speed[:] = np.nan
    
    new_file = np.zeros((n_pts_data, 1))
    new_file[0] = 1
    
    # Construct dataframe and store
    arr_all_data = np.concatenate((labels_dummy, gait_phase, gait_speed, new_file, arr_mot, arr_imu, arr_gon, arr_emg), axis=1)
    df_all_data = pd.DataFrame(arr_all_data, columns=COLS_TRAINING_DATA)
    df_all_data["Label"] = labels
    df_all_data.to_csv(fileTrainingData, header=False, index=False, mode="a")
    return

# Append data to the training data file for transitional activities
def append_transitional_data(fileTrainingData, activity, arr_bones, arr_powers, arr_mot, arr_imu, arr_gon, arr_emg):
    # Define the included activities
    if activity == "Calf raises":
        first_act = "Calf raise up"
        last_act = "Calf raise down"
        
    elif activity == "Sit-to-stand and stand-to-sit":
        first_act = "Sit-to-stand"
        last_act = "Stand-to-sit"
        
    elif activity == "Squat":
        first_act = "Squat down"
        last_act = "Squat up"

    elif activity == "Calf raises weight":
        first_act = "Calf raise up weight"
        last_act = "Calf raise down weight"
        
    elif activity == "Sit-to-stand and stand-to-sit weight":
        first_act = "Sit-to-stand weight"
        last_act = "Stand-to-sit weight"
        
    elif activity == "Squat weight":
        first_act = "Squat down weight"
        last_act = "Squat up weight"

    elif activity == "Calf raises obstacle":
        first_act = "Calf raise up obstacle"
        last_act = "Calf raise down obstacle"
        
    elif activity == "Sit-to-stand and stand-to-sit obstacle":
        first_act = "Sit-to-stand obstacle"
        last_act = "Stand-to-sit obstacle"
        
    elif activity == "Squat obstacle":
        first_act = "Squat down obstacle"
        last_act = "Squat up obstacle"
        
    else:
        raise ValueError("Something went wrong defining the transitional activities")
    
    # Find stationary points to base segmentation on
    
    # Start by looking at the velocity of the modeled bones in the lower body
    arr_bones = np.power(np.diff(arr_bones, axis=0)*100, 2)
    n_bones = int(np.round(np.shape(arr_bones)[1]/3))
    n_samples = np.shape(arr_mot)[0]
    bone_speeds = np.zeros((n_samples-1, n_bones))
    for j in np.arange(n_bones):
        bone_speeds[:, j] = np.sqrt(np.sum(arr_bones[:, 3*j:3*(j+1)], axis=1))
    bone_speeds = fill_filter_robust(bone_speeds, 100, 3)
    avg_bone_speed = np.nanmean(bone_speeds, axis=1)

    # Find where the velocity is the lowest
    bone_speed_peaks, _ = find_peaks(-avg_bone_speed, height=-150, distance=50, prominence=40, width=20, rel_height=0.5)
    if len(bone_speed_peaks) > 1:
        max_inter_peak_dist = np.max(np.diff(bone_speed_peaks))
        first = max((bone_speed_peaks[0] - max_inter_peak_dist), 0)
        last = min((bone_speed_peaks[-1] + max_inter_peak_dist), n_samples-2)
    else:
        first = 0
        last = len(avg_bone_speed)-1
    segment_points = np.concatenate(([first], bone_speed_peaks, [last])).tolist()

    # Manual correction for some errors
    if (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [25, 260, 409, 560, 702, 847, 1002, 1140, 1219, 1291, 1432, 1667, 1787]):
        segment_points = [25, 260, 409, 560, 702, 847, 1002, 1140, 1291, 1432, 1667]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [39, 342, 519, 690, 888, 1075, 1263, 1470, 1572, 1719, 2022, 2325]):
        segment_points = [39, 342, 519, 690, 888, 1075, 1263, 1470, 1719, 2022, 2325]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [73, 309, 545, 702, 896, 980, 1078, 1274, 1439, 1621, 1782, 2018]):
        segment_points = [73, 309, 545, 702, 896, 1078, 1274, 1439, 1621, 1782, 2018]
    elif (activity == "Squat") and (segment_points == [99, 318, 375, 480, 699, 824, 1035, 1163, 1364, 1511, 1702, 1921]):
        segment_points = [99, 318, 480, 699, 824, 1035, 1163, 1364, 1511, 1702, 1921]
    elif (activity == "Calf raises") and (segment_points == [82, 185, 284, 363, 461, 549, 645, 748, 839, 942]):
        segment_points = [82, 185, 284, 363, 461, 549, 645, 748, 839, 942, 1080]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [53, 280, 485, 654, 881, 1025, 1230, 1407, 1599, 1754, 1952, 2148, 2375]):
        segment_points = [53, 280, 485, 654, 881, 1025, 1230, 1407, 1599, 1754, 1952, 2148, 2580]


        
    # Refine this first guess by looking at joint powers
    segment_points_fine = []
    joint_power_thresh = 0.05
    for i in np.arange(len(segment_points)-1):
        arr_powers_seg = arr_powers[segment_points[i]:segment_points[i+1]]
        start_seg = segment_points[i]
        end_seg = segment_points[i+1]
        
        # Forward sweep
        j = 0
        while (j < (end_seg-start_seg)) and ((np.max(np.abs(arr_powers_seg[j, :])) < joint_power_thresh) or (np.sum(np.isnan(arr_powers_seg[j, :])) == 3)):
            j += 1
        if j != (end_seg-start_seg):
            start_seg += j
        
        # Backward sweep
        j = 1
        while (j < (end_seg-start_seg)) and ((np.max(np.abs(arr_powers_seg[-j, :])) < joint_power_thresh) or (np.sum(np.isnan(arr_powers_seg[-j, :])) == 3)):
            j += 1
        if j != (end_seg-start_seg):
            end_seg = end_seg - (j-1)

        # Update
        segment_points_fine.append(start_seg)
        segment_points_fine.append(end_seg)

    # Manual correction for some errors
    if (activity == "Squat") and (segment_points == [30, 244, 376, 504, 650, 776, 948, 1085, 1226, 1358, 1535, 1649, 1812, 1953, 2092, 2257, 2471, 2584, 2798]):
        segment_points_fine = [137, 244, 244, 357, 392, 504, 504, 609, 673, 776, 776, 890, 984, 1085, 1085, 1197, 1250, 1358, 1359, 1473, 1548, 1649, 1649, 1764, 1848, 1953, 1953, 2068, 2141, 2257, 2259, 2370, 2490, 2584, 2584, 2706]
    elif (activity == "Squat") and (1985 in segment_points):
        segment_points.remove(1985)
        segment_points_fine.remove(1985)
        segment_points_fine.remove(1985)
    elif (activity == "Sit-to-stand and stand-to-sit") and (2675 in segment_points):
        segment_points.remove(2675)
        segment_points_fine.remove(2675)
        segment_points_fine.remove(2620)
    elif (activity == "Calf raises") and (1483 in segment_points):
        segment_points.remove(1483)
        segment_points_fine.remove(1483)
        segment_points_fine.remove(1327)
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [0, 301, 774, 1248, 1736, 2265, 2647]):
        segment_points = [64, 301, 538, 774, 1011, 1248, 1492, 1736, 2001, 2265, 2529]
        segment_points_fine = [131, 287, 316, 510, 602, 746, 829, 983, 1078, 1235, 1277, 1464, 1560, 1721, 1779, 1973, 2069, 2229, 2289, 2501]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [47, 266, 485, 630, 822, 977, 1164, 1341, 1530, 1708, 1880, 2048, 2243, 2401, 2620]):
        segment_points_fine = [118, 254, 280, 485, 485, 613, 654, 769, 828, 973, 1005, 1138, 1193, 1323, 1368, 1504, 1554, 1690, 1721, 1848, 1903, 2030, 2065, 2190, 2268, 2375, 2416, 2529]
    elif (activity == "Squat") and (segment_points == [117, 231, 345, 435, 544, 630, 742, 826, 940, 1028, 1141, 1226, 1340]):
        segment_points_fine = [117, 244, 260, 340, 346, 448, 459, 539, 544, 645, 659, 737, 742, 841, 846, 935, 945, 1043, 1048, 1136, 1142, 1239, 1260, 1338]
    elif (activity == "Calf raises") and (segment_points == [81, 222, 363, 496, 618, 749, 879, 1015, 1137, 1272, 1413]):
        segment_points_fine = [145, 195, 229, 316, 416, 467, 504, 577, 676, 724, 749, 823, 938, 981, 1041, 1100, 1202, 1239, 1272, 1349]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [116, 298, 479, 616, 776, 906, 1088, 1232, 1405, 1550, 1732]):
        segment_points_fine = [173, 288, 301, 455, 490, 607, 640, 750, 787, 906, 906, 1071, 1095, 1212, 1239, 1387, 1416, 1533, 1563, 1686]
    elif (activity == "Squat") and (segment_points == [246, 370, 447, 571, 663, 744, 843, 958, 1056, 1148, 1253, 1351, 1475]):
        segment_points_fine = [246, 347, 367, 434, 474, 554, 567, 663, 663, 750, 762, 833, 874, 958, 958, 1036, 1072, 1154, 1158, 1237, 1272, 1357, 1368, 1437]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [26, 308, 590, 770, 1036, 1240, 1498, 1699, 1972, 2186, 2462, 2667, 2927, 3176, 3458]):
        segment_points = [590, 770, 1036, 1240, 1498, 1699, 1972, 2186, 2462, 2667, 2927]
        segment_points_fine = [590, 769, 774, 996, 1036, 1228, 1259, 1458, 1498, 1688, 1716, 1937, 1976, 2177, 2192, 2378, 2464, 2657, 2679, 2894]
    elif (activity == "Squat") and (segment_points == [77, 234, 384, 494, 645, 755, 911, 1021, 1178, 1297, 1443, 1556, 1713]):
        segment_points_fine = [142, 264, 264, 371, 401, 524, 524, 635, 662, 785, 785, 900, 925, 1051, 1051, 1171, 1187, 1327, 1327, 1439, 1455, 1586, 1586, 1709]
    elif (activity == "Calf raises") and (segment_points == [89, 232, 375, 486, 612, 747, 867, 979, 1116, 1249, 1392]):
        segment_points_fine = [128, 232, 232, 345, 409, 486, 486, 590, 661, 747, 747, 836, 902, 972, 979, 1082, 1164, 1249, 1249, 1372]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [57, 204, 324, 425, 572, 670, 789, 885, 996, 1106, 1253]):
        segment_points_fine = [96, 194, 209, 306, 328, 412, 435, 562, 572, 665, 670, 779, 789, 880, 885, 980, 1001, 1101, 1106, 1220]
    elif (activity == "Calf raises") and (segment_points == [107, 189, 271, 336, 417, 495, 575, 647, 721, 795, 877]):
        segment_points_fine = [137, 184, 189, 257, 288, 336, 335, 393, 447, 492, 499, 551, 601, 647, 647, 702, 749, 795, 795, 853]
    elif (activity == "Squat") and (segment_points == [90, 193, 283, 386, 472, 571, 660, 754, 840, 933, 1036]):
        segment_points_fine = [100, 196, 196, 280, 290, 389, 389, 468, 477, 574, 574, 659, 667, 757, 757, 835, 845, 936, 936, 1013]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [152, 460, 768, 913, 1094, 1246, 1393, 1581, 1766, 1929, 2126, 2299, 2607]):
        segment_points_fine = [232, 356, 547, 694, 790, 898, 935, 1044, 1108, 1237, 1266, 1376, 1445, 1564, 1608, 1730, 1789, 1917, 1943, 2069, 2147, 2294, 2309, 2433]
    elif (activity == "Squat") and (segment_points == [151, 261, 362, 472, 572, 678, 775, 881, 986, 1090, 1200]):
        segment_points_fine = [173, 263, 263, 350, 367, 474, 474, 572, 572, 680, 680, 768, 780, 883, 883, 986, 986, 1092, 1092, 1200]
    elif (activity == "Calf raises") and (segment_points == [44, 214, 384, 534, 671, 813, 960, 1105, 1262, 1420, 1590]):
        segment_points_fine = [152, 214, 214, 329, 441, 534, 534, 625, 711, 807, 826, 954, 1006, 1105, 1105, 1214, 1328, 1412, 1426, 1522]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [25, 260, 409, 560, 702, 847, 1002, 1140, 1291, 1432, 1667]):
        segment_points_fine = [117, 260, 264, 383, 420, 559, 565, 696, 708, 847, 847, 987, 1007, 1140, 1140, 1286, 1298, 1432, 1432, 1560]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [39, 342, 519, 690, 888, 1075, 1263, 1470, 1719, 2022, 2325]):
        segment_points_fine = [103, 252, 349, 478, 542, 683, 694, 845, 904, 1067, 1090, 1217, 1290, 1465, 1470, 1689, 1719, 1891, 2067, 2217]
    elif (activity == "Calf raises") and (segment_points == [115, 200, 282, 351, 436, 511, 582, 664, 730, 804, 889]):
        segment_points_fine = [148, 195, 200, 261, 300, 346, 351, 416, 456, 501, 511, 560, 599, 654, 664, 705, 747, 799, 804, 854]
    elif (activity == "Calf raises") and (segment_points == [0, 187, 324, 447, 595, 812, 945, 1079, 1230, 1391, 1662, 1860, 1975]):
        segment_points_fine = [103, 177, 197, 265, 353, 433, 500, 573, 643, 705, 838, 884, 982, 1042, 1108, 1166, 1280, 1342, 1413, 1473, 1769, 1820, 1870, 1912]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [92, 308, 524, 686, 902, 1062, 1227, 1428, 1612, 1799, 2001, 2176, 2392]):
        segment_points_fine = [146, 295, 325, 462, 532, 676, 719, 861, 906, 1044, 1066, 1227, 1277, 1419, 1445, 1597, 1643, 1776, 1830, 1961, 2001, 2163, 2220, 2337]
    elif (activity == "Calf raises") and (segment_points == [101, 247, 380, 497, 643, 758, 856, 1000, 1114, 1226, 1372]):
        segment_points_fine = [176, 231, 279, 317, 431, 490, 534, 580, 698, 748, 776, 831, 928, 979, 1001, 1055, 1165, 1216, 1237, 1291]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [53, 280, 485, 654, 881, 1025, 1230, 1407, 1599, 1754, 1952, 2148, 2580]):
        segment_points_fine = [134, 265, 297, 444, 524, 642, 682, 828, 890, 1007, 1056, 1187, 1238, 1366, 1433, 1584, 1607, 1745, 1776, 1919, 1984, 2121, 2341, 2517]
    elif (activity == "Calf raises") and (segment_points == [82, 185, 284, 363, 461, 549, 645, 748, 839, 942, 1080]):
        segment_points_fine = [135, 182, 182, 242, 322, 360, 360, 427, 505, 544, 546, 603, 696, 745, 745, 802, 884, 938, 940, 995]
    elif (activity == "Squat") and (segment_points == [132, 255, 366, 470, 577, 670, 793, 903, 1012, 1118, 1241]):
        segment_points_fine = [150, 258, 258, 362, 372, 473, 473, 576, 585, 673, 673, 791, 802, 906, 906, 1012, 1012, 1131, 1131, 1233]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [10, 292, 574, 748, 991, 1194, 1416, 1662, 1919, 2166, 2448]):
        segment_points_fine = [133, 268, 353, 456, 596, 720, 784, 913, 1021, 1175, 1228, 1380, 1480, 1641, 1686, 1842, 1967, 2132, 2219, 2355]
    elif (activity == "Calf raises") and (segment_points == [0, 127, 231, 346, 475, 612, 751, 860, 985, 1108, 1247]):
        segment_points_fine = [76, 125, 127, 175, 281, 346, 346, 434, 532, 609, 616, 674, 798, 860, 860, 930, 1045, 1100, 1114, 1180]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [61, 261, 441, 600, 800, 945, 1137, 1275, 1475]):
        segment_points_fine = [125, 257, 267, 384, 447, 592, 613, 745, 800, 942, 954, 1077, 1137, 1268, 1285, 1414]
    elif (activity == "Squat") and (segment_points == [93, 209, 325, 416, 527, 628, 725, 828, 934, 1036, 1152]):
        segment_points_fine = [112, 210, 210, 318, 330, 417, 417, 519, 535, 629, 629, 725, 725, 829, 829, 931, 938, 1037, 1037, 1145]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [22, 206, 390, 511, 666, 781, 945, 1074, 1195, 1338, 1522]):
        segment_points_fine = [82, 202, 208, 367, 406, 508, 516, 653, 686, 776, 790, 914, 953, 1058, 1078, 1191, 1221, 1316, 1343, 1470]
    elif (activity == "Calf raises") and (segment_points == [77, 201, 317, 396, 505, 602, 720, 826, 950, 1055, 1179]):
        segment_points_fine = [138, 201, 201, 271, 348, 396, 396, 461, 539, 598, 611, 677, 757, 813, 828, 895, 986, 1044, 1060, 1130]
    elif (activity == "Calf raises") and (segment_points == [0, 248, 451, 600, 821, 1052, 1250, 1506, 1767, 2011, 2295, 2536, 2792, 3061, 3299, 3469, 3737, 3937, 4183]):
        segment_points_fine = [140, 235, 281, 380, 499, 592, 622, 742, 898, 997, 1058, 1143, 1338, 1417, 1586, 1688, 1856, 1944, 2061, 2170, 2379, 2463, 2601, 2707, 2895, 3011, 3078, 3167, 3366, 3450, 3543, 3662, 3811, 3906, 3974, 4084]
    elif (activity == "Sit-to-stand and stand-to-sit") and (segment_points == [0, 412, 673, 906, 1146, 1368, 1839, 2193, 2502, 2751, 3167, 3418, 3644, 4096, 4571, 4781, 5188, 5584, 5850, 6325]):
        segment_points = [0, 412, 673, 906, 1146, 1368, 1839, 2193, 2502, 2751, 3167, 3644, 4096, 4571, 4781, 5188, 5584, 5850, 6325]
        segment_points_fine = [233, 384, 460, 638, 707, 845, 961, 1122, 1177, 1339, 1419, 1581, 1909, 2063, 2250, 2449, 2540, 2699, 2867, 3097, 3212, 3390, 3689, 3867, 4192, 4370, 4588, 4748, 4837, 5011, 5363, 5522, 5638, 5819, 5956, 6131]

    # Fix test data
    if (activity == "Calf raises weight") and (segment_points == [62, 220, 328, 486, 644]):
        segment_points_fine = [112, 200, 220, 328, 393, 486, 490, 577]
    elif (activity == "Calf raises obstacle") and (segment_points == [12, 210, 362, 560, 758]):
        segment_points_fine = [82, 210, 210, 362, 362, 560, 560, 688]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [14, 266, 518, 726, 978]):
        segment_points_fine = [80, 207, 281, 422, 567, 707, 833, 958]
        
    elif (activity == "Calf raises weight") and (segment_points == [0, 283, 541]):
        segment_points_fine = [125, 188, 317, 418]
    elif (activity == "Calf raises obstacle") and (segment_points == [0, 268, 463]):
        segment_points_fine = [120, 188, 318, 393]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [0, 554, 1456]):
        segment_points_fine = [373, 482, 610, 779]
    elif (activity == "Squat weight") and (segment_points == [0, 209, 689]):
        segment_points_fine = [99, 209, 319, 447]
    elif (activity == "Squat obstacle") and (segment_points == [0, 162, 572]):
        segment_points_fine = [43, 162, 292, 412]

    elif (activity == "Calf raises obstacle") and (segment_points == [241, 337, 433, 529, 625]):
        segment_points_fine = [261, 337, 337, 433, 433, 529, 529, 595]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [236, 410, 584, 696, 870]):
        segment_points_fine = [290, 400, 424, 544, 584, 685, 712, 827]

    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [131, 333, 533, 735, 868, 1070]):
        segment_points == [333, 533, 735, 868, 1070]
        segment_points_fine = [390, 505, 555, 675, 745, 855, 886, 1004]

    elif (activity == "Calf raises obstacle") and (segment_points == [0, 206, 379]):
        segment_points_fine = [105, 156, 216, 274]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [0, 307, 733]):
        segment_points_fine = [169, 279, 336, 449]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [0, 275, 526]):
        segment_points_fine = [129, 253, 302, 447]
    elif (activity == "Squat weight") and (segment_points == [0, 355, 853]):
        segment_points_fine = [200, 325, 365, 494]

    elif (activity == "Calf raises obstacle") and (segment_points == [0, 238, 486]):
        segment_points_fine = [120, 210, 250, 330]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [0, 340, 553]):
        segment_points_fine = [188, 324, 354, 470]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [0, 313, 600]):
        segment_points_fine = [156, 281, 330, 450]
    elif (activity == "Squat obstacle") and (segment_points == [0, 240, 522]):
        segment_points_fine = [120, 234, 245, 350]

    elif (activity == "Calf raises obstacle") and (segment_points == [197, 367, 506, 676, 846]):
        segment_points_fine = [257, 367, 367, 506, 506, 676, 676, 786]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [336, 563, 681, 860, 1087, 1314]):
        segment_points == [336, 563, 860, 1087, 1314]
        segment_points_fine = [402, 548, 598, 772, 910, 1072, 1119, 1285]

    elif (activity == "Calf raises obstacle") and (segment_points == [103, 216, 326, 439, 552]):
        segment_points_fine = [123, 186, 226, 276, 346, 409, 449, 522]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [83, 247, 411, 564, 728]):
        segment_points_fine = [110, 231, 270, 385, 434, 541, 584, 698]

    elif (activity == "Calf raises obstacle") and (segment_points == [117, 214, 311, 390, 487]):
        segment_points_fine = [134, 197, 221, 275, 339, 388, 400, 458]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [128, 284, 437, 593, 749]):
        segment_points_fine = [128, 268, 296, 416, 444, 581, 607, 720]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [0, 184, 272, 458, 621, 807]):
        segment_points = [0, 272, 458, 621, 807]
        segment_points_fine = [120, 272, 272, 427, 471, 615, 632, 762]
    elif (activity == "Squat obstacle") and (segment_points == [115, 249, 383, 493, 613, 733]):
        segment_points = [115, 249, 383, 493, 613]
        segment_points_fine = [127, 246, 259, 361, 393, 490, 500, 592]

    elif (activity == "Calf raises obstacle") and (segment_points == [90, 176, 257, 343, 429]):
        segment_points_fine = [122, 166, 184, 242, 276, 331, 351, 394]
    elif (activity == "Calf raises weight") and (segment_points == [103, 196, 289, 376, 469]):
        segment_points_fine = [145, 196, 196, 266, 308, 376, 376, 433]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [76, 284, 492, 663, 871]):
        segment_points_fine = [144, 279, 304, 422, 498, 645, 687, 805]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [76, 309, 542, 742, 975]):
        segment_points_fine = [152, 309, 309, 454, 562, 725, 767, 904]
    elif (activity == "Squat obstacle") and (segment_points == [95, 201, 291, 397, 503]):
        segment_points_fine = [113, 201, 205, 286, 320, 391, 402, 474]
    elif (activity == "Squat weight") and (segment_points == [335, 520, 705, 842, 958, 1143]):
        segment_points = [520, 705, 842, 958, 1143]
        segment_points_fine = [608, 705, 705, 810, 860, 958, 958, 1062]

    elif (activity == "Calf raises obstacle") and (segment_points == [105, 212, 297, 404, 511]):
        segment_points_fine = [128, 186, 218, 284, 341, 394, 429, 486]
    elif (activity == "Calf raises weight") and (segment_points == [110, 205, 300, 386, 481]):
        segment_points_fine = [122, 175, 215, 268, 322, 378, 412, 475]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [95, 286, 477, 640, 831]):
        segment_points_fine = [114, 247, 326, 442, 497, 595, 686, 806]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [2, 317, 517, 721, 1006, 1321, 1576]):
        segment_points = [2, 317, 517, 721, 1006]
        segment_points_fine = [152, 305, 351, 487, 534, 695, 743, 900]

    elif (activity == "Calf raises obstacle") and (segment_points == [76, 206, 308, 438, 568]):
        segment_points_fine = [106, 176, 216, 278, 328, 433, 443, 538]
    elif (activity == "Calf raises weight") and (segment_points == [1612, 1728, 1844, 1955, 2071]):
        segment_points_fine = [1605, 1708, 1738, 1832, 1852, 1941, 1972, 2076]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [40, 297, 554, 770, 1027]):
        segment_points_fine = [123, 258, 331, 491, 602, 732, 804, 969]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [66, 290, 470, 694, 918]):
        segment_points_fine = [130, 254, 301, 460, 530, 670, 721, 880]

    elif (activity == "Calf raises obstacle") and (segment_points == [12, 141, 270, 390, 519]):
        segment_points_fine = [52, 136, 151, 230, 305, 385, 395, 469]
    elif (activity == "Calf raises weight") and (segment_points == [24, 190, 330, 496, 650]):
        segment_points_fine = [85, 185, 195, 302, 391, 491, 499, 578]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [119, 333, 547, 748, 962]):
        segment_points_fine = [179, 301, 366, 494, 570, 715, 790, 912]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [0, 211, 458, 685, 932]):
        segment_points_fine = [18, 197, 224, 399, 522, 676, 710, 855]

    elif (activity == "Calf raises obstacle") and (segment_points == [14, 110, 201, 297, 393]):
        segment_points_fine = [44, 110, 110, 181, 221, 297, 297, 373]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [109, 256, 403, 537, 684]):
        segment_points_fine = [143, 242, 267, 352, 427, 525, 553, 654]
    elif (activity == "Sit-to-stand and stand-to-sit weight") and (segment_points == [96, 246, 396, 536, 686]):
        segment_points_fine = [133, 244, 249, 342, 432, 529, 550, 641]

    elif (activity == "Calf raises obstacle") and (segment_points == [110, 224, 335, 449, 563]):
        segment_points_fine = [160, 214, 234, 315, 385, 439, 459, 533]
    elif (activity == "Sit-to-stand and stand-to-sit obstacle") and (segment_points == [77, 222, 367, 507, 652]):
        segment_points_fine = [107, 211, 240, 357, 372, 487, 527, 632]
    elif (activity == "Squat obstacle") and (segment_points == [93, 194, 295, 396]):
        segment_points = [93, 194, 295]
        segment_points_fine = [99, 179, 203, 288]
        

    # If something goes wrong, uncomment this plot to troubleshoot
    print(segment_points)
    print(segment_points_fine)
    plt.figure()
    plt.title(activity)
    plt.plot(avg_bone_speed, label="avg_bone_speed")
    plt.scatter(segment_points, avg_bone_speed[segment_points], label="segment_points")
    plt.scatter(segment_points_fine, avg_bone_speed[segment_points_fine], marker="+", label="segment_points_fine")
    plt.legend()
    plt.show()

    # Segment and update
    for i in np.arange(int(np.round(len(segment_points_fine)/2))):
        # Label correctly
        if (i%2) == 0:
            activity = first_act
        else:
            activity = last_act

        # Segmentation stuff
        start_seg = segment_points_fine[2*i]
        end_seg = segment_points_fine[2*i+1]
        n_pts_data = end_seg - start_seg # not +1 because last point will not be included
        
        # Construct "Label", "GaitPhase", "GaitSpeed" and "NewFile" column
        labels_dummy = np.zeros((n_pts_data, 1))
        labels = [activity for i in np.arange(n_pts_data)]
        
        gait_phase = np.linspace(0, 1, num=n_pts_data, endpoint=False).reshape((-1, 1))

        gait_speed = np.ones((n_pts_data, 1))*(100/n_pts_data)
        
        new_file = np.zeros((n_pts_data, 1))
        new_file[0] = 1
        
        # Construct dataframe and store
        arr_all_data = np.concatenate((labels_dummy, gait_phase, gait_speed, new_file, 
                                       arr_mot[start_seg:end_seg, :], 
                                       arr_imu[start_seg:end_seg, :], 
                                       arr_gon[start_seg:end_seg, :], 
                                       arr_emg[start_seg:end_seg, :]), axis=1)
        df_all_data = pd.DataFrame(arr_all_data, columns=COLS_TRAINING_DATA)
        df_all_data["Label"] = labels
        df_all_data.to_csv(fileTrainingData, header=False, index=False, mode="a")

    return

# Append data to the training data file for cyclic activities
def append_cyclic_data(fileTrainingData, activity, arr_foot_markers, arr_mot, arr_imu, arr_gon, arr_emg):
    # Change name of some activities
    if (activity == "Stair ascent (SOS)(L first)") or (activity == "Stair ascent (SOS)(R first)"):
        activity = "Stair ascent (SOS)"
    elif (activity == "Stair descent (SOS)(L first)") or (activity == "Stair descent (SOS)(R first)"):
        activity = "Stair descent (SOS)"
    elif (activity == "Stair ascent (SOS)(L first) weight") or (activity == "Stair ascent (SOS)(R first) weight"):
        activity = "Stair ascent (SOS) weight"
    elif (activity == "Stair descent (SOS)(L first) weight") or (activity == "Stair descent (SOS)(R first) weight"):
        activity = "Stair descent (SOS) weight"
    elif (activity == "Stair ascent (SOS)(L first) obstacle") or (activity == "Stair ascent (SOS)(R first) obstacle"):
        activity = "Stair ascent (SOS) obstacle"
    elif (activity == "Stair descent (SOS)(L first) obstacle") or (activity == "Stair descent (SOS)(R first) obstacle"):
        activity = "Stair descent (SOS) obstacle"
    
    # Determine usable ranges of data, as segmentation is based on RTOE and RHEE markers
    n_pts = np.shape(arr_foot_markers)[0]
    usable_rows = np.logical_not(np.any(np.isnan(arr_foot_markers), axis=1))
    usable_blocks = []
    in_block = False
    for i, e in enumerate(usable_rows):
        if (not in_block) and e:
            start_block = i
            in_block = True
        elif in_block and (not e):
            end_block = i-1
            n_block = end_block-start_block+1
            if n_block >= 100:
                usable_blocks.append((start_block, end_block))
            in_block = False
        else:
            pass
            
    # Check last point, outside for-loop to speed things up
    if in_block and e:
        end_block = i
        n_block = end_block-start_block+1
        if n_block >= 100:
            usable_blocks.append((start_block, end_block))

    # Treat every block of possibly usable data separately
    for block in usable_blocks:
        start_ind_block, stop_ind_block = block
        n_block = stop_ind_block-start_ind_block+1
        
        # Get footstrikes for each block
        foot_strikes = get_foot_strikes_train(arr_foot_markers[start_ind_block:stop_ind_block+1, 0:3], 
                                              arr_foot_markers[start_ind_block:stop_ind_block+1, 3:6], 
                                              f_sample=100, threshold_speed=0.2, filt_freq=None)

        # Check if at least one full stride is available and calculate accompanying gait phase and speed
        gait_phase = np.empty((n_block, 1))
        gait_phase[:] = np.nan
        gait_speed = np.empty((n_block, 1))
        gait_speed[:] = np.nan
        n_foot_strikes = len(foot_strikes)
        
        # Check
        if n_foot_strikes >= 2:
            # Gait phase
            for i in np.arange(n_foot_strikes-1):
                start_ind_stride = foot_strikes[i]
                stop_ind_stride = foot_strikes[i+1]-1
                n_pts_stride = stop_ind_stride-start_ind_stride+1
                gait_phase[start_ind_stride:stop_ind_stride+1] = np.linspace(0, 1, num=n_pts_stride, endpoint=False).reshape((-1, 1))
                gait_speed[start_ind_stride:stop_ind_stride+1] = np.ones((n_pts_stride, 1))*(100/n_pts_stride)
            gait_phase = gait_phase[foot_strikes[0]:foot_strikes[-1]]
            gait_speed = gait_speed[foot_strikes[0]:foot_strikes[-1]]

            # Construct "Label" and "NewFile" column
            start_ind_data = start_ind_block + foot_strikes[0]
            stop_ind_data = start_ind_block + foot_strikes[-1] - 1
            n_pts_data = stop_ind_data - start_ind_data + 1
            
            new_file = np.zeros((n_pts_data, 1))
            new_file[0] = 1
            
            labels_dummy = np.zeros((n_pts_data, 1))
            labels = [activity for i in np.arange(n_pts_data)]
            
            # Construct dataframe and store
            arr_all_data = np.concatenate((labels_dummy, gait_phase, gait_speed, new_file, 
                                           arr_mot[start_ind_data:stop_ind_data+1], 
                                           arr_imu[start_ind_data:stop_ind_data+1], 
                                           arr_gon[start_ind_data:stop_ind_data+1],
                                           arr_emg[start_ind_data:stop_ind_data+1]), axis=1)
            df_all_data = pd.DataFrame(arr_all_data, columns=COLS_TRAINING_DATA)
            df_all_data["Label"] = labels
            df_all_data.to_csv(fileTrainingData, header=False, index=False, mode="a")
    return

# Append data to the training data file for the combinations of activities
def append_combination_data(fileTrainingData, activity, arr_mot, arr_imu, arr_gon, arr_emg):
    
    # Construct "Label", "GaitPhase", "GaitSpeed" and "NewFile" column
    n_pts_data = np.shape(arr_mot)[0]
    labels_dummy = np.zeros((n_pts_data, 1))
    labels = [activity for i in np.arange(n_pts_data)]
    
    gait_phase = np.empty((n_pts_data, 1))
    gait_phase[:] = np.nan

    gait_speed = np.empty((n_pts_data, 1))
    gait_speed[:] = np.nan
    
    new_file = np.zeros((n_pts_data, 1))
    new_file[0] = 1
    
    # Construct dataframe and store
    arr_all_data = np.concatenate((labels_dummy, gait_phase, gait_speed, new_file, arr_mot, arr_imu, arr_gon, arr_emg), axis=1)
    df_all_data = pd.DataFrame(arr_all_data, columns=COLS_TRAINING_DATA)
    df_all_data["Label"] = labels
    df_all_data.to_csv(fileTrainingData, header=False, index=False, mode="a")
    return


###################
# Data management #
###################

# Get list of all the folders with data to read and the corresponding subject info files
def get_data_locations(dataFolder, subjectSessions, dataToRead):
    foldersToRead = []
    subjectInfoFiles = []
    subjects = subjectSessions.keys()
    
    for subject in subjects:
        subjectFolder = "/".join((dataFolder, subject))
        subjectFolderDirList = os.listdir(subjectFolder)
        sessions = subjectSessions[subject]
        for session in sessions:
            sessionFolder = "/".join((dataFolder, subject, session))
            sessionDataToRead = sessionFolder + dataToRead
            if (session + dataToRead) in subjectFolderDirList:
                foldersToRead.append(sessionDataToRead)
                subjectInfoFiles.append(sessionFolder + "/" + subject + " - info.xlsx")

    return foldersToRead, subjectInfoFiles

# Get list of tuples with every file to read and its corresponding activity, as well as a list with corresponding subject information
def get_filenames_ref(foldersToRead, subjectInfoFiles):
    filesToRead = []
    subjectInfo = []
    
    # Loop over all the measurement sessions
    for i, folder in enumerate(foldersToRead):
        fileList = os.listdir(folder)
        subjectInfoFile = pd.read_excel(subjectInfoFiles[i])
        filesChecked = []
        
        # Loop over all the different activities and their iterations
        for file in fileList:
            splitFileName = file.split(" - ")
            activity = splitFileName[2]
            
            # To not repeat same activity
            try:
                if (activity[-3] == "(") and (activity[-1] == ")") and (int(activity[-2]) in range(10)):
                    activity = activity[:-3]
            except:
                pass
                
            splitFileName[-1] = splitFileName[-1].split(".")[0] # Remove extension
            fileNameBase = " - ".join(splitFileName[0:4])
    
            # Update list of checked sessions
            if (fileNameBase not in filesChecked) and (activity not in IGNORE_REF):
                filesChecked.append(fileNameBase)
        
                # Find the right flags
                flags = ""
                if " - badMarker" in file:
                    flags = flags + " - badMarker"
                if " - badForcePlate" in file:
                    flags = flags + " - badForcePlate"
                if " - badStep" in file:
                    flags = flags + " - badStep"
                
                # Get all five files for one step
                fileGRF = folder + "/" + fileNameBase + flags + ".csv"
                fileFootMarkers = folder + "/" + fileNameBase + " - FootMarkers" + flags + ".csv"
                fileModelOutputs = folder + "/" + fileNameBase + " - ModelOutputs" + flags + ".csv"
                fileEMGFeat = folder + "/" + fileNameBase + " - EMG_features" + flags + ".csv"
                fileGlobalAngles = folder + "/" + fileNameBase + " - GlobalAngles" + flags + ".csv"
    
                # Update
                filesToRead.append((activity, fileGRF, fileFootMarkers, fileModelOutputs, fileEMGFeat, fileGlobalAngles))
                subjectInfo.append(subjectInfoFile)

    return filesToRead, subjectInfo

# Get list of tuples with every file to read and its corresponding activity, as well as a list with corresponding subject information
def get_filenames_train(foldersToRead, subjectInfoFiles, test_data=False):
    filesToRead = []
    subjectInfo = []
    
    # Loop over all the measurement sessions
    for i, folder in enumerate(foldersToRead):
        fileList = os.listdir(folder)
        subjectInfoFile = pd.read_excel(subjectInfoFiles[i])
        filesChecked = []
        
        # Loop over all the different activities and their iterations
        for file in fileList:
            splitFileName = file.split(" - ")
            activity = splitFileName[2]
            
            # To not repeat same activity
            try:
                if (activity[-3] == "(") and (activity[-1] == ")") and (int(activity[-2]) in range(10)):
                    activity = activity[:-3]
            except:
                pass
                
            splitFileName[-1] = splitFileName[-1].split(".")[0] # Remove extension
            fileNameBase = " - ".join(splitFileName[0:4])
    
            # Update list of checked sessions
            if (fileNameBase not in filesChecked) and ((activity not in IGNORE_REF) or test_data):
                filesChecked.append(fileNameBase)
                
                # Get all six files for one step
                fileAllTraj = folder + "/" + fileNameBase + " - AllTraj.csv"
                fileMOTNorm = folder + "/" + fileNameBase + " - MotNorm.csv"
                fileIMUFeat = folder + "/" + fileNameBase + " - IMU_features.csv"
                fileGONFeat = folder + "/" + fileNameBase + " - GON_features.csv"
                fileEMGFeat = folder + "/" + fileNameBase + " - EMG_features.csv"
                fileModOut = folder + "/" + fileNameBase + " - ModOut.csv"
    
                # Update
                filesToRead.append((activity, fileAllTraj, fileMOTNorm, fileIMUFeat, fileGONFeat, fileEMGFeat, fileModOut))
                subjectInfo.append(subjectInfoFile)

    return filesToRead, subjectInfo

# Create activity dictionaries
def initialize_activity_dictionaries(foldersToRead):
    # Get dictionary with every type of activity that was performed as a key
    # The value accompanying each key is a numpy array of size (N_PTS_REPR, N_COLS_REFERENCE_DATA)
    # There are N_PTS_REPR rows representing one full gait cycle for each of the N_COLS_REFERENCE_DATA signals recorded
    # The mean of each of those signals will be calculated in this dictionary
    # This way it contains the reference trajectories for each of these signals for each different activity
    
    activityDict = {} # To calculate sum
    activityDictCounter = {} # To keep track of amount of terms added together
    activitiesChecked = []
    
    # Get activity names
    for folder in foldersToRead:
        fileList = os.listdir(folder)    
        
        for file in fileList:
            splitFileName = file.split(" - ")
            activity = splitFileName[2]
            
            # To not repeat same activity
            try:
                if (activity[-3] == "(") and (activity[-1] == ")") and (int(activity[-2]) in range(10)):
                    activity = activity[:-3]
            except:
                pass
                
            # Update dictionaries
            if (activity not in activitiesChecked) and (activity not in IGNORE_REF):
                if activity == "Calf raises":
                    activityDict["Calf raise up"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Calf raise up"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDict["Calf raise down"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Calf raise down"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif activity == "Sit-to-stand and stand-to-sit":
                    activityDict["Sit-to-stand"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Sit-to-stand"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDict["Stand-to-sit"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Stand-to-sit"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif activity == "Squat":
                    activityDict["Squat up"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Squat up"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDict["Squat down"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Squat down"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif (activity == "Circular walking (diam 2m)(CCW)") or (activity == "Circular walking (diam 2m)(CW)"):
                    activityDict["Circular walking (outer leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Circular walking (outer leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDict["Circular walking (inner leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Circular walking (inner leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif (activity == "Side stepping (L)") or (activity == "Side stepping (R)"):
                    activityDict["Side stepping (leading leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Side stepping (leading leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDict["Side stepping (following leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Side stepping (following leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif (activity == "Stair ascent (SBS)(L first)") or (activity == "Stair ascent (SBS)(R first)"):
                    activityDict["Stair ascent (SBS)(leading leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Stair ascent (SBS)(leading leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDict["Stair ascent (SBS)(following leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Stair ascent (SBS)(following leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif (activity == "Stair descent (SBS)(L first)") or (activity == "Stair descent (SBS)(R first)"):
                    activityDict["Stair descent (SBS)(leading leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Stair descent (SBS)(leading leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDict["Stair descent (SBS)(following leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Stair descent (SBS)(following leg)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif (activity == "Stair ascent (SOS)(L first)") or (activity == "Stair ascent (SOS)(R first)"):
                    activityDict["Stair ascent (SOS)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Stair ascent (SOS)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                elif (activity == "Stair descent (SOS)(L first)") or (activity == "Stair descent (SOS)(R first)"):
                    activityDict["Stair descent (SOS)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter["Stair descent (SOS)"] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                else:
                    activityDict[activity] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                    activityDictCounter[activity] = np.zeros((N_PTS_REPR, N_COLS_REFERENCE_DATA))
                activitiesChecked.append(activity)

    return activityDict, activityDictCounter