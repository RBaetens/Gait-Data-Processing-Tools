#######################
# Visualization stuff #
#######################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import use as use_fig_manager


### Visuals to store

def set_font_size_ticks(axs, size):
    xticks = axs.get_xticks()
    xticks[0] = -5
    xticks[-1] = 105
    xticklabels = axs.get_xticklabels()
    xticklabels[0] = None
    xticklabels[-1] = None
    yticks = axs.get_yticks()
    yticklabels = axs.get_yticklabels()
    axs.set_xticks(xticks, labels=xticklabels, size=size)
    axs.set_yticks(yticks, labels=yticklabels, size=size)
    return

def get_nchar_yticks(axs):
    yticklabels = axs.get_yticklabels()
    nchar = 0
    for label in yticklabels:
        nchar_i = len(label.get_text())
        if nchar_i > nchar:
            nchar = nchar_i
        else:
            pass
    return nchar
        
def classic_gait_plot(processed_data, activity, style="shaded", filtered=False, offset=False, filename=None):
    ### Select right data
    x_axis = np.linspace(0, 99, num=100)
    meanDF = processed_data.meanDict[activity]
    sdsDF = processed_data.sdsDict[activity]

    if offset:
        # Joint angles
        ankleAngle = meanDF["AnkleAnglesX_offset"].to_numpy()
        kneeAngle = meanDF["KneeAnglesX_offset"].to_numpy()
        hipAngle = meanDF["HipAnglesX_offset"].to_numpy()
        ankleAngleSds = sdsDF["AnkleAnglesX_offset"].to_numpy()
        kneeAngleSds = sdsDF["KneeAnglesX_offset"].to_numpy()
        hipAngleSds = sdsDF["HipAnglesX_offset"].to_numpy()
    else:
        # Joint angles
        ankleAngle = meanDF["AnkleAnglesX"].to_numpy()
        kneeAngle = meanDF["KneeAnglesX"].to_numpy()
        hipAngle = meanDF["HipAnglesX"].to_numpy()
        ankleAngleSds = sdsDF["AnkleAnglesX"].to_numpy()
        kneeAngleSds = sdsDF["KneeAnglesX"].to_numpy()
        hipAngleSds = sdsDF["HipAnglesX"].to_numpy()

    if filtered:
        # Joint moments
        ankleMoment = meanDF["AnkleMomentX_filt"].to_numpy()
        kneeMoment = meanDF["KneeMomentX_filt"].to_numpy()
        hipMoment = meanDF["HipMomentX_filt"].to_numpy()
        ankleMomentSds = sdsDF["AnkleMomentX_filt"].to_numpy()
        kneeMomentSds = sdsDF["KneeMomentX_filt"].to_numpy()
        hipMomentSds = sdsDF["HipMomentX_filt"].to_numpy()
    
        # Joint powers
        anklePower = meanDF["AnklePowerZ_filt"].to_numpy()
        kneePower = meanDF["KneePowerZ_filt"].to_numpy()
        hipPower = meanDF["HipPowerZ_filt"].to_numpy()
        anklePowerSds = sdsDF["AnklePowerZ_filt"].to_numpy()
        kneePowerSds = sdsDF["KneePowerZ_filt"].to_numpy()
        hipPowerSds = sdsDF["HipPowerZ_filt"].to_numpy()
    else:
        # Joint moments
        ankleMoment = meanDF["AnkleMomentX"].to_numpy()
        kneeMoment = meanDF["KneeMomentX"].to_numpy()
        hipMoment = meanDF["HipMomentX"].to_numpy()
        ankleMomentSds = sdsDF["AnkleMomentX"].to_numpy()
        kneeMomentSds = sdsDF["KneeMomentX"].to_numpy()
        hipMomentSds = sdsDF["HipMomentX"].to_numpy()
    
        # Joint powers
        anklePower = meanDF["AnklePowerZ"].to_numpy()
        kneePower = meanDF["KneePowerZ"].to_numpy()
        hipPower = meanDF["HipPowerZ"].to_numpy()
        anklePowerSds = sdsDF["AnklePowerZ"].to_numpy()
        kneePowerSds = sdsDF["KneePowerZ"].to_numpy()
        hipPowerSds = sdsDF["HipPowerZ"].to_numpy()

    ### Define style
    suptitle_size = 18
    suptitle_weight = "bold"
    
    title_size = 13
    title_weight = "normal"
    titlepad = 12
    
    ylabel_size = 13
    ylabel_weight = "normal"
    ylabelpad_1 = 45
    ylabelpad_2 = 45
    ylabelpad_3 = 45
    char_comp = -4
    
    xlabel_size = 13
    xlabel_weight = "normal"
    xlabelpad = 12
    
    tick_size = 8
    
    linestyle_mean = "-"
    linewidth_mean = 2
    linestyle_sds = "--"
    linewidth_sds = 2
    
    fig, axs = plt.subplots(3, 3, figsize=(24, 18), dpi=300)
    fig.suptitle(activity, size=suptitle_size, weight=suptitle_weight)

    ### Plot means
    # Joint angles
    axs[0, 0].plot(x_axis, ankleAngle, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[0, 1].plot(x_axis, kneeAngle, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[0, 2].plot(x_axis, hipAngle, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    # Joint moments
    axs[1, 0].plot(ankleMoment, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[1, 1].plot(kneeMoment, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[1, 2].plot(hipMoment, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    # Joint powers
    axs[2, 0].plot(anklePower, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[2, 1].plot(kneePower, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[2, 2].plot(hipPower, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    ### Add standard deviations
    if style == "shaded":
        # Angles
        axs[0, 0].fill_between(x_axis, ankleAngle + ankleAngleSds, ankleAngle - ankleAngleSds, alpha=0.2)
        axs[0, 1].fill_between(x_axis, kneeAngle + kneeAngleSds, kneeAngle - kneeAngleSds, alpha=0.2)
        axs[0, 2].fill_between(x_axis, hipAngle + hipAngleSds, hipAngle - hipAngleSds, alpha=0.2)

        # Moments
        axs[1, 0].fill_between(x_axis, ankleMoment + ankleMomentSds, ankleMoment - ankleMomentSds, alpha=0.2)
        axs[1, 1].fill_between(x_axis, kneeMoment + kneeMomentSds, kneeMoment - kneeMomentSds, alpha=0.2)
        axs[1, 2].fill_between(x_axis, hipMoment + hipMomentSds, hipMoment - hipMomentSds, alpha=0.2)

        # Powers
        axs[2, 0].fill_between(x_axis, anklePower + anklePowerSds, anklePower - anklePowerSds, alpha=0.2)
        axs[2, 1].fill_between(x_axis, kneePower + kneePowerSds, kneePower - kneePowerSds, alpha=0.2)
        axs[2, 2].fill_between(x_axis, hipPower + hipPowerSds, hipPower - hipPowerSds, alpha=0.2)
        
    elif style == "dotted":
        # Angles
        axs[0, 0].plot(x_axis, ankleAngle + ankleAngleSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 0].plot(x_axis, ankleAngle - ankleAngleSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 1].plot(x_axis, kneeAngle + kneeAngleSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 1].plot(x_axis, kneeAngle - kneeAngleSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 2].plot(x_axis, hipAngle + hipAngleSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 2].plot(x_axis, hipAngle - hipAngleSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)

        # Moments
        axs[1, 0].plot(x_axis, ankleMoment + ankleMomentSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 0].plot(x_axis, ankleMoment - ankleMomentSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 1].plot(x_axis, kneeMoment + kneeMomentSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 1].plot(x_axis, kneeMoment - kneeMomentSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 2].plot(x_axis, hipMoment + hipMomentSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 2].plot(x_axis, hipMoment - hipMomentSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)

        # Powers
        axs[2, 0].plot(x_axis, anklePower + anklePowerSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 0].plot(x_axis, anklePower - anklePowerSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 1].plot(x_axis, kneePower + kneePowerSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 1].plot(x_axis, kneePower - kneePowerSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 2].plot(x_axis, hipPower + hipPowerSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 2].plot(x_axis, hipPower - hipPowerSds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        
    else:
        raise ValueError("'style' should be defined as 'shaded' or 'dotted'.")

    ### Set labels and grids
    # Angles
    axs[0, 0].set_title("Ankle", size=title_size, weight=title_weight, pad=titlepad)
    axs[0, 0].set_ylabel("FE angle [°]", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_1+char_comp*get_nchar_yticks(axs[0, 0])))
    set_font_size_ticks(axs[0, 0], tick_size)
    axs[0, 0].grid()
    
    axs[0, 1].set_title("Knee", size=title_size, weight=title_weight, pad=titlepad)
    set_font_size_ticks(axs[0, 1], tick_size)
    axs[0, 1].grid()

    axs[0, 2].set_title("Hip", size=title_size, weight=title_weight, pad=titlepad)
    set_font_size_ticks(axs[0, 2], tick_size)
    axs[0, 2].grid()

    # Moments
    axs[1, 0].set_ylabel("FE moment [Nmm/kg]", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_2+char_comp*get_nchar_yticks(axs[1, 0])))
    set_font_size_ticks(axs[1, 0], tick_size)
    axs[1, 0].grid()

    set_font_size_ticks(axs[1, 1], tick_size)
    axs[1, 1].grid()

    set_font_size_ticks(axs[1, 2], tick_size)
    axs[1, 2].grid()

    # Powers
    axs[2, 0].set_xlabel("Gait phase [%]", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    axs[2, 0].set_ylabel("FE power [W/kg]", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_3+char_comp*get_nchar_yticks(axs[2, 0])))
    set_font_size_ticks(axs[2, 0], tick_size)
    axs[2, 0].grid()

    axs[2, 1].set_xlabel("Gait phase [%]", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    set_font_size_ticks(axs[2, 1], tick_size)
    axs[2, 1].grid()

    axs[2, 2].set_xlabel("Gait phase [%]", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    set_font_size_ticks(axs[2, 2], tick_size)
    axs[2, 2].grid()

    ### Done!
    if filename:
        fig.savefig(filename)
        plt.close(fig)
    else:
        plt.show()
    
    return

def emg_plot(processed_data, activity, style="shaded", signal="ACT", filename=None):
    ### Select right data
    x_axis = np.linspace(0, 99, num=100)
    meanDF = processed_data.meanDict[activity]
    sdsDF = processed_data.sdsDict[activity]

    if signal == "ACT":
        # Activation signal means
        act_ta = meanDF["ACT TA"].to_numpy()
        act_gm = meanDF["ACT GM"].to_numpy()
        act_gl = meanDF["ACT GL"].to_numpy()
        act_vm = meanDF["ACT VM"].to_numpy()
        act_vl = meanDF["ACT VL"].to_numpy()
        act_rf = meanDF["ACT RF"].to_numpy()
        act_st = meanDF["ACT ST"].to_numpy()
        act_bf = meanDF["ACT BF"].to_numpy()
    
        # Activation signals standard deviations
        act_ta_sds = sdsDF["ACT TA"].to_numpy()
        act_gm_sds = sdsDF["ACT GM"].to_numpy()
        act_gl_sds = sdsDF["ACT GL"].to_numpy()
        act_vm_sds = sdsDF["ACT VM"].to_numpy()
        act_vl_sds = sdsDF["ACT VL"].to_numpy()
        act_rf_sds = sdsDF["ACT RF"].to_numpy()
        act_st_sds = sdsDF["ACT ST"].to_numpy()
        act_bf_sds = sdsDF["ACT BF"].to_numpy()

    elif signal == "RMS":
        # Root mean square signal means
        act_ta = meanDF["RMS TA"].to_numpy()
        act_gm = meanDF["RMS GM"].to_numpy()
        act_gl = meanDF["RMS GL"].to_numpy()
        act_vm = meanDF["RMS VM"].to_numpy()
        act_vl = meanDF["RMS VL"].to_numpy()
        act_rf = meanDF["RMS RF"].to_numpy()
        act_st = meanDF["RMS ST"].to_numpy()
        act_bf = meanDF["RMS BF"].to_numpy()
    
        # Root mean square signal standard deviations
        act_ta_sds = sdsDF["RMS TA"].to_numpy()
        act_gm_sds = sdsDF["RMS GM"].to_numpy()
        act_gl_sds = sdsDF["RMS GL"].to_numpy()
        act_vm_sds = sdsDF["RMS VM"].to_numpy()
        act_vl_sds = sdsDF["RMS VL"].to_numpy()
        act_rf_sds = sdsDF["RMS RF"].to_numpy()
        act_st_sds = sdsDF["RMS ST"].to_numpy()
        act_bf_sds = sdsDF["RMS BF"].to_numpy()

    else:
        raise ValueError("'signal' should be 'ACT' or 'RMS'")

    ### Define style
    suptitle_size = 18
    suptitle_weight = "bold"
    
    title_size = 13
    title_weight = "normal"
    titlepad = 10
    
    ylabel_size = 13
    ylabel_weight = "normal"
    ylabelpad_1 = 25
    ylabelpad_2 = 25
    ylabelpad_3 = 25
    char_comp = -4
    
    xlabel_size = 13
    xlabel_weight = "normal"
    xlabelpad = 12
    
    tick_size = 8
    
    linestyle_mean = "-"
    linewidth_mean = 2
    linestyle_sds = "--"
    linewidth_sds = 2
    
    fig, axs = plt.subplots(4, 2, figsize=(24, 18), dpi=300)
    fig.suptitle(activity, size=suptitle_size, weight=suptitle_weight)
    plt.subplots_adjust(hspace=0.7)

    ### Plot means
    axs[0, 0].plot(x_axis, act_ta, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[0, 1].plot(x_axis, act_gm, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    axs[1, 0].plot(x_axis, act_gl, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[1, 1].plot(x_axis, act_vm, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    axs[2, 0].plot(x_axis, act_vl, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[2, 1].plot(x_axis, act_rf, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    axs[3, 0].plot(x_axis, act_st, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[3, 1].plot(x_axis, act_bf, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    ### Add standard deviations
    if style == "shaded":
        axs[0, 0].fill_between(x_axis, act_ta + 1.5*act_ta_sds, act_ta - 0.5*act_ta_sds, alpha=0.2)
        axs[0, 1].fill_between(x_axis, act_gm + 1.5*act_gm_sds, act_gm - 0.5*act_gm_sds, alpha=0.2)

        axs[1, 0].fill_between(x_axis, act_gl + 1.5*act_gl_sds, act_gl - 0.5*act_gl_sds, alpha=0.2)
        axs[1, 1].fill_between(x_axis, act_vm + 1.5*act_vm_sds, act_vm - 0.5*act_vm_sds, alpha=0.2)

        axs[2, 0].fill_between(x_axis, act_vl + 1.5*act_vl_sds, act_vl - 0.5*act_vl_sds, alpha=0.2)
        axs[2, 1].fill_between(x_axis, act_rf + 1.5*act_rf_sds, act_rf - 0.5*act_rf_sds, alpha=0.2)

        axs[3, 0].fill_between(x_axis, act_st + 1.5*act_st_sds, act_st - 0.5*act_st_sds, alpha=0.2)
        axs[3, 1].fill_between(x_axis, act_bf + 1.5*act_bf_sds, act_bf - 0.5*act_bf_sds, alpha=0.2)
        
    elif style == "dotted":
        axs[0, 0].plot(x_axis, act_ta + 1.5*act_ta_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 0].plot(x_axis, act_ta - 0.5*act_ta_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 1].plot(x_axis, act_gm + 1.5*act_gm_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 1].plot(x_axis, act_gm - 0.5*act_gm_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)

        axs[1, 0].plot(x_axis, act_gl + 1.5*act_gl_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 0].plot(x_axis, act_gl - 0.5*act_gl_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 1].plot(x_axis, act_vm + 1.5*act_vm_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 1].plot(x_axis, act_vm - 0.5*act_vm_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)

        axs[2, 0].plot(x_axis, act_vl + 1.5*act_vl_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 0].plot(x_axis, act_vl - 0.5*act_vl_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 1].plot(x_axis, act_rf + 1.5*act_rf_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[2, 1].plot(x_axis, act_rf - 0.5*act_rf_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)

        axs[3, 0].plot(x_axis, act_st + 1.5*act_st_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[3, 0].plot(x_axis, act_st - 0.5*act_st_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[3, 1].plot(x_axis, act_bf + 1.5*act_bf_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[3, 1].plot(x_axis, act_bf - 0.5*act_bf_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
    
    elif style == "line":
        pass
        
    else:
        raise ValueError("'style' should be defined as 'shaded', 'line' or 'dotted'.")

    ### Set labels and grids
    axs[0, 0].set_title("Tibialis anterior", size=title_size, weight=title_weight, pad=titlepad)
    axs[0, 0].set_ylabel("Activation", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_1+char_comp*get_nchar_yticks(axs[0, 0])))
    set_font_size_ticks(axs[0, 0], tick_size)
    axs[0, 0].grid()
    
    axs[0, 1].set_title("Gastrocnemius medialis", size=title_size, weight=title_weight, pad=titlepad)
    set_font_size_ticks(axs[0, 1], tick_size)
    axs[0, 1].grid()

    axs[1, 0].set_title("Gastrocnemius lateralis", size=title_size, weight=title_weight, pad=titlepad)
    axs[1, 0].set_ylabel("Activation", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_2+char_comp*get_nchar_yticks(axs[1, 0])))
    set_font_size_ticks(axs[1, 0], tick_size)
    axs[1, 0].grid()

    axs[1, 1].set_title("Vastus medialis", size=title_size, weight=title_weight, pad=titlepad)
    set_font_size_ticks(axs[1, 1], tick_size)
    axs[1, 1].grid()

    axs[2, 0].set_title("Vastus lateralis", size=title_size, weight=title_weight, pad=titlepad)
    axs[2, 0].set_ylabel("Activation", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_3+char_comp*get_nchar_yticks(axs[2, 0])))
    set_font_size_ticks(axs[2, 0], tick_size)
    axs[2, 0].grid()

    axs[2, 1].set_title("Rectus femoris", size=title_size, weight=title_weight, pad=titlepad)
    set_font_size_ticks(axs[2, 1], tick_size)
    axs[2, 1].grid()

    axs[3, 0].set_title("Semitendinosus", size=title_size, weight=title_weight, pad=titlepad)
    axs[3, 0].set_xlabel("Gait phase [%]", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    axs[3, 0].set_ylabel("Activation", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_3+char_comp*get_nchar_yticks(axs[3, 0])))
    set_font_size_ticks(axs[3, 0], tick_size)
    axs[3, 0].grid()

    axs[3, 1].set_title("Biceps femoris", size=title_size, weight=title_weight, pad=titlepad)
    axs[3, 1].set_xlabel("Gait phase [%]", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    set_font_size_ticks(axs[3, 1], tick_size)
    axs[3, 1].grid()

    ### Done!
    if filename:
        fig.savefig(filename)
        plt.close(fig)
    else:
        plt.show()
    
    return

def grf_plot(processed_data, activity, style="shaded", filename=None):
    ### Select right data
    x_axis = np.linspace(0, 99, num=100)
    meanDF = processed_data.meanDict[activity]
    sdsDF = processed_data.sdsDict[activity]

    # Ground reaction force means
    grfx = meanDF["GroundReactionForceX"].to_numpy()
    grfy = meanDF["GroundReactionForceY"].to_numpy()
    grfz = meanDF["GroundReactionForceZ"].to_numpy()
    grmz = meanDF["GroundReactionMomentZ"].to_numpy()

    # Ground reaction force standard deviations
    grfx_sds = sdsDF["GroundReactionForceX"].to_numpy()
    grfy_sds = sdsDF["GroundReactionForceY"].to_numpy()
    grfz_sds = sdsDF["GroundReactionForceZ"].to_numpy()
    grmz_sds = sdsDF["GroundReactionMomentZ"].to_numpy()
    
    ### Define style
    suptitle_size = 18
    suptitle_weight = "bold"
    
    title_size = 13
    title_weight = "normal"
    titlepad = 10
    
    ylabel_size = 13
    ylabel_weight = "normal"
    ylabelpad_1 = 25
    ylabelpad_2 = 25
    ylabelpad_3 = 25
    char_comp = -4
    
    xlabel_size = 13
    xlabel_weight = "normal"
    xlabelpad = 12
    
    tick_size = 8
    
    linestyle_mean = "-"
    linewidth_mean = 2
    linestyle_sds = "--"
    linewidth_sds = 2
    
    fig, axs = plt.subplots(2, 2, figsize=(24, 18), dpi=300)
    fig.suptitle(activity, size=suptitle_size, weight=suptitle_weight)
    plt.subplots_adjust(hspace=0.3)

    ### Plot means
    axs[0, 0].plot(x_axis, grfx, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[0, 1].plot(x_axis, grfy, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    axs[1, 0].plot(x_axis, grfz, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)
    axs[1, 1].plot(x_axis, grmz, color=(0, 0, 0), linestyle=linestyle_mean, linewidth=linewidth_mean)

    ### Add standard deviations
    if style == "shaded":
        axs[0, 0].fill_between(x_axis, grfx + grfx_sds, grfx - grfx_sds, alpha=0.2)
        axs[0, 1].fill_between(x_axis, grfy + grfy_sds, grfy - grfy_sds, alpha=0.2)

        axs[1, 0].fill_between(x_axis, grfz + grfz_sds, grfz - grfz_sds, alpha=0.2)
        axs[1, 1].fill_between(x_axis, grmz + grmz_sds, grmz - grmz_sds, alpha=0.2)
        
    elif style == "dotted":
        axs[0, 0].plot(x_axis, grfx + grfx_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 0].plot(x_axis, grfx - grfx_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 1].plot(x_axis, grfy + grfy_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[0, 1].plot(x_axis, grfy - grfy_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)

        axs[1, 0].plot(x_axis, grfz + grfz_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 0].plot(x_axis, grfz - grfz_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 1].plot(x_axis, grmz + grmz_sds, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        axs[1, 1].plot(x_axis, grmz - grmz_sdss, color=(1, 0, 0), linestyle=linestyle_sds, linewidth=linewidth_sds)
        
    else:
        raise ValueError("'style' should be defined as 'shaded' or 'dotted'.")

    ### Set labels and grids
    axs[0, 0].set_title("Anteroposterior GRF", size=title_size, weight=title_weight, pad=titlepad)
    axs[0, 0].set_ylabel("Force [N/kg]", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_1+char_comp*get_nchar_yticks(axs[0, 0])))
    set_font_size_ticks(axs[0, 0], tick_size)
    axs[0, 0].grid()
    
    axs[0, 1].set_title("Mediolateral GRF", size=title_size, weight=title_weight, pad=titlepad)
    set_font_size_ticks(axs[0, 1], tick_size)
    axs[0, 1].grid()

    axs[1, 0].set_title("Vertical GRF", size=title_size, weight=title_weight, pad=titlepad)
    axs[1, 0].set_xlabel("Gait phase [%]", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    axs[1, 0].set_ylabel("Force [N/kg] - Moment [Nm/kg]", size=ylabel_size, weight=ylabel_weight, 
                         labelpad=(ylabelpad_2+char_comp*get_nchar_yticks(axs[1, 0])))
    set_font_size_ticks(axs[1, 0], tick_size)
    axs[1, 0].grid()

    axs[1, 1].set_title("Frictional GRM", size=title_size, weight=title_weight, pad=titlepad)
    axs[1, 1].set_xlabel("Gait phase [%]", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    set_font_size_ticks(axs[1, 1], tick_size)
    axs[1, 1].grid()

    ### Done!
    if filename:
        fig.savefig(filename)
        plt.close(fig)
    else:
        plt.show()
    
    return

def cop_traj_plot(processed_data, activity, filename=None):
    ### Select right data
    meanDF = processed_data.meanDict[activity]
    sdsDF = processed_data.sdsDict[activity]

    # Center of pressure coordinate means
    cop_x = meanDF["COPX"].to_numpy()
    cop_y = meanDF["COPY"].to_numpy()

    # Center of pressure coordinate standard deviations
    cop_x_sds = sdsDF["COPX"].to_numpy()
    cop_y_sds = sdsDF["COPY"].to_numpy()

    ### Create plot
    # Create ellipses
    ells = [Ellipse(xy=[cop_x[i], cop_y[i]], width=2*cop_x_sds[i], height=2*cop_y_sds[i], alpha=0.2, color=(i/99, 0, (99-i)/99)) for i in np.arange(100)]

    # Create canvas
    fig, ax = plt.subplots(figsize=(24, 18), dpi=300)
    ax.set(xlim=(np.min(cop_x)-np.max(cop_x_sds), np.max(cop_x)+np.max(cop_x_sds)), 
           ylim=(np.min(cop_y)-np.max(cop_y_sds), np.max(cop_y)+np.max(cop_y_sds)), aspect=1)
    ax.set_axisbelow(True)
    ax.grid()

    # Add ellipses
    for e in ells:
        ax.add_artist(e)

    # Set titles and labels
    title_size = 18
    title_weight = "bold"
    titlepad = 20

    ylabel_size = 13
    ylabel_weight = "normal"
    ylabelpad = 25
    char_comp = -4
    
    xlabel_size = 13
    xlabel_weight = "normal"
    xlabelpad = 12
    
    ax.set_title(activity, size=title_size, weight=title_weight, pad=titlepad)
    ax.set_xlabel("Normalized x-coordinate CoP", size=xlabel_size, weight=xlabel_weight, labelpad=xlabelpad)
    ax.set_ylabel("Normalized y-coordinate CoP", size=ylabel_size, weight=ylabel_weight, labelpad=ylabelpad)
    
    ### Done!
    if filename:
        fig.savefig(filename)
        plt.close(fig)
    else:
        plt.show()
        
    return
    

### Quick visuals in notebook

def visualize_features(features, title=None, labels=None, emg=False):
    """Intended to visualize a set of features belonging to the same activity.
    
    Arguments
    features: [(mean_1: np.array, sds_1: np.array), (mean_2: np.array, sds_2: np.array), ...]
    title: str
    labels: [feature_title1: str, feature_title2: str, ...]
    """

    n_feat = len(features)
    n_rows = int(np.ceil(n_feat/2))
    if n_rows > 5:
        raise ValueError("Plotting more than 10 features at once is not supported.")

    fig, axs = plt.subplots(n_rows, 2, squeeze=False)
    if title:
        fig.suptitle(title)
    fig.set_figwidth(10)
    fig.set_figheight(7)
    fig.tight_layout(pad=2)

    i_feat = 0
    while i_feat < n_feat:
        i_row = i_feat // 2
        i_col = i_feat % 2
        
        axs[i_row, i_col].plot(features[i_feat][0], color=(0, 0, 0), linestyle="-")
        if emg:
            axs[i_row, i_col].plot(features[i_feat][0] + 1.5*features[i_feat][1], color=(1, 0, 0), linestyle="--")
            axs[i_row, i_col].plot(features[i_feat][0] - 0.5*features[i_feat][1], color=(1, 0, 0), linestyle="--")
        else:
            axs[i_row, i_col].plot(features[i_feat][0] + features[i_feat][1], color=(1, 0, 0), linestyle="--")
            axs[i_row, i_col].plot(features[i_feat][0] - features[i_feat][1], color=(1, 0, 0), linestyle="--")
        axs[i_row, i_col].grid()

        if labels:
            axs[i_row, i_col].set_title(labels[i_feat])
            
        i_feat += 1

    if (i_row*2 + i_col + 1) < n_rows*2:
        fig.delaxes(axs[i_row, 1])
    plt.show()
    return

def visualize_cop_traj(COPX, COPY, COPXSds, COPYSds):
    ells = [Ellipse(xy=[COPX[i], COPY[i]], width=2*COPXSds[i], height=2*COPYSds[i], alpha=0.5, color=(i/99, 0, (99-i)/99)) for i in np.arange(100)]
    
    fig, ax = plt.subplots()
    ax.set(xlim=(np.min(COPX)-np.max(COPXSds), np.max(COPX)+np.max(COPXSds)), 
           ylim=(np.min(COPY)-np.max(COPYSds), np.max(COPY)+np.max(COPYSds)), aspect=1)
    ax.set_axisbelow(True)
    ax.grid()
    for e in ells:
        ax.add_artist(e)
    ax.set_title("CoP trajectories")
    ax.set_xlabel("X-coordinate")
    ax.set_ylabel("Y-coordinate")
    plt.show()
    return