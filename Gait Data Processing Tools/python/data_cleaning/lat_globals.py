###################################
### Global variable definitions ###
###################################

# Number of points to represent one gait cycle
N_PTS_REPR = 100


###################################
### Reference data header stuff ###
###################################

# The names of the signals to be stored and some derived variables for indexing
COLS_REFERENCE_DATA = ["AnkleAnglesX", "AnkleAnglesY", "AnkleAnglesZ", "AnkleMomentX", "AnkleMomentY", "AnkleMomentZ", 
                       "AnklePowerX", "AnklePowerY", "AnklePowerZ", "KneeAnglesX", "KneeAnglesY", "KneeAnglesZ", 
                       "KneeMomentX", "KneeMomentY", "KneeMomentZ", "KneePowerX", "KneePowerY", "KneePowerZ", 
                       "HipAnglesX", "HipAnglesY", "HipAnglesZ", "HipMomentX", "HipMomentY", "HipMomentZ", 
                       "HipPowerX", "HipPowerY", "HipPowerZ", 
                       "GroundReactionForceX", "GroundReactionForceY", "GroundReactionForceZ", "GroundReactionMomentZ", 
                       "AnkleAnglesX_offset", "KneeAnglesX_offset", "HipAnglesX_offset", 
                       "AnkleMomentX_filt", "AnkleMomentY_filt", "AnkleMomentZ_filt", "KneeMomentX_filt", "KneeMomentY_filt", "KneeMomentZ_filt", 
                       "HipMomentX_filt", "HipMomentY_filt", "HipMomentZ_filt", "AnklePowerX_filt", "AnklePowerY_filt", "AnklePowerZ_filt", 
                       "KneePowerX_filt", "KneePowerY_filt", "KneePowerZ_filt", "HipPowerX_filt", "HipPowerY_filt", "HipPowerZ_filt", 
                       "RMS TA", "RMS GM", "RMS GL", "RMS VM", "RMS VL", "RMS RF", "RMS ST", "RMS BF", 
                       "MNF TA", "MNF GM", "MNF GL", "MNF VM", "MNF VL", "MNF RF", "MNF ST", "MNF BF", 
                       "ACT TA", "ACT GM", "ACT GL", "ACT VM", "ACT VL", "ACT RF", "ACT ST", "ACT BF", 
                       "COPX", "COPY", "PelvicTilt"]
# COLS_REFERENCE_DATA = ["AnkleAnglesX", "AnkleAnglesY", "AnkleAnglesZ", "AnkleMomentX", "AnkleMomentY", "AnkleMomentZ", 
#                        "AnklePowerX", "AnklePowerY", "AnklePowerZ", "KneeAnglesX", "KneeAnglesY", "KneeAnglesZ", 
#                        "KneeMomentX", "KneeMomentY", "KneeMomentZ", "KneePowerX", "KneePowerY", "KneePowerZ", 
#                        "HipAnglesX", "HipAnglesY", "HipAnglesZ", "HipMomentX", "HipMomentY", "HipMomentZ", 
#                        "HipPowerX", "HipPowerY", "HipPowerZ", 
#                        "GroundReactionForceX", "GroundReactionForceY", "GroundReactionForceZ", "GroundReactionMomentZ", 
#                        "RMS TA", "RMS GM", "RMS GL", "RMS VM", "RMS VL", "RMS RF", "RMS ST", "RMS BF", 
#                        "MNF TA", "MNF GM", "MNF GL", "MNF VM", "MNF VL", "MNF RF", "MNF ST", "MNF BF", 
#                        "ACT TA", "ACT GM", "ACT GL", "ACT VM", "ACT VL", "ACT RF", "ACT ST", "ACT BF", 
#                        "COPX", "COPY", "PelvicTilt"]

FORCES_SLICE = [COLS_REFERENCE_DATA.index("GroundReactionForceX"), COLS_REFERENCE_DATA.index("GroundReactionForceY")]
FRICTION_IND = COLS_REFERENCE_DATA.index("GroundReactionMomentZ")
ALL_FORCES_SLICE = [COLS_REFERENCE_DATA.index("GroundReactionForceX"), COLS_REFERENCE_DATA.index("GroundReactionForceY"), 
                    COLS_REFERENCE_DATA.index("GroundReactionForceZ"), COLS_REFERENCE_DATA.index("GroundReactionMomentZ")]
COP_SLICE = [COLS_REFERENCE_DATA.index("COPX"), COLS_REFERENCE_DATA.index("COPY")]
N_COLS_REFERENCE_DATA = len(COLS_REFERENCE_DATA)
FP_RELATED_HEADERS = ["AnkleMomentX", "AnkleMomentY", "AnkleMomentZ", "AnklePowerX", "AnklePowerY", "AnklePowerZ", 
                      "KneeMomentX", "KneeMomentY", "KneeMomentZ", "KneePowerX", "KneePowerY", "KneePowerZ", 
                      "HipMomentX", "HipMomentY", "HipMomentZ", "HipPowerX", "HipPowerY", "HipPowerZ", 
                      "GroundReactionForceX", "GroundReactionForceY", "GroundReactionForceZ", "GroundReactionMomentZ", 
                      "COPX", "COPY"]
IDX_FP_RELATED_HEADERS = [COLS_REFERENCE_DATA.index(el) for el in FP_RELATED_HEADERS]

# The names of the model output data columns for the left leg
COLS_REFERENCE_DATA_LEFT = ["LAnkleAnglesX", "LAnkleAnglesY", "LAnkleAnglesZ", "LAnkleMomentX", "LAnkleMomentY", "LAnkleMomentZ", 
                            "LAnklePowerX", "LAnklePowerY", "LAnklePowerZ", "LKneeAnglesX", "LKneeAnglesY", "LKneeAnglesZ", 
                            "LKneeMomentX", "LKneeMomentY", "LKneeMomentZ", "LKneePowerX", "LKneePowerY", "LKneePowerZ", 
                            "LHipAnglesX", "LHipAnglesY", "LHipAnglesZ", "LHipMomentX", "LHipMomentY", "LHipMomentZ", 
                            "LHipPowerX", "LHipPowerY", "LHipPowerZ", 
                            "LGroundReactionForceX", "LGroundReactionForceY", "LGroundReactionForceZ", "LGroundReactionMomentZ", 
                            "LAnkleAnglesX_offset", "LKneeAnglesX_offset", "LHipAnglesX_offset", 
                            "LAnkleMomentX_filt", "LAnkleMomentY_filt", "LAnkleMomentZ_filt", "LKneeMomentX_filt", "LKneeMomentY_filt", "LKneeMomentZ_filt", 
                            "LHipMomentX_filt", "LHipMomentY_filt", "LHipMomentZ_filt", "LAnklePowerX_filt", "LAnklePowerY_filt", "LAnklePowerZ_filt", 
                            "LKneePowerX_filt", "LKneePowerY_filt", "LKneePowerZ_filt", "LHipPowerX_filt", "LHipPowerY_filt", "LHipPowerZ_filt"]

# The names of the model output data columns for the right leg
COLS_REFERENCE_DATA_RIGHT = ["RAnkleAnglesX", "RAnkleAnglesY", "RAnkleAnglesZ", "RAnkleMomentX", "RAnkleMomentY", "RAnkleMomentZ", 
                             "RAnklePowerX", "RAnklePowerY", "RAnklePowerZ", "RKneeAnglesX", "RKneeAnglesY", "RKneeAnglesZ", 
                             "RKneeMomentX", "RKneeMomentY", "RKneeMomentZ", "RKneePowerX", "RKneePowerY", "RKneePowerZ", 
                             "RHipAnglesX", "RHipAnglesY", "RHipAnglesZ", "RHipMomentX", "RHipMomentY", "RHipMomentZ", 
                             "RHipPowerX", "RHipPowerY", "RHipPowerZ", 
                             "RGroundReactionForceX", "RGroundReactionForceY", "RGroundReactionForceZ", "RGroundReactionMomentZ", 
                             "RAnkleAnglesX_offset", "RKneeAnglesX_offset", "RHipAnglesX_offset", 
                             "RAnkleMomentX_filt", "RAnkleMomentY_filt", "RAnkleMomentZ_filt", "RKneeMomentX_filt", "RKneeMomentY_filt", "RKneeMomentZ_filt", 
                             "RHipMomentX_filt", "RHipMomentY_filt", "RHipMomentZ_filt", "RAnklePowerX_filt", "RAnklePowerY_filt", "RAnklePowerZ_filt", 
                             "RKneePowerX_filt", "RKneePowerY_filt", "RKneePowerZ_filt", "RHipPowerX_filt", "RHipPowerY_filt", "RHipPowerZ_filt"]

# # The names of the model output data columns for the left leg
# COLS_REFERENCE_DATA_LEFT = ["LAnkleAnglesX", "LAnkleAnglesY", "LAnkleAnglesZ", "LAnkleMomentX", "LAnkleMomentY", "LAnkleMomentZ", 
#                             "LAnklePowerX", "LAnklePowerY", "LAnklePowerZ", "LKneeAnglesX", "LKneeAnglesY", "LKneeAnglesZ", 
#                             "LKneeMomentX", "LKneeMomentY", "LKneeMomentZ", "LKneePowerX", "LKneePowerY", "LKneePowerZ", 
#                             "LHipAnglesX", "LHipAnglesY", "LHipAnglesZ", "LHipMomentX", "LHipMomentY", "LHipMomentZ", 
#                             "LHipPowerX", "LHipPowerY", "LHipPowerZ", 
#                             "LGroundReactionForceX", "LGroundReactionForceY", "LGroundReactionForceZ", "LGroundReactionMomentZ"]

# # The names of the model output data columns for the right leg
# COLS_REFERENCE_DATA_RIGHT = ["RAnkleAnglesX", "RAnkleAnglesY", "RAnkleAnglesZ", "RAnkleMomentX", "RAnkleMomentY", "RAnkleMomentZ", 
#                              "RAnklePowerX", "RAnklePowerY", "RAnklePowerZ", "RKneeAnglesX", "RKneeAnglesY", "RKneeAnglesZ", 
#                              "RKneeMomentX", "RKneeMomentY", "RKneeMomentZ", "RKneePowerX", "RKneePowerY", "RKneePowerZ", 
#                              "RHipAnglesX", "RHipAnglesY", "RHipAnglesZ", "RHipMomentX", "RHipMomentY", "RHipMomentZ", 
#                              "RHipPowerX", "RHipPowerY", "RHipPowerZ", 
#                              "RGroundReactionForceX", "RGroundReactionForceY", "RGroundReactionForceZ", "RGroundReactionMomentZ"]

# Header for the raw EMG files + short version + activation headers
EMG_HEADER = ["tibialis anterior", "gastrocnemius medialis", "gastrocnemius lateralis", "vastus medialis",  
              "vastus lateralis", "rectus femoris", "semitendinosus", "biceps femoris"]
EMG_HEADER_SHORT = ["TA", "GM", "GL", "VM", "VL", "RF", "ST", "BF"]
EMG_HEADER_ACT = ["ACT TA", "ACT GM", "ACT GL", "ACT VM", "ACT VL", "ACT RF", "ACT ST", "ACT BF"]

# The names of the bones in the lower body
COLS_BONES = ["LFE_TX", "LFE_TY", "LFE_TZ", "LFO_TX", "LFO_TY", "LFO_TZ", "LTI_TX", "LTI_TY", "LTI_TZ", 
              "LTO_TX", "LTO_TY", "LTO_TZ", "PEL_TX", "PEL_TY", "PEL_TZ", "RFE_TX", "RFE_TY", "RFE_TZ", 
              "RFO_TX", "RFO_TY", "RFO_TZ", "RTI_TX", "RTI_TY", "RTI_TZ", "RTO_TX", "RTO_TY", "RTO_TZ"]

# The names of the right heel and toe marker columns
R_FOOT_MARKERS = ["RHEEX", "RHEEY", "RHEEZ", "RTOEX", "RTOEY", "RTOEZ"]


##################################
### Simulated IMU header stuff ###
##################################

# Headers of available data columns related to IMU features (in global angle files), list of lists (per body segment)
COLS_SEGMENTS = [["Root_RX", "Root_RY", "Root_RZ", "Root_TX", "Root_TY", "Root_TZ"], 
                 ["L_Femur_RX", "L_Femur_RY", "L_Femur_RZ", "L_Femur_TX", "L_Femur_TY", "L_Femur_TZ"], 
                 ["L_Tibia_RX", "L_Tibia_RY", "L_Tibia_RZ", "L_Tibia_TX", "L_Tibia_TY", "L_Tibia_TZ"], 
                 ["L_Foot_RX", "L_Foot_RY", "L_Foot_RZ", "L_Foot_TX", "L_Foot_TY", "L_Foot_TZ"], 
                 ["R_Femur_RX", "R_Femur_RY", "R_Femur_RZ", "R_Femur_TX", "R_Femur_TY", "R_Femur_TZ"],
                 ["R_Tibia_RX", "R_Tibia_RY", "R_Tibia_RZ", "R_Tibia_TX", "R_Tibia_TY", "R_Tibia_TZ"], 
                 ["R_Foot_RX", "R_Foot_RY", "R_Foot_RZ", "R_Foot_TX", "R_Foot_TY", "R_Foot_TZ"]]
COLS_SEGMENTS_FLAT = [
    COL
    for COLS in COLS_SEGMENTS
    for COL in COLS
]

# Headers of calculated/simulated IMU features
COLS_IMU_FEAT = [["Root_DRX", "Root_DRY", "Root_DRZ", "Root_DDTX", "Root_DDTY", "Root_DDTZ"], 
                 ["L_Femur_DRX", "L_Femur_DRY", "L_Femur_DRZ", "L_Femur_DDTX", "L_Femur_DDTY", "L_Femur_DDTZ"], 
                 ["L_Tibia_DRX", "L_Tibia_DRY", "L_Tibia_DRZ", "L_Tibia_DDTX", "L_Tibia_DDTY", "L_Tibia_DDTZ"], 
                 ["L_Foot_DRX", "L_Foot_DRY", "L_Foot_DRZ", "L_Foot_DDTX", "L_Foot_DDTY", "L_Foot_DDTZ"], 
                 ["R_Femur_DRX", "R_Femur_DRY", "R_Femur_DRZ", "R_Femur_DDTX", "R_Femur_DDTY", "R_Femur_DDTZ"],
                 ["R_Tibia_DRX", "R_Tibia_DRY", "R_Tibia_DRZ", "R_Tibia_DDTX", "R_Tibia_DDTY", "R_Tibia_DDTZ"], 
                 ["R_Foot_DRX", "R_Foot_DRY", "R_Foot_DRZ", "R_Foot_DDTX", "R_Foot_DDTY", "R_Foot_DDTZ"]]
COLS_IMU_FEAT_FLAT = [
    COL
    for COLS in COLS_IMU_FEAT
    for COL in COLS
]
N_IMU_FEAT = len(COLS_IMU_FEAT_FLAT)


#############################################
### Simulated encoder/goniometer features ###
#############################################

# How they are named in the model output files
COLS_JOINT_ANGLES = ["LAnkleAnglesX", "LAnkleAnglesY", "LAnkleAnglesZ", "LKneeAnglesX", "LKneeAnglesY", "LKneeAnglesZ", "LHipAnglesX", "LHipAnglesY", 
                     "LHipAnglesZ", "RAnkleAnglesX", "RAnkleAnglesY", "RAnkleAnglesZ", "RKneeAnglesX", "RKneeAnglesY", "RKneeAnglesZ", "RHipAnglesX", 
                     "RHipAnglesY", "RHipAnglesZ"]

# New names: Angles -> Angle
COLS_JOINT_ANGLE = ["LAnkleAngleX", "LAnkleAngleY", "LAnkleAngleZ", "LKneeAngleX", "LKneeAngleY", "LKneeAngleZ", "LHipAngleX", "LHipAngleY", 
                     "LHipAngleZ", "RAnkleAngleX", "RAnkleAngleY", "RAnkleAngleZ", "RKneeAngleX", "RKneeAngleY", "RKneeAngleZ", "RHipAngleX", 
                     "RHipAngleY", "RHipAngleZ"]
COLS_JOINT_VELOCITY = ["LAnkleAngleDX", "LAnkleAngleDY", "LAnkleAngleDZ", "LKneeAngleDX", "LKneeAngleDY", "LKneeAngleDZ", "LHipAngleDX", 
                       "LHipAngleDY", "LHipAngleDZ", "RAnkleAngleDX", "RAnkleAngleDY", "RAnkleAngleDZ", "RKneeAngleDX", "RKneeAngleDY", 
                       "RKneeAngleDZ", "RHipAngleDX", "RHipAngleDY", "RHipAngleDZ"]


############################
### TRAINING DATA HEADER ###
############################

COLS_TRAINING_DATA = ["Label", "GaitPhase", "GaitSpeed", "NewFile", 
                      
                      "LTOTF", "LCOPX", "LCOPY", "RTOTF", "RCOPX", "RCOPY", "LP1", "LP2", "LP3", "LP4", "LP5", "LP6", "LP7", "LP8", 
                      "LP9", "LP10", "LP11", "LP12", "LP13", "LP14", "LP15", "LP16", "LACCX", "LACCY", "LACCZ", "LANGX", "LANGY", "LANGZ", 
                      "RP1", "RP2", "RP3", "RP4", "RP5", "RP6", "RP7", "RP8", "RP9", "RP10", "RP11", "RP12", "RP13", "RP14", "RP15", "RP16", 
                      "RACCX", "RACCY", "RACCZ", "RANGX", "RANGY", "RANGZ", 
                      "LTOTF 50ms", "LCOPX 50ms", "LCOPY 50ms", "RTOTF 50ms", "RCOPX 50ms", "RCOPY 50ms", "LP1 50ms", "LP2 50ms", "LP3 50ms", "LP4 50ms", 
                      "LP5 50ms", "LP6 50ms", "LP7 50ms", "LP8 50ms", "LP9 50ms", "LP10 50ms", "LP11 50ms", "LP12 50ms", "LP13 50ms", "LP14 50ms", 
                      "LP15 50ms", "LP16 50ms", "LACCX 50ms", "LACCY 50ms", "LACCZ 50ms", "LANGX 50ms", "LANGY 50ms", "LANGZ 50ms", "RP1 50ms", "RP2 50ms", 
                      "RP3 50ms", "RP4 50ms", "RP5 50ms", "RP6 50ms", "RP7 50ms", "RP8 50ms", "RP9 50ms", "RP10 50ms", "RP11 50ms", "RP12 50ms", "RP13 50ms", 
                      "RP14 50ms", "RP15 50ms", "RP16 50ms", "RACCX 50ms", "RACCY 50ms", "RACCZ 50ms", "RANGX 50ms", "RANGY 50ms", "RANGZ 50ms", 
                      "LTOTF 100ms", "LCOPX 100ms", "LCOPY 100ms", "RTOTF 100ms", "RCOPX 100ms", "RCOPY 100ms", "LP1 100ms", "LP2 100ms", "LP3 100ms", 
                      "LP4 100ms", "LP5 100ms", "LP6 100ms", "LP7 100ms", "LP8 100ms", "LP9 100ms", "LP10 100ms", "LP11 100ms", "LP12 100ms", "LP13 100ms", 
                      "LP14 100ms", "LP15 100ms", "LP16 100ms", "LACCX 100ms", "LACCY 100ms", "LACCZ 100ms", "LANGX 100ms", "LANGY 100ms", "LANGZ 100ms", 
                      "RP1 100ms", "RP2 100ms", "RP3 100ms", "RP4 100ms", "RP5 100ms", "RP6 100ms", "RP7 100ms", "RP8 100ms", "RP9 100ms", "RP10 100ms", 
                      "RP11 100ms", "RP12 100ms", "RP13 100ms", "RP14 100ms", "RP15 100ms", "RP16 100ms", "RACCX 100ms", "RACCY 100ms", "RACCZ 100ms", 
                      "RANGX 100ms", "RANGY 100ms", "RANGZ 100ms", 
                      "LTOTF 150ms", "LCOPX 150ms", "LCOPY 150ms", "RTOTF 150ms", "RCOPX 150ms", "RCOPY 150ms", "LP1 150ms", "LP2 150ms", "LP3 150ms", 
                      "LP4 150ms", "LP5 150ms", "LP6 150ms", "LP7 150ms", "LP8 150ms", "LP9 150ms", "LP10 150ms", "LP11 150ms", "LP12 150ms", "LP13 150ms", 
                      "LP14 150ms", "LP15 150ms", "LP16 150ms", "LACCX 150ms", "LACCY 150ms", "LACCZ 150ms", "LANGX 150ms", "LANGY 150ms", "LANGZ 150ms", 
                      "RP1 150ms", "RP2 150ms", "RP3 150ms", "RP4 150ms", "RP5 150ms", "RP6 150ms", "RP7 150ms", "RP8 150ms", "RP9 150ms", "RP10 150ms", 
                      "RP11 150ms", "RP12 150ms", "RP13 150ms", "RP14 150ms", "RP15 150ms", "RP16 150ms", "RACCX 150ms", "RACCY 150ms", "RACCZ 150ms", 
                      "RANGX 150ms", "RANGY 150ms", "RANGZ 150ms", 
                      "LTOTF 200ms", "LCOPX 200ms", "LCOPY 200ms", "RTOTF 200ms", "RCOPX 200ms", "RCOPY 200ms", "LP1 200ms", "LP2 200ms", "LP3 200ms", 
                      "LP4 200ms", "LP5 200ms", "LP6 200ms", "LP7 200ms", "LP8 200ms", "LP9 200ms", "LP10 200ms", "LP11 200ms", "LP12 200ms", "LP13 200ms", 
                      "LP14 200ms", "LP15 200ms", "LP16 200ms", "LACCX 200ms", "LACCY 200ms", "LACCZ 200ms", "LANGX 200ms", "LANGY 200ms", "LANGZ 200ms", 
                      "RP1 200ms", "RP2 200ms", "RP3 200ms", "RP4 200ms", "RP5 200ms", "RP6 200ms", "RP7 200ms", "RP8 200ms", "RP9 200ms", "RP10 200ms", 
                      "RP11 200ms", "RP12 200ms", "RP13 200ms", "RP14 200ms", "RP15 200ms", "RP16 200ms", "RACCX 200ms", "RACCY 200ms", "RACCZ 200ms", 
                      "RANGX 200ms", "RANGY 200ms", "RANGZ 200ms", 
                      
                      "Root_DRX", "Root_DRY", 
                      "Root_DRZ", "Root_DDTX", "Root_DDTY", "Root_DDTZ", "L_Femur_DRX", "L_Femur_DRY", "L_Femur_DRZ", "L_Femur_DDTX", "L_Femur_DDTY", 
                      "L_Femur_DDTZ", "L_Tibia_DRX", "L_Tibia_DRY", "L_Tibia_DRZ", "L_Tibia_DDTX", "L_Tibia_DDTY", "L_Tibia_DDTZ", "L_Foot_DRX", "L_Foot_DRY", 
                      "L_Foot_DRZ", "L_Foot_DDTX", "L_Foot_DDTY", "L_Foot_DDTZ", "R_Femur_DRX", "R_Femur_DRY", "R_Femur_DRZ", "R_Femur_DDTX", "R_Femur_DDTY", 
                      "R_Femur_DDTZ", "R_Tibia_DRX", "R_Tibia_DRY", "R_Tibia_DRZ", "R_Tibia_DDTX", "R_Tibia_DDTY", "R_Tibia_DDTZ", "R_Foot_DRX", "R_Foot_DRY", 
                      "R_Foot_DRZ", "R_Foot_DDTX", "R_Foot_DDTY", "R_Foot_DDTZ", 
                      "Root_DRX 50ms", "Root_DRY 50ms", "Root_DRZ 50ms", "Root_DDTX 50ms", "Root_DDTY 50ms", "Root_DDTZ 50ms", "L_Femur_DRX 50ms", 
                      "L_Femur_DRY 50ms", "L_Femur_DRZ 50ms", "L_Femur_DDTX 50ms", "L_Femur_DDTY 50ms", "L_Femur_DDTZ 50ms", "L_Tibia_DRX 50ms", 
                      "L_Tibia_DRY 50ms", "L_Tibia_DRZ 50ms", "L_Tibia_DDTX 50ms", "L_Tibia_DDTY 50ms", "L_Tibia_DDTZ 50ms", "L_Foot_DRX 50ms", 
                      "L_Foot_DRY 50ms", "L_Foot_DRZ 50ms", "L_Foot_DDTX 50ms", "L_Foot_DDTY 50ms", "L_Foot_DDTZ 50ms", "R_Femur_DRX 50ms", 
                      "R_Femur_DRY 50ms", "R_Femur_DRZ 50ms", "R_Femur_DDTX 50ms", "R_Femur_DDTY 50ms", "R_Femur_DDTZ 50ms", "R_Tibia_DRX 50ms", 
                      "R_Tibia_DRY 50ms", "R_Tibia_DRZ 50ms", "R_Tibia_DDTX 50ms", "R_Tibia_DDTY 50ms", "R_Tibia_DDTZ 50ms", "R_Foot_DRX 50ms", 
                      "R_Foot_DRY 50ms", "R_Foot_DRZ 50ms", "R_Foot_DDTX 50ms", "R_Foot_DDTY 50ms", "R_Foot_DDTZ 50ms", 
                      "Root_DRX 100ms", "Root_DRY 100ms", "Root_DRZ 100ms", "Root_DDTX 100ms", "Root_DDTY 100ms", "Root_DDTZ 100ms", 
                      "L_Femur_DRX 100ms", "L_Femur_DRY 100ms", "L_Femur_DRZ 100ms", "L_Femur_DDTX 100ms", "L_Femur_DDTY 100ms", "L_Femur_DDTZ 100ms", 
                      "L_Tibia_DRX 100ms", "L_Tibia_DRY 100ms", "L_Tibia_DRZ 100ms", "L_Tibia_DDTX 100ms", "L_Tibia_DDTY 100ms", "L_Tibia_DDTZ 100ms", 
                      "L_Foot_DRX 100ms", "L_Foot_DRY 100ms", "L_Foot_DRZ 100ms", "L_Foot_DDTX 100ms", "L_Foot_DDTY 100ms", "L_Foot_DDTZ 100ms", 
                      "R_Femur_DRX 100ms", "R_Femur_DRY 100ms", "R_Femur_DRZ 100ms", "R_Femur_DDTX 100ms", "R_Femur_DDTY 100ms", "R_Femur_DDTZ 100ms", 
                      "R_Tibia_DRX 100ms", "R_Tibia_DRY 100ms", "R_Tibia_DRZ 100ms", "R_Tibia_DDTX 100ms", "R_Tibia_DDTY 100ms", "R_Tibia_DDTZ 100ms", 
                      "R_Foot_DRX 100ms", "R_Foot_DRY 100ms", "R_Foot_DRZ 100ms", "R_Foot_DDTX 100ms", "R_Foot_DDTY 100ms", "R_Foot_DDTZ 100ms", 
                      "Root_DRX 150ms", "Root_DRY 150ms", "Root_DRZ 150ms", "Root_DDTX 150ms", "Root_DDTY 150ms", "Root_DDTZ 150ms", 
                      "L_Femur_DRX 150ms", "L_Femur_DRY 150ms", "L_Femur_DRZ 150ms", "L_Femur_DDTX 150ms", "L_Femur_DDTY 150ms", "L_Femur_DDTZ 150ms", 
                      "L_Tibia_DRX 150ms", "L_Tibia_DRY 150ms", "L_Tibia_DRZ 150ms", "L_Tibia_DDTX 150ms", "L_Tibia_DDTY 150ms", "L_Tibia_DDTZ 150ms", 
                      "L_Foot_DRX 150ms", "L_Foot_DRY 150ms", "L_Foot_DRZ 150ms", "L_Foot_DDTX 150ms", "L_Foot_DDTY 150ms", "L_Foot_DDTZ 150ms", 
                      "R_Femur_DRX 150ms", "R_Femur_DRY 150ms", "R_Femur_DRZ 150ms", "R_Femur_DDTX 150ms", "R_Femur_DDTY 150ms", "R_Femur_DDTZ 150ms", 
                      "R_Tibia_DRX 150ms", "R_Tibia_DRY 150ms", "R_Tibia_DRZ 150ms", "R_Tibia_DDTX 150ms", "R_Tibia_DDTY 150ms", "R_Tibia_DDTZ 150ms", 
                      "R_Foot_DRX 150ms", "R_Foot_DRY 150ms", "R_Foot_DRZ 150ms", "R_Foot_DDTX 150ms", "R_Foot_DDTY 150ms", "R_Foot_DDTZ 150ms", 
                      "Root_DRX 200ms", "Root_DRY 200ms", "Root_DRZ 200ms", "Root_DDTX 200ms", "Root_DDTY 200ms", "Root_DDTZ 200ms", 
                      "L_Femur_DRX 200ms", "L_Femur_DRY 200ms", "L_Femur_DRZ 200ms", "L_Femur_DDTX 200ms", "L_Femur_DDTY 200ms", "L_Femur_DDTZ 200ms", 
                      "L_Tibia_DRX 200ms", "L_Tibia_DRY 200ms", "L_Tibia_DRZ 200ms", "L_Tibia_DDTX 200ms", "L_Tibia_DDTY 200ms", "L_Tibia_DDTZ 200ms", 
                      "L_Foot_DRX 200ms", "L_Foot_DRY 200ms", "L_Foot_DRZ 200ms", "L_Foot_DDTX 200ms", "L_Foot_DDTY 200ms", "L_Foot_DDTZ 200ms", 
                      "R_Femur_DRX 200ms", "R_Femur_DRY 200ms", "R_Femur_DRZ 200ms", "R_Femur_DDTX 200ms", "R_Femur_DDTY 200ms", "R_Femur_DDTZ 200ms", 
                      "R_Tibia_DRX 200ms", "R_Tibia_DRY 200ms", "R_Tibia_DRZ 200ms", "R_Tibia_DDTX 200ms", "R_Tibia_DDTY 200ms", "R_Tibia_DDTZ 200ms", 
                      "R_Foot_DRX 200ms", "R_Foot_DRY 200ms", "R_Foot_DRZ 200ms", "R_Foot_DDTX 200ms", "R_Foot_DDTY 200ms", "R_Foot_DDTZ 200ms", 
                      
                      "LAnkleAngleX", "LAnkleAngleY", "LAnkleAngleZ", "LKneeAngleX", "LKneeAngleY", "LKneeAngleZ", "LHipAngleX", "LHipAngleY", "LHipAngleZ", 
                      "RAnkleAngleX", "RAnkleAngleY", "RAnkleAngleZ", "RKneeAngleX", "RKneeAngleY", "RKneeAngleZ", "RHipAngleX", "RHipAngleY", "RHipAngleZ", 
                      "LAnkleAngleDX", "LAnkleAngleDY", "LAnkleAngleDZ", "LKneeAngleDX", "LKneeAngleDY", "LKneeAngleDZ", "LHipAngleDX", "LHipAngleDY", 
                      "LHipAngleDZ", "RAnkleAngleDX", "RAnkleAngleDY", "RAnkleAngleDZ", "RKneeAngleDX", "RKneeAngleDY", "RKneeAngleDZ", "RHipAngleDX", 
                      "RHipAngleDY", "RHipAngleDZ", 
                      "LAnkleAngleX 50ms", "LAnkleAngleY 50ms", "LAnkleAngleZ 50ms", "LKneeAngleX 50ms", "LKneeAngleY 50ms", "LKneeAngleZ 50ms", 
                      "LHipAngleX 50ms", "LHipAngleY 50ms", "LHipAngleZ 50ms", "RAnkleAngleX 50ms", "RAnkleAngleY 50ms", "RAnkleAngleZ 50ms", 
                      "RKneeAngleX 50ms", "RKneeAngleY 50ms", "RKneeAngleZ 50ms", "RHipAngleX 50ms", "RHipAngleY 50ms", "RHipAngleZ 50ms", 
                      "LAnkleAngleDX 50ms", "LAnkleAngleDY 50ms", "LAnkleAngleDZ 50ms", "LKneeAngleDX 50ms", "LKneeAngleDY 50ms", "LKneeAngleDZ 50ms", 
                      "LHipAngleDX 50ms", "LHipAngleDY 50ms", "LHipAngleDZ 50ms", "RAnkleAngleDX 50ms", "RAnkleAngleDY 50ms", "RAnkleAngleDZ 50ms", 
                      "RKneeAngleDX 50ms", "RKneeAngleDY 50ms", "RKneeAngleDZ 50ms", "RHipAngleDX 50ms", "RHipAngleDY 50ms", "RHipAngleDZ 50ms", 
                      "LAnkleAngleX 100ms", "LAnkleAngleY 100ms", "LAnkleAngleZ 100ms", "LKneeAngleX 100ms", "LKneeAngleY 100ms", "LKneeAngleZ 100ms", 
                      "LHipAngleX 100ms", "LHipAngleY 100ms", "LHipAngleZ 100ms", "RAnkleAngleX 100ms", "RAnkleAngleY 100ms", "RAnkleAngleZ 100ms", 
                      "RKneeAngleX 100ms", "RKneeAngleY 100ms", "RKneeAngleZ 100ms", "RHipAngleX 100ms", "RHipAngleY 100ms", "RHipAngleZ 100ms", 
                      "LAnkleAngleDX 100ms", "LAnkleAngleDY 100ms", "LAnkleAngleDZ 100ms", "LKneeAngleDX 100ms", "LKneeAngleDY 100ms", "LKneeAngleDZ 100ms", 
                      "LHipAngleDX 100ms", "LHipAngleDY 100ms", "LHipAngleDZ 100ms", "RAnkleAngleDX 100ms", "RAnkleAngleDY 100ms", "RAnkleAngleDZ 100ms", 
                      "RKneeAngleDX 100ms", "RKneeAngleDY 100ms", "RKneeAngleDZ 100ms", "RHipAngleDX 100ms", "RHipAngleDY 100ms", "RHipAngleDZ 100ms", 
                      "LAnkleAngleX 150ms", "LAnkleAngleY 150ms", "LAnkleAngleZ 150ms", "LKneeAngleX 150ms", "LKneeAngleY 150ms", "LKneeAngleZ 150ms", 
                      "LHipAngleX 150ms", "LHipAngleY 150ms", "LHipAngleZ 150ms", "RAnkleAngleX 150ms", "RAnkleAngleY 150ms", "RAnkleAngleZ 150ms", 
                      "RKneeAngleX 150ms", "RKneeAngleY 150ms", "RKneeAngleZ 150ms", "RHipAngleX 150ms", "RHipAngleY 150ms", "RHipAngleZ 150ms", 
                      "LAnkleAngleDX 150ms", "LAnkleAngleDY 150ms", "LAnkleAngleDZ 150ms", "LKneeAngleDX 150ms", "LKneeAngleDY 150ms", "LKneeAngleDZ 150ms", 
                      "LHipAngleDX 150ms", "LHipAngleDY 150ms", "LHipAngleDZ 150ms", "RAnkleAngleDX 150ms", "RAnkleAngleDY 150ms", "RAnkleAngleDZ 150ms", 
                      "RKneeAngleDX 150ms", "RKneeAngleDY 150ms", "RKneeAngleDZ 150ms", "RHipAngleDX 150ms", "RHipAngleDY 150ms", "RHipAngleDZ 150ms", 
                      "LAnkleAngleX 200ms", "LAnkleAngleY 200ms", "LAnkleAngleZ 200ms", "LKneeAngleX 200ms", "LKneeAngleY 200ms", "LKneeAngleZ 200ms", 
                      "LHipAngleX 200ms", "LHipAngleY 200ms", "LHipAngleZ 200ms", "RAnkleAngleX 200ms", "RAnkleAngleY 200ms", "RAnkleAngleZ 200ms", 
                      "RKneeAngleX 200ms", "RKneeAngleY 200ms", "RKneeAngleZ 200ms", "RHipAngleX 200ms", "RHipAngleY 200ms", "RHipAngleZ 200ms", 
                      "LAnkleAngleDX 200ms", "LAnkleAngleDY 200ms", "LAnkleAngleDZ 200ms", "LKneeAngleDX 200ms", "LKneeAngleDY 200ms", "LKneeAngleDZ 200ms", 
                      "LHipAngleDX 200ms", "LHipAngleDY 200ms", "LHipAngleDZ 200ms", "RAnkleAngleDX 200ms", "RAnkleAngleDY 200ms", "RAnkleAngleDZ 200ms", 
                      "RKneeAngleDX 200ms", "RKneeAngleDY 200ms", "RKneeAngleDZ 200ms", "RHipAngleDX 200ms", "RHipAngleDY 200ms", "RHipAngleDZ 200ms", 
                      
                      "RMS TA", "RMS GM", "RMS GL", "RMS VM", "RMS VL", "RMS RF", "RMS ST", "RMS BF", "MNF TA", "MNF GM", "MNF GL", "MNF VM", "MNF VL", 
                      "MNF RF", "MNF ST", "MNF BF", "ACT TA", "ACT GM", "ACT GL", "ACT VM", "ACT VL", "ACT RF", "ACT ST", "ACT BF", 
                      "RMS TA 50ms", "RMS GM 50ms", "RMS GL 50ms", "RMS VM 50ms", "RMS VL 50ms", "RMS RF 50ms", "RMS ST 50ms", "RMS BF 50ms", "MNF TA 50ms", 
                      "MNF GM 50ms", "MNF GL 50ms", "MNF VM 50ms", "MNF VL 50ms", "MNF RF 50ms", "MNF ST 50ms", "MNF BF 50ms", "ACT TA 50ms", "ACT GM 50ms", 
                      "ACT GL 50ms", "ACT VM 50ms", "ACT VL 50ms", "ACT RF 50ms", "ACT ST 50ms", "ACT BF 50ms", 
                      "RMS TA 100ms", "RMS GM 100ms", "RMS GL 100ms", "RMS VM 100ms", "RMS VL 100ms", "RMS RF 100ms", "RMS ST 100ms", "RMS BF 100ms", 
                      "MNF TA 100ms", "MNF GM 100ms", "MNF GL 100ms", "MNF VM 100ms", "MNF VL 100ms", "MNF RF 100ms", "MNF ST 100ms", "MNF BF 100ms", 
                      "ACT TA 100ms", "ACT GM 100ms", "ACT GL 100ms", "ACT VM 100ms", "ACT VL 100ms", "ACT RF 100ms", "ACT ST 100ms", "ACT BF 100ms", 
                      "RMS TA 150ms", "RMS GM 150ms", "RMS GL 150ms", "RMS VM 150ms", "RMS VL 150ms", "RMS RF 150ms", "RMS ST 150ms", "RMS BF 150ms", 
                      "MNF TA 150ms", "MNF GM 150ms", "MNF GL 150ms", "MNF VM 150ms", "MNF VL 150ms", "MNF RF 150ms", "MNF ST 150ms", "MNF BF 150ms", 
                      "ACT TA 150ms", "ACT GM 150ms", "ACT GL 150ms", "ACT VM 150ms", "ACT VL 150ms", "ACT RF 150ms", "ACT ST 150ms", "ACT BF 150ms", 
                      "RMS TA 200ms", "RMS GM 200ms", "RMS GL 200ms", "RMS VM 200ms", "RMS VL 200ms", "RMS RF 200ms", "RMS ST 200ms", "RMS BF 200ms", 
                      "MNF TA 200ms", "MNF GM 200ms", "MNF GL 200ms", "MNF VM 200ms", "MNF VL 200ms", "MNF RF 200ms", "MNF ST 200ms", "MNF BF 200ms", 
                      "ACT TA 200ms", "ACT GM 200ms", "ACT GL 200ms", "ACT VM 200ms", "ACT VL 200ms", "ACT RF 200ms", "ACT ST 200ms", "ACT BF 200ms"]


######################
### Activity stuff ###
######################

# Static Activities
SA = ["Sitting", "Standing"]
SA_TEST = ["Sitting obstacle", "Standing weight", "Standing obstacle"]
# SA_TEST = ["Sitting", "Sitting obstacle", "Standing", "Standing weight", "Standing obstacle"]

# Transitional Activities
TA = ["Calf raises", "Squat", "Sit-to-stand and stand-to-sit"]
TA_TEST = ["Calf raises weight", "Calf raises obstacle", "Squat weight", "Squat obstacle", 
           "Sit-to-stand and stand-to-sit weight", "Sit-to-stand and stand-to-sit obstacle"]
# TA_TEST = ["Calf raises", "Calf raises weight", "Calf raises obstacle", "Squat", "Squat weight", "Squat obstacle", 
#            "Sit-to-stand and stand-to-sit", "Sit-to-stand and stand-to-sit weight", "Sit-to-stand and stand-to-sit obstacle"]

# Symmetric Rhythmic Activities
SRA = ["Normal level walking", "Slow level walking", "Fast level walking", "Jogging", "Walking backwards", 
       "Stair ascent (SOS)(L first)", "Stair descent (SOS)(L first)", "Stair ascent (SOS)(R first)", 
       "Stair descent (SOS)(R first)", "Ramp ascent (10deg)", "Ramp descent (10deg)"]
SRA_TEST = ["Normal level walking weight", "Normal level walking obstacle", "Slow level walking weight", 
            "Slow level walking obstacle", "Fast level walking weight", "Fast level walking obstacle", "Jogging weight", 
            "Jogging obstacle", "Walking backwards weight", "Walking backwards obstacle", 
            "Stair ascent (SOS)(L first) weight", "Stair ascent (SOS)(L first) obstacle", 
            "Stair descent (SOS)(L first) weight", "Stair descent (SOS)(L first) obstacle", 
            "Stair ascent (SOS)(R first) weight", "Stair ascent (SOS)(R first) obstacle", 
            "Stair descent (SOS)(R first) weight", "Stair descent (SOS)(R first) obstacle", "Ramp ascent (10deg) weight", 
            "Ramp ascent (10deg) obstacle", "Ramp descent (10deg) weight", "Ramp descent (10deg) obstacle"]
# SRA_TEST = ["Normal level walking", "Normal level walking weight", "Normal level walking obstacle", "Slow level walking", "Slow level walking weight", 
#             "Slow level walking obstacle", "Fast level walking", "Fast level walking weight", "Fast level walking obstacle", "Jogging", "Jogging weight", 
#             "Jogging obstacle", "Walking backwards", "Walking backwards weight", "Walking backwards obstacle", "Stair ascent (SOS)(L first)", 
#             "Stair ascent (SOS)(L first) weight", "Stair ascent (SOS)(L first) obstacle", "Stair descent (SOS)(L first)", 
#             "Stair descent (SOS)(L first) weight", "Stair descent (SOS)(L first) obstacle", "Stair ascent (SOS)(R first)", 
#             "Stair ascent (SOS)(R first) weight", "Stair ascent (SOS)(R first) obstacle", "Stair descent (SOS)(R first)", 
#             "Stair descent (SOS)(R first) weight", "Stair descent (SOS)(R first) obstacle", "Ramp ascent (10deg)", "Ramp ascent (10deg) weight", 
#             "Ramp ascent (10deg) obstacle", "Ramp descent (10deg)", "Ramp descent (10deg) weight", "Ramp descent (10deg) obstacle"]

# Asymmetric Rhythmic Activities
ARA = ["Circular walking (diam 2m)(CCW)", "Circular walking (diam 2m)(CW)", "Side stepping (L)", "Side stepping (R)", 
       "Stair ascent (SBS)(L first)", "Stair ascent (SBS)(R first)", "Stair descent (SBS)(L first)", "Stair descent (SBS)(R first)"]
ARA_TRAIN = ["Circular walking (diam 2m)(CCW)", "Circular walking (diam 2m)(CW)", "Turning on the spot (CW)", "Turning on the spot (CCW)", 
             "Side stepping (L)", "Side stepping (R)", "Stair ascent (SBS)(L first)", "Stair ascent (SBS)(R first)", 
             "Stair descent (SBS)(L first)", "Stair descent (SBS)(R first)"]
ARA_TEST = ["Turning on the spot (CW) weight", "Turning on the spot (CCW) weight", 
            "Side stepping (L) weight", "Side stepping (R) weight", "Circular walking (diam 2m)(CW) weight", "Circular walking (diam 2m)(CCW) weight", 
            "Side stepping (L) obstacle", "Side stepping (R) obstacle", "Circular walking (diam 2m)(CW) obstacle", "Circular walking (diam 2m)(CCW) obstacle", 
            "Stair ascent (SBS)(L first) weight", "Stair descent (SBS)(L first) weight", "Stair ascent (SBS)(R first) weight", 
            "Stair descent (SBS)(R first) weight", "Stair ascent (SBS)(L first) obstacle", "Stair descent (SBS)(L first) obstacle", 
            "Stair ascent (SBS)(R first) obstacle", "Stair descent (SBS)(R first) obstacle"]
# ARA_TEST = ["Circular walking (diam 2m)(CCW)", "Circular walking (diam 2m)(CW)", "Turning on the spot (CW)", "Turning on the spot (CCW)", 
#             "Side stepping (L)", "Side stepping (R)", "Stair ascent (SBS)(L first)", "Stair ascent (SBS)(R first)", 
#             "Stair descent (SBS)(L first)", "Stair descent (SBS)(R first)", "Turning on the spot (CW) weight", "Turning on the spot (CCW) weight", 
#             "Side stepping (L) weight", "Side stepping (R) weight", "Circular walking (diam 2m)(CW) weight", "Circular walking (diam 2m)(CCW) weight", 
#             "Side stepping (L) obstacle", "Side stepping (R) obstacle", "Circular walking (diam 2m)(CW) obstacle", "Circular walking (diam 2m)(CCW) obstacle", 
#             "Stair ascent (SBS)(L first) weight", "Stair descent (SBS)(L first) weight", "Stair ascent (SBS)(R first) weight", 
#             "Stair descent (SBS)(R first) weight", "Stair ascent (SBS)(L first) obstacle", "Stair descent (SBS)(L first) obstacle", 
#             "Stair ascent (SBS)(R first) obstacle", "Stair descent (SBS)(R first) obstacle"]

# Combination of activities in one file
COMBINATIONS = ["Sit-to-stand-to-walk", "Combination 1", "Combination 2", "Combination 3"]

# Activities that should not be processed for reference trajectory calculations
IGNORE_REF = ["Normal level walking weight", "Slow level walking weight", "Fast level walking weight", "Jogging weight", "Walking backwards weight", 
              "Turning on the spot (CW) weight", "Turning on the spot (CCW) weight", "Side stepping (L) weight", "Side stepping (R) weight", 
              "Circular walking (diam 2m)(CW) weight", "Circular walking (diam 2m)(CCW) weight", "Normal level walking obstacle", 
              "Slow level walking obstacle", "Fast level walking obstacle", "Jogging obstacle", "Walking backwards obstacle", "Side stepping (L) obstacle", 
              "Side stepping (R) obstacle", "Circular walking (diam 2m)(CW) obstacle", "Circular walking (diam 2m)(CCW) obstacle", 
              "Stair ascent (SOS)(L first) weight", "Stair descent (SOS)(L first) weight", "Stair ascent (SOS)(R first) weight", 
              "Stair descent (SOS)(R first) weight", "Stair ascent (SBS)(L first) weight", "Stair descent (SBS)(L first) weight", 
              "Stair ascent (SBS)(R first) weight", "Stair descent (SBS)(R first) weight", "Ramp ascent (10deg) weight", "Ramp descent (10deg) weight",
              "Standing weight", "Squat weight", "Calf raises weight", "Sit-to-stand and stand-to-sit weight", "Stair ascent (SOS)(L first) obstacle", 
              "Stair descent (SOS)(L first) obstacle", "Stair ascent (SOS)(R first) obstacle", "Stair descent (SOS)(R first) obstacle", 
              "Stair ascent (SBS)(L first) obstacle", "Stair descent (SBS)(L first) obstacle", "Stair ascent (SBS)(R first) obstacle", 
              "Stair descent (SBS)(R first) obstacle", "Ramp ascent (10deg) obstacle", "Ramp descent (10deg) obstacle", "Standing obstacle", "Squat obstacle", 
              "Calf raises obstacle", "Sitting obstacle", "Sit-to-stand and stand-to-sit obstacle", 
              "Sit-to-stand-to-walk", "Combination 1", "Combination 2", "Combination 3", "Forward walking with 360 turn", 
              "Forward walking x stop x backward walking", "Forward walking with 180 turn", "Backward walking x stop x forward walking",
              "Dynamic side stepping", "Diagonally forward walking (strafing)", "Diagonally backward walking (strafing)", "Stair ascent backward", 
              "Stair ascent forward x stop x stair descent backward", "Stair ascent forward x 180 turn x stair descent forward", 
              "Stair ascent backward x stop x stair descent forward", "Stair descent backward", "Stair descent forward x stop x stair ascent backward", 
              "Stair descent forward x 180 turn x stair ascent forward", "Stair descent backward x stop x stair ascent forward", 
              "Ramp ascent x 360 turn x ramp ascent", "Ramp ascent forward x stop x ramp descent backward", 
              "Ramp ascent forward x 180 turn x ramp descent forward", "Ramp ascent backward x stop x ramp descent forward", "Ramp ascent backward", 
              "Ramp descent x 360 turn x ramp descent", "Ramp descent forward x stop x ramp ascent backward", 
              "Ramp descent forward x 180 turn x ramp ascent forward", "Ramp descent backward x stop x ramp ascent forward", "Ramp descent backward"]