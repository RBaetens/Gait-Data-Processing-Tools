%% Start with a clean workspace and complete
%  -----------------------------------------

% Clean
close all
clear
clc

% Add toolbox
addpath(".\EZC3D")


%% Read a C3D file and extract some data in a clear format
%  -------------------------------------------------------

% Filename must be a character string and not an actual string
test_file = 'C:\Users\robaeten\OneDrive - UGent\Full GRF estimation\Gait classification\Subject1\Session2\Recording2.c3d';
test_c3d_struct = ezc3dRead(test_file);

% Extract some subject metadata
subj_name = test_c3d_struct.parameters.SUBJECTS.NAMES.DATA;
subj_mass = test_c3d_struct.parameters.PROCESSING.Bodymass.DATA;
subj_height = test_c3d_struct.parameters.PROCESSING.Height.DATA; % All other subj parameters in PROCESSING as well
marker_set = test_c3d_struct.parameters.SUBJECTS.MARKER_SETS.DATA;

% Extract some marker metadata
f_markers = test_c3d_struct.parameters.POINT.RATE.DATA;
first_frame_markers = test_c3d_struct.header.points.firstFrame;
last_frame_markers = test_c3d_struct.header.points.lastFrame;
n_frames_markers = last_frame_markers-first_frame_markers+1;
n_markers = test_c3d_struct.parameters.POINT.USED.DATA;
marker_labels = test_c3d_struct.parameters.POINT.LABELS.DATA;

% Extract some analog metadata
f_analog = test_c3d_struct.parameters.ANALOG.RATE.DATA;
first_frame_analog = test_c3d_struct.header.analogs.firstFrame;
last_frame_analog = test_c3d_struct.header.analogs.lastFrame;
n_frames_analog = last_frame_analog-first_frame_analog+1;
n_analog_signals = test_c3d_struct.parameters.ANALOG.USED.DATA;
analog_labels = test_c3d_struct.parameters.ANALOG.LABELS.DATA;
analog_descriptions = test_c3d_struct.parameters.ANALOG.DESCRIPTIONS.DATA;

% Extract some force plate metadata
n_force_plates = test_c3d_struct.parameters.FORCE_PLATFORM.USED.DATA;
corners_force_plates = test_c3d_struct.parameters.FORCE_PLATFORM.CORNERS.DATA;
cal_matrix_force_plates = test_c3d_struct.parameters.FORCE_PLATFORM.CAL_MATRIX.DATA; % Empty now !!!

% Extract marker data
markers = test_c3d_struct.data.points; % Of size (3, n_markers, n_frames_markers)

% Extract analog data
analog_signals = test_c3d_struct.data.analogs; % Of size (n_frames_analog, n_analog_signals)

% Construct clear analog labels
clear_analog_labels = analog_descriptions;
for i = 1:size(analog_labels, 1)
    if startsWith(clear_analog_labels{i}, 'Analog EMG::Voltage') || startsWith(clear_analog_labels{i}, 'Generic Analog::Electric Potential')
        clear_analog_labels{i} = '';
    end
    
    clear_analog_labels{i} = replace(clear_analog_labels{i}, 'Kistler Force Plate 2.0.0.0', 'Force plate');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '::Raw', '');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[16]', 'A');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[17]', 'B');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[18]', 'C');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[19]', 'D');
    
    clear_analog_labels{i} = append(clear_analog_labels{i}, ' ', analog_labels{i});
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '.', ' ');
end