%% Start with a clean workspace and complete
%  -----------------------------------------

% Clean
close all
clear
clc

% Add toolbox
addpath(".\EZC3D")
addpath(".\fp_cal_tools")


%% Data locations
%  --------------

dataFolder = '..\Calibration Data\Final measurements';
dataRec1 = strcat(dataFolder, '\FP1B_cal_data_01.c3d');
dataRec2 = strcat(dataFolder, '\FP1B_cal_data_02.c3d'); % validation


%% dataRec1
%  --------

% ---
% Read a C3D file and extract some data in a clear format
% ---

Rec1_c3d_struct = ezc3dRead(dataRec1);

% Extract marker data
markers = Rec1_c3d_struct.data.points; % Of size (3, n_markers, n_frames_markers)
%markers = markers(:, :, 1:29000);

% Extract analog data
analog_signals = Rec1_c3d_struct.data.analogs; % Of size (n_frames_analog, n_analog_signals)
%analog_signals = analog_signals(1:290000, :);

% Construct clear analog labels
analog_labels = Rec1_c3d_struct.parameters.ANALOG.LABELS.DATA;
analog_descriptions = Rec1_c3d_struct.parameters.ANALOG.DESCRIPTIONS.DATA;
clear_analog_labels = analog_descriptions;
for i = 1:size(analog_labels, 1)
    if startsWith(clear_analog_labels{i}, 'Analog EMG::Voltage') || startsWith(clear_analog_labels{i}, 'Generic Analog::Electric Potential')
        clear_analog_labels{i} = '';
    end

    clear_analog_labels{i} = replace(clear_analog_labels{i}, 'Kistler Force Plate 2.0.0.0', 'Force plate');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '::Raw', '');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[15]', 'D');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[12]', 'A');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[13]', 'B');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[14]', 'C');
    
    clear_analog_labels{i} = append(clear_analog_labels{i}, ' ', analog_labels{i});
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '.', ' ');
end

% Get force plate corners
corners_force_plates = Rec1_c3d_struct.parameters.FORCE_PLATFORM.CORNERS.DATA;

% ---
% Get the data you actually need
% ---

% Force plates
PILS_meas_mat = get_fp_fm_meas(analog_signals, 1);
cop_plate = get_fp_cop_meas(analog_signals, corners_force_plates, 1);
Ftot_plate = sqrt(PILS_meas_mat(:, 1).^2 + PILS_meas_mat(:, 2).^2 + PILS_meas_mat(:, 3).^2);

% Rod
PILS_ref_mat = get_pole_fm_meas(analog_signals, markers, corners_force_plates, 1);
tip = get_pole_tip_meas(markers);
Ftot_rod = sqrt(PILS_ref_mat(:, 1).^2 + PILS_ref_mat(:, 2).^2 + PILS_ref_mat(:, 3).^2);

% % ---
% % Compare with plots
% % ---
% 
% figure;
% plot(Ftot_rod, "blue");
% hold on
% plot(Ftot_plate, "red");
% hold off
% 
% figure;
% scatter(tip(:, 1), tip(:, 2), "blue");
% hold on
% scatter(cop_plate(:, 1), cop_plate(:, 2), "red");
% hold off

% ---
% Threshold
% ---

F_thresh = 100;
mask = (Ftot_rod >= F_thresh) & (Ftot_plate >= F_thresh);
n_pts = sum(mask);

R = zeros(6, n_pts);
S = zeros(6, n_pts);
cop_plate_masked = zeros(n_pts, 3);
tip_masked = zeros(n_pts, 3);

j = 1;
for i=1:size(mask, 1)
    if mask(i)
        R(:, j) = PILS_ref_mat(i, :)';
        S(:, j) = PILS_meas_mat(i, :)';
        cop_plate_masked(j, :) = cop_plate(i, :);
        tip_masked(j, :) = tip(i, :);
        j = j+1;
    end
end

% ---
% Compare with plots
% ---

figure;
scatter(tip_masked(:, 1), tip_masked(:, 2), "blue");
hold on
scatter(cop_plate_masked(:, 1), cop_plate_masked(:, 2), "red");
hold off


%% Calculate calibration matrix
%  ----------------------------

C = R * S' * inv(S * S');


%% dataRec2: validation
%  --------------------

% ---
% Read a C3D file and extract some data in a clear format
% ---

Rec2_c3d_struct = ezc3dRead(dataRec2);

% Extract marker data
markers = Rec2_c3d_struct.data.points; % Of size (3, n_markers, n_frames_markers)

% Extract analog data
analog_signals = Rec2_c3d_struct.data.analogs; % Of size (n_frames_analog, n_analog_signals)

% Construct clear analog labels
analog_labels = Rec2_c3d_struct.parameters.ANALOG.LABELS.DATA;
analog_descriptions = Rec2_c3d_struct.parameters.ANALOG.DESCRIPTIONS.DATA;
clear_analog_labels = analog_descriptions;
for i = 1:size(analog_labels, 1)
    if startsWith(clear_analog_labels{i}, 'Analog EMG::Voltage') || startsWith(clear_analog_labels{i}, 'Generic Analog::Electric Potential')
        clear_analog_labels{i} = '';
    end

    clear_analog_labels{i} = replace(clear_analog_labels{i}, 'Kistler Force Plate 2.0.0.0', 'Force plate');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '::Raw', '');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[15]', 'D');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[12]', 'A');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[13]', 'B');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[14]', 'C');
    
    clear_analog_labels{i} = append(clear_analog_labels{i}, ' ', analog_labels{i});
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '.', ' ');
end

% Get force plate corners
corners_force_plates = Rec2_c3d_struct.parameters.FORCE_PLATFORM.CORNERS.DATA;

% ---
% Get the data you actually need
% ---

% Force plates
PILS_meas_mat = get_fp_fm_meas(analog_signals, 1);
cop_plate = get_fp_cop_meas(analog_signals, corners_force_plates, 1);
Ftot_plate = sqrt(PILS_meas_mat(:, 1).^2 + PILS_meas_mat(:, 2).^2 + PILS_meas_mat(:, 3).^2);

% Rod
PILS_ref_mat = get_pole_fm_meas(analog_signals, markers, corners_force_plates, 1);
tip = get_pole_tip_meas(markers);
Ftot_rod = sqrt(PILS_ref_mat(:, 1).^2 + PILS_ref_mat(:, 2).^2 + PILS_ref_mat(:, 3).^2);

% ---
% Threshold
% ---

F_thresh = 100;
mask = (Ftot_rod >= F_thresh) & (Ftot_plate >= F_thresh);
n_pts = sum(mask);

R_val = zeros(6, n_pts);
S_val = zeros(6, n_pts);
cop_plate_masked = zeros(n_pts, 3);
tip_masked = zeros(n_pts, 3);

j = 1;
for i=1:size(mask, 1)
    if mask(i)
        R_val(:, j) = PILS_ref_mat(i, :)';
        S_val(:, j) = PILS_meas_mat(i, :)';
        cop_plate_masked(j, :) = cop_plate(i, :);
        tip_masked(j, :) = tip(i, :);
        j = j+1;
    end
end

% ---
% Get cops from uncorrected force plate data, rod data and corrected force plate data
% ---

S_val_corr = C * S_val;
cop_fp = get_cop_from_fm(S_val, corners_force_plates, 1);
cop_rod = get_cop_from_fm(R_val, corners_force_plates, 1);
cop_fp_corr = get_cop_from_fm(S_val_corr, corners_force_plates, 1);

% ---
% Compare with plots
% ---

figure;
scatter(cop_fp(:, 1), cop_fp(:, 2), "red");
hold on
scatter(cop_plate_masked(:, 1), cop_plate_masked(:, 2), "green");
hold off

figure;
scatter(tip_masked(:, 1), tip_masked(:, 2), "blue");
hold on
scatter(cop_plate_masked(:, 1), cop_plate_masked(:, 2), "red");
hold off

figure;
scatter(cop_rod(:, 1), cop_rod(:, 2), "blue"); % Don't know why there is so much variation in this cop, figure this out!!!
hold on
scatter(cop_fp(:, 1), cop_fp(:, 2), "red");
scatter(cop_fp_corr(:, 1), cop_fp_corr(:, 2), "green"); % Don't know why there is so much variation in this cop, figure this out!!!
hold off


%% Calculate RMSES
%  ---------------

errors_uncorr = zeros(size(cop_fp, 1), 1);
for i=1:size(cop_fp, 1)
    errors_uncorr(i) = norm(cop_fp(i, :) - cop_rod(i, :));
end
RMSE_uncorr = sqrt(mean(errors_uncorr.^2));

errors_corr = zeros(size(cop_fp_corr, 1), 1);
for i=1:size(cop_fp_corr, 1)
    errors_corr(i) = norm(cop_fp_corr(i, :) - cop_rod(i, :));
end
RMSE_corr = sqrt(mean(errors_corr.^2));

RMSE_uncorr_comp = rmse(cop_fp, cop_rod);
RMSE_corr_comp = rmse(cop_fp_corr, cop_rod);


ampl_FM = max(R_val') - min(R_val');
RMSE_FM_uncorr = rmse(S_val', R_val')./ampl_FM;
RMSE_FM_corr = rmse(S_val_corr', R_val')./ampl_FM;


%% Store calibration matrix
%  ------------------------

saveFolder = '..\Calibration Data\Calibration results';
writematrix(C, strcat(saveFolder, '\C_matrix_FP1B.txt'))

