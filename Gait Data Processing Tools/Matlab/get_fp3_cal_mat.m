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
dataRec1 = strcat(dataFolder, '\FP3D_cal_data_01.c3d');
dataRec2 = strcat(dataFolder, '\FP3D_cal_data_02.c3d');
dataRec3 = strcat(dataFolder, '\FP3D_cal_data_03.c3d');
dataRec4 = strcat(dataFolder, '\FP3D_cal_data_04.c3d'); % validation


%% dataRec1
%  --------

% ---
% Read a C3D file and extract some data in a clear format
% ---

Rec1_c3d_struct = ezc3dRead(dataRec1);

% Extract marker data
markers = Rec1_c3d_struct.data.points; % Of size (3, n_markers, n_frames_markers)

% Extract analog data
analog_signals = Rec1_c3d_struct.data.analogs; % Of size (n_frames_analog, n_analog_signals)

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
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[15]', 'A');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[12]', 'B');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[13]', 'C');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[14]', 'D');
    
    clear_analog_labels{i} = append(clear_analog_labels{i}, ' ', analog_labels{i});
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '.', ' ');
end

% Get force plate corners
corners_force_plates = Rec1_c3d_struct.parameters.FORCE_PLATFORM.CORNERS.DATA;

% ---
% Get the data you actually need
% ---

% Force plates
PILS_meas_mat = get_fp_fm_meas(analog_signals, 3);
cop_plate = get_fp_cop_meas(analog_signals, corners_force_plates, 3);
Ftot_plate = sqrt(PILS_meas_mat(:, 1).^2 + PILS_meas_mat(:, 2).^2 + PILS_meas_mat(:, 3).^2);

% Rod
PILS_ref_mat = get_pole_fm_meas(analog_signals, markers, corners_force_plates, 3);
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


% %% dataRec2
% %  --------
% 
% % ---
% % Read a C3D file and extract some data in a clear format
% % ---
% 
% Rec2_c3d_struct = ezc3dRead(dataRec2);
% 
% % Extract marker data
% markers = Rec2_c3d_struct.data.points; % Of size (3, n_markers, n_frames_markers)
% 
% % Extract analog data
% analog_signals = Rec2_c3d_struct.data.analogs; % Of size (n_frames_analog, n_analog_signals)
% 
% % Construct clear analog labels
% analog_labels = Rec2_c3d_struct.parameters.ANALOG.LABELS.DATA;
% analog_descriptions = Rec2_c3d_struct.parameters.ANALOG.DESCRIPTIONS.DATA;
% clear_analog_labels = analog_descriptions;
% for i = 1:size(analog_labels, 1)
%     if startsWith(clear_analog_labels{i}, 'Analog EMG::Voltage') || startsWith(clear_analog_labels{i}, 'Generic Analog::Electric Potential')
%         clear_analog_labels{i} = '';
%     end
% 
%     clear_analog_labels{i} = replace(clear_analog_labels{i}, 'Kistler Force Plate 2.0.0.0', 'Force plate');
%     clear_analog_labels{i} = replace(clear_analog_labels{i}, '::Raw', '');
%     clear_analog_labels{i} = replace(clear_analog_labels{i}, '[15]', 'A');
%     clear_analog_labels{i} = replace(clear_analog_labels{i}, '[12]', 'B');
%     clear_analog_labels{i} = replace(clear_analog_labels{i}, '[13]', 'C');
%     clear_analog_labels{i} = replace(clear_analog_labels{i}, '[14]', 'D');
% 
%     clear_analog_labels{i} = append(clear_analog_labels{i}, ' ', analog_labels{i});
%     clear_analog_labels{i} = replace(clear_analog_labels{i}, '.', ' ');
% end
% 
% % Get force plate corners
% corners_force_plates = Rec2_c3d_struct.parameters.FORCE_PLATFORM.CORNERS.DATA;
% 
% % ---
% % Get the data you actually need
% % ---
% 
% % Force plates
% PILS_meas_mat = get_fp_fm_meas(analog_signals, 3);
% cop_plate = get_fp_cop_meas(analog_signals, corners_force_plates, 3);
% Ftot_plate = sqrt(PILS_meas_mat(:, 1).^2 + PILS_meas_mat(:, 2).^2 + PILS_meas_mat(:, 3).^2);
% 
% % Rod
% PILS_ref_mat = get_pole_fm_meas(analog_signals, markers, corners_force_plates, 3);
% tip = get_pole_tip_meas(markers);
% Ftot_rod = sqrt(PILS_ref_mat(:, 1).^2 + PILS_ref_mat(:, 2).^2 + PILS_ref_mat(:, 3).^2);
% 
% % % ---
% % % Compare with plots
% % % ---
% % 
% % figure;
% % plot(Ftot_rod, "blue");
% % hold on
% % plot(Ftot_plate, "red");
% % hold off
% % 
% % figure;
% % scatter(tip(:, 1), tip(:, 2), "blue");
% % hold on
% % scatter(cop_plate(:, 1), cop_plate(:, 2), "red");
% % hold off
% 
% % ---
% % Threshold
% % ---
% 
% F_thresh = 100;
% mask = (Ftot_rod >= F_thresh) & (Ftot_plate >= F_thresh);
% n_pts = sum(mask);
% 
% R_2 = zeros(6, n_pts);
% S_2 = zeros(6, n_pts);
% cop_plate_masked = zeros(n_pts, 3);
% tip_masked = zeros(n_pts, 3);
% 
% j = 1;
% for i=1:size(mask, 1)
%     if mask(i)
%         R_2(:, j) = PILS_ref_mat(i, :)';
%         S_2(:, j) = PILS_meas_mat(i, :)';
%         cop_plate_masked(j, :) = cop_plate(i, :);
%         tip_masked(j, :) = tip(i, :);
%         j = j+1;
%     end
% end
% 
% R = [R, R_2];
% S = [S, S_2];
% 
% % ---
% % Compare with plots
% % ---
% 
% figure;
% scatter(tip_masked(:, 1), tip_masked(:, 2), "blue");
% hold on
% scatter(cop_plate_masked(:, 1), cop_plate_masked(:, 2), "red");
% hold off


%% dataRec3
%  --------

% ---
% Read a C3D file and extract some data in a clear format
% ---

Rec3_c3d_struct = ezc3dRead(dataRec3);

% Extract marker data
markers = Rec3_c3d_struct.data.points; % Of size (3, n_markers, n_frames_markers)

% Extract analog data
analog_signals = Rec3_c3d_struct.data.analogs; % Of size (n_frames_analog, n_analog_signals)

% Construct clear analog labels
analog_labels = Rec3_c3d_struct.parameters.ANALOG.LABELS.DATA;
analog_descriptions = Rec3_c3d_struct.parameters.ANALOG.DESCRIPTIONS.DATA;
clear_analog_labels = analog_descriptions;
for i = 1:size(analog_labels, 1)
    if startsWith(clear_analog_labels{i}, 'Analog EMG::Voltage') || startsWith(clear_analog_labels{i}, 'Generic Analog::Electric Potential')
        clear_analog_labels{i} = '';
    end

    clear_analog_labels{i} = replace(clear_analog_labels{i}, 'Kistler Force Plate 2.0.0.0', 'Force plate');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '::Raw', '');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[15]', 'A');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[12]', 'B');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[13]', 'C');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[14]', 'D');
    
    clear_analog_labels{i} = append(clear_analog_labels{i}, ' ', analog_labels{i});
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '.', ' ');
end

% Get force plate corners
corners_force_plates = Rec3_c3d_struct.parameters.FORCE_PLATFORM.CORNERS.DATA;

% ---
% Get the data you actually need
% ---

% Force plates
PILS_meas_mat = get_fp_fm_meas(analog_signals, 3);
cop_plate = get_fp_cop_meas(analog_signals, corners_force_plates, 3);
Ftot_plate = sqrt(PILS_meas_mat(:, 1).^2 + PILS_meas_mat(:, 2).^2 + PILS_meas_mat(:, 3).^2);

% Rod
PILS_ref_mat = get_pole_fm_meas(analog_signals, markers, corners_force_plates, 3);
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

R_3 = zeros(6, n_pts);
S_3 = zeros(6, n_pts);
cop_plate_masked = zeros(n_pts, 3);
tip_masked = zeros(n_pts, 3);

j = 1;
for i=1:size(mask, 1)
    if mask(i)
        R_3(:, j) = PILS_ref_mat(i, :)';
        S_3(:, j) = PILS_meas_mat(i, :)';
        cop_plate_masked(j, :) = cop_plate(i, :);
        tip_masked(j, :) = tip(i, :);
        j = j+1;
    end
end

R = [R, R_3];
S = [S, S_3];

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


% add the following:
% rms cop error (total and components)


%% Calculate az0
%  -------------

% As validation for given az0 = 0.41: is correct!

% a = 0.21;
% b = 0.26;
% 
% fc = 25;
% fs = 1000;
% [but_b, but_a] = butter(2, fc/(fs/2), 'low');
% 
% Fx = -(analog_signals(:, 41) + analog_signals(:, 42));
% Fx = filtfilt(but_b, but_a, Fx);
% Fx = Fx(1:10:size(Fx, 1));
% My = -a * (-analog_signals(:, 45) + analog_signals(:, 46) + analog_signals(:, 47) - analog_signals(:, 48));
% My = filtfilt(but_b, but_a, My);
% My = My(1:10:size(My, 1));
% My_acc = filtfilt(but_b, but_a, analog_signals(:, 17)/1000);
% My_acc = My_acc(1:10:size(My_acc, 1));
% 
% az0 = (My - My_acc)./Fx;
% az0_masked = zeros(n_pts, 1);
% 
% j = 1;
% for i=1:size(mask, 1)
%     if mask(i)
%         az0_masked(j) = az0(i);
%         j = j+1;
%     end
% end

% a = 0.21;
% b = 0.26;
% 
% fc = 25;
% fs = 1000;
% [but_b, but_a] = butter(2, fc/(fs/2), 'low');
% 
% Fy = -(analog_signals(:, 43) + analog_signals(:, 44));
% Fy = filtfilt(but_b, but_a, Fy);
% Fy = Fy(1:10:size(Fy, 1));
% Mx = -b * (analog_signals(:, 45) + analog_signals(:, 46) - analog_signals(:, 47) - analog_signals(:, 48));
% Mx = filtfilt(but_b, but_a, Mx);
% Mx = Mx(1:10:size(Mx, 1));
% Mx_acc = filtfilt(but_b, but_a, analog_signals(:, 16)/1000);
% Mx_acc = Mx_acc(1:10:size(Mx_acc, 1));
% 
% az0 = (Mx_acc - Mx)./Fy;
% az0_masked = zeros(n_pts, 1);
% 
% j = 1;
% for i=1:size(mask, 1)
%     if mask(i)
%         az0_masked(j) = az0(i);
%         j = j+1;
%     end
% end


%% Plot rod f.o.r.
%  ---------------

% side_marker = squeeze(markers(:, 1, :))';
% up_marker = squeeze(markers(:, 2, :))';
% center_marker = squeeze(markers(:, 3, :))';
% 
% x_axis_rod = side_marker - center_marker;
% for i=1:size(x_axis_rod, 1)
%     x_axis_rod(i, :) = x_axis_rod(i, :) / norm(x_axis_rod(i, :));
% end
% 
% z_axis_rod = (up_marker - 4.1*x_axis_rod) - center_marker; % correction to make axis parallel with bar
% for i=1:size(z_axis_rod, 1)
%     z_axis_rod(i, :) = z_axis_rod(i, :) / norm(z_axis_rod(i, :));
% end
% 
% y_axis_rod = zeros(size(x_axis_rod, 1), 3);
% for i=1:size(y_axis_rod)
%     y_axis_rod(i, :) = cross(z_axis_rod(i, :), x_axis_rod(i, :));
% end
% 
% % 1, 2, 3, is for x-, y-, z-axis
% x = zeros(size(x_axis_rod, 1), 3);
% y = zeros(size(x_axis_rod, 1), 3);
% z = zeros(size(x_axis_rod, 1), 3);
% u = zeros(size(x_axis_rod, 1), 3);
% v = zeros(size(x_axis_rod, 1), 3);
% w = zeros(size(x_axis_rod, 1), 3);
% 
% for i = 1:size(x_axis_rod, 1)
%     % Origin is the same for each axis
%     x(i, 1) = center_marker(i, 1);
%     y(i, 1) = center_marker(i, 2);
%     z(i, 1) = center_marker(i, 3);
% 
%     x(i, 2) = center_marker(i, 1);
%     y(i, 2) = center_marker(i, 2);
%     z(i, 2) = center_marker(i, 3);
% 
%     x(i, 3) = center_marker(i, 1);
%     y(i, 3) = center_marker(i, 2);
%     z(i, 3) = center_marker(i, 3);
% 
%     % Direction
%     u(i, 1) = x_axis_rod(i, 1);
%     v(i, 1) = x_axis_rod(i, 2);
%     w(i, 1) = x_axis_rod(i, 3);
% 
%     u(i, 2) = y_axis_rod(i, 1);
%     v(i, 2) = y_axis_rod(i, 2);
%     w(i, 2) = y_axis_rod(i, 3);
% 
%     u(i, 3) = z_axis_rod(i, 1);
%     v(i, 3) = z_axis_rod(i, 2);
%     w(i, 3) = z_axis_rod(i, 3);
% end
% 
% figure
% quiver3(x(:, 1), y(:, 1), z(:, 1), u(:, 1), v(:, 1), w(:, 1), "red")
% hold on
% quiver3(x(:, 2), y(:, 2), z(:, 2), u(:, 2), v(:, 2), w(:, 2), "green")
% quiver3(x(:, 3), y(:, 3), z(:, 3), u(:, 3), v(:, 3), w(:, 3), "blue")
% axis equal
% hold off


%% dataRec4: validation
%  --------------------

% ---
% Read a C3D file and extract some data in a clear format
% ---

Rec4_c3d_struct = ezc3dRead(dataRec4);

% Extract marker data
markers = Rec4_c3d_struct.data.points; % Of size (3, n_markers, n_frames_markers)

% Extract analog data
analog_signals = Rec4_c3d_struct.data.analogs; % Of size (n_frames_analog, n_analog_signals)

% Construct clear analog labels
analog_labels = Rec4_c3d_struct.parameters.ANALOG.LABELS.DATA;
analog_descriptions = Rec4_c3d_struct.parameters.ANALOG.DESCRIPTIONS.DATA;
clear_analog_labels = analog_descriptions;
for i = 1:size(analog_labels, 1)
    if startsWith(clear_analog_labels{i}, 'Analog EMG::Voltage') || startsWith(clear_analog_labels{i}, 'Generic Analog::Electric Potential')
        clear_analog_labels{i} = '';
    end

    clear_analog_labels{i} = replace(clear_analog_labels{i}, 'Kistler Force Plate 2.0.0.0', 'Force plate');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '::Raw', '');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[15]', 'A');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[12]', 'B');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[13]', 'C');
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '[14]', 'D');
    
    clear_analog_labels{i} = append(clear_analog_labels{i}, ' ', analog_labels{i});
    clear_analog_labels{i} = replace(clear_analog_labels{i}, '.', ' ');
end

% Get force plate corners
corners_force_plates = Rec4_c3d_struct.parameters.FORCE_PLATFORM.CORNERS.DATA;

% ---
% Get the data you actually need
% ---

% Force plates
PILS_meas_mat = get_fp_fm_meas(analog_signals, 3);
cop_plate = get_fp_cop_meas(analog_signals, corners_force_plates, 3);
Ftot_plate = sqrt(PILS_meas_mat(:, 1).^2 + PILS_meas_mat(:, 2).^2 + PILS_meas_mat(:, 3).^2);

% Rod
PILS_ref_mat = get_pole_fm_meas(analog_signals, markers, corners_force_plates, 3);
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
cop_fp = get_cop_from_fm(S_val, corners_force_plates, 3);
cop_rod = get_cop_from_fm(R_val, corners_force_plates, 3);
cop_fp_corr = get_cop_from_fm(S_val_corr, corners_force_plates, 3);

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
writematrix(C, strcat(saveFolder, '\C_matrix_FP3D.txt'))

