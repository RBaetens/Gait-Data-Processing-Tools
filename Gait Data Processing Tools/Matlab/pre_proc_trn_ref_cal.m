%% Start with a clean workspace and complete
%  -----------------------------------------

% Clean
close all
clear
clc

% Add toolbox
addpath(".\pre_proc_tools")


%% Declare variables
%  -----------------

% Original data location and file names
dataFolder = "..\..\Gait Classification and References\Data\Subject20\Session1";

dataMVC_TibAnt = dataFolder + "\MVC_TibAnt.csv";
dataMVC_Gast = dataFolder + "\MVC_Gast.csv";
dataMVC_Quad = dataFolder + "\MVC_Quad.csv";
dataMVC_Ham = dataFolder + "\MVC_Ham.csv";

dataAllTraj = dataFolder + "\spatCal_AllTraj.csv";
dataEMG = dataFolder + "\spatCal_EMG.csv";
dataFP = dataFolder + "\spatCal_ForcePlates.csv";
dataBut = dataFolder + "\spatCal_Buttons.csv";
dataMot = dataFolder + "\spatCal.txt";

% Writing location for processed data
syncDataFolder = "..\..\Gait Classification and References\Data\Subject20\Session1_Sync";

% Subject information
BW = 68*9.81; % body weight in Newton
AW = 65; % ankle width in mm
% insole_length_x = 288.0; % size 7, in mm
% insole_length_y = 101.4; % size 7, in mm
% insole_length_x = 274.2; % size 6, in mm
% insole_length_y = 97.5; % size 6, in mm
% insole_length_x = 261.1; % size 5, in mm
% insole_length_y = 93.8; % size 5, in mm
insole_length_x = 248.6; % size 4, in mm
insole_length_y = 90.2; % size 4, in mm

% Trial information
mark_diam = 14; % marker diameter in mm
mark_base = 2; % height of marker base in mm
marker_mode = "full"; % "full" or "lower"

% Sampling information
f_AllTraj = 100;
%f_FP = 1000;
f_FP = 200;
%f_FP = 100;
f_EMG = 1000;
f_But = 1000;
f_Mot = 100;


%% Read data from csv files
%  ------------------------

tabMVC_TibAnt = readmatrix(dataMVC_TibAnt); % Tibialis anterior, gastrocnemius medialis, gastrocnemius lateralis, vastus medialis, vastus lateralis, rectus femoris, semitendinosus, biceps femoris
tabMVC_Gast = readmatrix(dataMVC_Gast);
tabMVC_Quad = readmatrix(dataMVC_Quad);
tabMVC_Ham = readmatrix(dataMVC_Ham);

tabAllTraj = readmatrix(dataAllTraj); % LASI (xyz), RASI (xyz), LPSI (...), RPSI, LTHI, LKNE, LTIB, LANK, LHEE, LTOE, RTHI, RKNE, RTIB, RANK, RHEE, RTOE
tabEMG = readmatrix(dataEMG); % Tibialis anterior, gastrocnemius medialis, gastrocnemius lateralis, vastus medialis, vastus lateralis, rectus femoris, semitendinosus, biceps femoris
tabFP = readmatrix(dataFP); % Plate A: Fx, Fy, Fz, Mx, My, Mz, Cx, Cy, Cz, FX12, FX34, FY14, FY23, FZ1, FZ2, FZ3, FZ4; Plate B: ... ; Plate C: ... ; Plate D: ...
tabBut = readmatrix(dataBut); % Error button, rough segmentation button, fine segmentation button, Moticon sync signal
tabMot = readmatrix(dataMot); % LP1, LP2, LP3, LP4, LP5, LP6, LP7, LP8, LP9, LP10, LP11, LP12, LP13, LP14, LP15, LP16, LACCX, LACCY, LACCZ, LANGX, LANGY, LANGZ, LTF, LCOPX, LCOPY, RP1, RP2, RP3, RP4, RP5, RP6, RP7, RP8, RP9, RP10, RP11, RP12, RP13, RP14, RP15, RP16, RACCX, RACCY, RACCZ, RANGX, RANGY, RANGZ, RTF, RCOPX, RCOPY


%% Check if all vicon data covers the same time segment
%  ----------------------------------------------------

[tabAllTraj, tabBut, tabEMG, tabFP] = check_measurement_lengths(tabAllTraj, tabBut, tabEMG, tabFP, f_AllTraj, f_But, f_EMG, f_FP);


%% Restructure data
%  ----------------

% MVC callibration
% ================

matMVC_TibAnt = tabMVC_TibAnt(4:end, 3:end);
matMVC_Gast = tabMVC_Gast(4:end, 3:end);
matMVC_Quad = tabMVC_Quad(4:end, 3:end);
matMVC_Ham = tabMVC_Ham(4:end, 3:end);


% Spatial callibration
% ====================

% trajectories
if marker_mode == "full"
    matAllTraj = tabAllTraj(4:end, 72:end);
elseif marker_mode == "lower"
    matAllTraj = tabAllTraj(4:end, 3:end);
end

% EMG
matEMG = tabEMG(4:end, 3:10);

% Force plates
% fp1 force xyz, fp1 cop xyz, fp1 moment z, fp2 force xyz, fp2 cop xyz , fp2 moment z, ...
length_tabFP = size(tabFP, 1);
ncols_tabFP = size(tabFP, 2);
if ncols_tabFP == 70
    matFP = zeros(length_tabFP-3, 4*7);
    matFP(:, 1:3) = tabFP(4:end, 3:5);
    matFP(:, 4:6) = tabFP(4:end, 9:11);
    matFP(:, 7) = tabFP(4:end, 8);
    matFP(:, 8:10) = tabFP(4:end, 20:22);
    matFP(:, 11:13) = tabFP(4:end, 26:28);
    matFP(:, 14) = tabFP(4:end, 25);
    matFP(:, 15:17) = tabFP(4:end, 37:39);
    matFP(:, 18:20) = tabFP(4:end, 43:45);
    matFP(:, 21) = tabFP(4:end, 42);
    matFP(:, 22:24) = tabFP(4:end, 54:56);
    matFP(:, 25:27) = tabFP(4:end, 60:62);
    matFP(:, 28) = tabFP(4:end, 59);
elseif ncols_tabFP == 38
    matFP = zeros(length_tabFP-3, 4*7);
    matFP(:, 1:3) = tabFP(4:end, 3:5);
    matFP(:, 4:6) = tabFP(4:end, 9:11);
    matFP(:, 7) = tabFP(4:end, 8);
    matFP(:, 8:10) = tabFP(4:end, 12:14);
    matFP(:, 11:13) = tabFP(4:end, 18:20);
    matFP(:, 14) = tabFP(4:end, 17);
    matFP(:, 15:17) = tabFP(4:end, 21:23);
    matFP(:, 18:20) = tabFP(4:end, 27:29);
    matFP(:, 21) = tabFP(4:end, 26);
    matFP(:, 22:24) = tabFP(4:end, 30:32);
    matFP(:, 25:27) = tabFP(4:end, 36:38);
    matFP(:, 28) = tabFP(4:end, 35);
else
    error('Force plate data has invalid size.')
end

% Buttons
matBut = tabBut(4:end, 3:end);

% Moticon
% left total force, left cop x, left cop y, right total force, right cop x, right cop y
length_tabMot = size(tabMot, 1);
matMotPrim = zeros(length_tabMot-3, 6);
matMotPrim(:, 1:3) = tabMot(4:end, 24:26);
matMotPrim(:, 4:6) = tabMot(4:end, 49:51);

% 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 laccX laccY laccZ langX langY langZ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 raccX raccY raccZ rangX rangY rangZ
matMotSecu = zeros(length_tabMot-3, 44);
matMotSecu(:, 1:22) = tabMot(4:end, 2:23);
matMotSecu(:, 23:44) = tabMot(4:end, 27:48);


%% Make some room in the memory/clean up workspace
%  -----------------------------------------------

clearvars tabMVC_TibAnt tabMVC_Gast tabMVC_Quad tabMVC_Ham
clearvars tabAllTraj tabEMG tabFP tabBut tabMot
clearvars dataMVC_TibAnt dataMVC_Gast dataMVC_Quad dataMVC_Ham
clearvars dataAllTraj dataEMG dataFP dataBut dataMot

%% Fill in isolated missing data
%  -----------------------------

matMVC_TibAnt = simple_impute(matMVC_TibAnt);
matMVC_Gast = simple_impute(matMVC_Gast);
matMVC_Quad = simple_impute(matMVC_Quad);
matMVC_Ham = simple_impute(matMVC_Ham);

matAllTraj = simple_impute(matAllTraj);
matEMG = simple_impute(matEMG);
matFP = simple_impute(matFP);
matBut = simple_impute(matBut);
matMotPrim = simple_impute(matMotPrim);
matMotSecu = simple_impute(matMotSecu);


%% Resample
%  --------

matAllTraj = simple_resample(matAllTraj, f_AllTraj, f_Mot);
matFP = simple_resample(matFP, f_FP, f_Mot);
matBut = simple_resample(matBut, f_But, f_Mot);

% Time axes etc.
n_down = size(matAllTraj, 1);
dt_down = 1/f_Mot;
T_down = dt_down*(n_down-1);
t_down = linspace(0, T_down, n_down);

n_Mot = size(matMotPrim, 1);
dt_Mot = 1/f_Mot;
T_Mot = dt_Mot*(n_Mot-1);
t_Mot = linspace(0, T_Mot, n_Mot);

% Compensation for discretisation
n_cutoff_EMG = int32( ( (n_down*f_Mot/f_AllTraj) - double(int32(floor(n_down*f_Mot/f_AllTraj))) )*f_EMG/f_Mot );
matEMG = matEMG(1:end-n_cutoff_EMG, :);


%% Visualise signals before synchronisation
%  ----------------------------------------

% vertical forces from force plates and insoles
Fv_FP1 = -matFP(:, 3);
Fv_FP2 = -matFP(:, 10);
Fv_FP3 = -matFP(:, 17);
Fv_FP4 = -matFP(:, 24);

Fv_INL = matMotPrim(:, 1);
Fv_INR = matMotPrim(:, 4);

% visualise before further synchronisation
figure
p1 = plot(t_down, Fv_FP1, 'red');
t1 = "FP1";
hold on
p2 = plot(t_down, Fv_FP2, 'green');
t2 = "FP2";
p3 = plot(t_down, Fv_FP3, 'blue');
t3 = "FP3";
p4 = plot(t_down, Fv_FP4, 'yellow');
t4 = "FP4";
p5 = plot(t_Mot, Fv_INL, 'magenta');
t5 = "INL";
p6 = plot(t_Mot, Fv_INR, 'cyan');
t6 = "INR";
legend([p1, p2, p3, p4, p5, p6], [t1, t2, t3, t4, t5, t6]);
title("Unsynchronised signals")
hold off

% Plot error button, rough segmentation button, fine segmentation button and moticon sync signals
figure
subplot(4, 1, 1);
plot(matBut(:, 1));
hold on
subplot(4, 1, 2);
plot(matBut(:, 2));
subplot(4, 1, 3);
plot(matBut(:, 3));
subplot(4, 1, 4);
plot(matBut(:, 4));
hold off

% Manual correction
% % Subject 19, spatCal
% matBut(175:1500, 4) = 2*ones(1326, 1);


%% Crude synchronisation based on synchronisation signal
%  -----------------------------------------------------

% Find index of start and end of block-'wave'
start_stop_indices = block_wave_edges((matBut(:, 4) > 1));
start_ind = start_stop_indices(1, 1);
stop_ind = start_stop_indices(1, 2);

% Cut
matAllTraj = matAllTraj(start_ind:stop_ind, :);
matFP = matFP(start_ind:stop_ind, :);
matBut = matBut(start_ind:stop_ind, :);
matEMG = matEMG(int32((start_ind-1)*(f_EMG/f_Mot))+1:int32(stop_ind*(f_EMG/f_Mot)), :);

t_down = t_down(start_ind:stop_ind);
t_down = t_down-t_down(1);

% Correction if signals are not exactly as long
n_Mot = size(matMotPrim, 1);
n_down = size(matAllTraj, 1);

if n_Mot > n_down
    matMotPrim = matMotPrim(1:n_down, :);
    matMotSecu = matMotSecu(1:n_down, :);
    t_Mot = t_Mot(1:n_down);

elseif n_Mot < n_down
    matAllTraj = matAllTraj(1:n_Mot, :);
    matFP = matFP(1:n_Mot, :);
    matBut = matBut(1:n_Mot, :);
    t_down = t_down(1:n_Mot);

    matEMG = matEMG(1:int32(n_Mot*f_EMG/f_Mot), :);
end


% Visualise crudely synchronised signals
% --------------------------------------

% Vertical forces from force plates and insoles
Fv_FP1 = -matFP(:, 3);
Fv_FP2 = -matFP(:, 10);
Fv_FP3 = -matFP(:, 17);
Fv_FP4 = -matFP(:, 24);

Fv_INL = matMotPrim(:, 1);
Fv_INR = matMotPrim(:, 4);

% Visualise before further synchronisation
figure
p1 = plot(t_down, Fv_FP1, 'red');
t1 = "FP1";
hold on
p2 = plot(t_down, Fv_FP2, 'green');
t2 = "FP2";
p3 = plot(t_down, Fv_FP3, 'blue');
t3 = "FP3";
p4 = plot(t_down, Fv_FP4, 'yellow');
t4 = "FP4";
p5 = plot(t_Mot, Fv_INL, 'magenta');
t5 = "INL";
p6 = plot(t_Mot, Fv_INR, 'cyan');
t6 = "INR";
legend([p1, p2, p3, p4, p5, p6], [t1, t2, t3, t4, t5, t6]);
title("Crudely synchronised signals")
hold off


%% Label: foot over force plate
%  ----------------------------

% Force plate activation threshold
ThLOW = 0.05*BW;

% Calculate which foot is on which FP at any given time
foot_on_fp_cols = foot_on_fp(-[matFP(:, 3), matFP(:, 10), matFP(:, 17), matFP(:, 24)], ... % Fz data
    [matFP(:, 4:6), matFP(:, 11:13), matFP(:, 18:20), matFP(:, 25:27)], ... % CoP data
    matAllTraj(:, 25:27), matAllTraj(:, 28:30), matAllTraj(:, 43:45), matAllTraj(:, 46:48), ... % LHEE, LTOE, RHEE, RTOE
    ThLOW);


% Plot foot on FP signals
figure
p1_L = plot(foot_on_fp_cols(:, 1), 'red');
t1_L = "FP1_L";
hold on
p1_R = plot(foot_on_fp_cols(:, 2), 'Color', 'red', 'LineStyle', '--');
t1_R = "FP1_R";
p2_L = plot(foot_on_fp_cols(:, 3), 'green');
t2_L = "FP2_L";
p2_R = plot(foot_on_fp_cols(:, 4), 'Color', 'green', 'LineStyle', '--');
t2_R = "FP2_R";
p3_L = plot(foot_on_fp_cols(:, 5), 'blue');
t3_L = "FP3_L";
p3_R = plot(foot_on_fp_cols(:, 6), 'Color', 'blue', 'LineStyle', '--');
t3_R = "FP3_R";
p4_L = plot(foot_on_fp_cols(:, 7), 'yellow');
t4_L = "FP4_L";
p4_R = plot(foot_on_fp_cols(:, 8), 'Color', 'yellow', 'LineStyle', '--');
t4_R = "FP4_R";
legend([p1_L, p2_L, p3_L, p4_L, p1_R, p2_R, p3_R, p4_R], [t1_L, t2_L, t3_L, t4_L, t1_R, t2_R, t3_R, t4_R]);
title("Foot on FP")
hold off


%% Fine synchronisation based on correlation
%  -----------------------------------------

% Maximally allowed shift
max_shift = 50;

% Make a wider version of 'foot_on_fp_cols': this expands the region where
% we expect to be able to meaningfully look for correlation
foot_on_fp_cols_wide = widen_blocks(foot_on_fp_cols, max_shift);

% Find globally best shift
n_shifts = 2*max_shift+1;
shifts = linspace(-max_shift, max_shift, n_shifts);
measures = zeros(n_shifts, 1);

best_shift = NaN;
best_measure = 0;
for i=1:n_shifts
    measures(i) = nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP1(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 1));
    measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP1(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 2));
    measures(i) = measures(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP2(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 3));
    measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP2(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 4));
    measures(i) = measures(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP3(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 5));
    measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP3(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 6));
    measures(i) = measures(i) + nansum(Fv_INL(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP4(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 7));
    measures(i) = measures(i) + nansum(Fv_INR(1+max_shift+shifts(i):end-max_shift+shifts(i)).*Fv_FP4(max_shift+1:end-max_shift).*foot_on_fp_cols_wide(max_shift+1:end-max_shift, 8));

    if measures(i) > best_measure
        best_shift = shifts(i);
        best_measure = measures(i);
    end
end

%figure
%plot(shifts, measures)

% Apply globally best time shift
if best_shift > 0
    Fv_INL = Fv_INL(1+best_shift:end);
    Fv_INR = Fv_INR(1+best_shift:end);
    t_Mot = t_Mot(1+best_shift:end)-t_Mot(1+best_shift);
    matMotPrim = matMotPrim(1+best_shift:end, :);
    matMotSecu = matMotSecu(1+best_shift:end, :);

    Fv_FP1 = Fv_FP1(1:end-best_shift);
    Fv_FP2 = Fv_FP2(1:end-best_shift);
    Fv_FP3 = Fv_FP3(1:end-best_shift);
    Fv_FP4 = Fv_FP4(1:end-best_shift);
    t_down = t_down(1:end-best_shift);
    matAllTraj = matAllTraj(1:end-best_shift, :);
    matFP = matFP(1:end-best_shift, :);
    matBut = matBut(1:end-best_shift, :);
    matEMG = matEMG(1:end-int32(best_shift*f_EMG/f_Mot), :);

elseif best_shift < 0
    Fv_INL = Fv_INL(1:end+best_shift);
    Fv_INR = Fv_INR(1:end+best_shift);
    t_Mot = t_Mot(1:end+best_shift);
    matMotPrim = matMotPrim(1:end+best_shift, :);
    matMotSecu = matMotSecu(1:end+best_shift, :);

    Fv_FP1 = Fv_FP1(1-best_shift:end);
    Fv_FP2 = Fv_FP2(1-best_shift:end);
    Fv_FP3 = Fv_FP3(1-best_shift:end);
    Fv_FP4 = Fv_FP4(1-best_shift:end);
    t_down = t_down(1-best_shift:end)-t_down(1-best_shift);
    matAllTraj = matAllTraj(1-best_shift:end, :);
    matFP = matFP(1-best_shift:end, :);
    matBut = matBut(1-best_shift:end, :);
    matEMG = matEMG(1-int32(best_shift*f_EMG/f_Mot):end, :);
end

n_Mot = size(matMotPrim, 1);
n_down = n_Mot;


% Visualise fine synchronisation
% ------------------------------

figure
p1 = plot(t_down, Fv_FP1, 'red');
t1 = "FP1";
hold on
p2 = plot(t_down, Fv_FP2, 'green');
t2 = "FP2";
p3 = plot(t_down, Fv_FP3, 'blue');
t3 = "FP3";
p4 = plot(t_down, Fv_FP4, 'yellow');
t4 = "FP4";
p5 = plot(t_Mot, Fv_INL, 'magenta');
t5 = "INL";
p6 = plot(t_Mot, Fv_INR, 'cyan');
t6 = "INR";
legend([p1, p2, p3, p4, p5, p6], [t1, t2, t3, t4, t5, t6]);
title("Finely synchronised signals")
hold off


%% Make some room in the memory/clean up workspace
%  -----------------------------------------------

clearvars Fv_FP1 Fv_FP2 Fv_FP3 Fv_FP4 Fv_INL Fv_INR
clearvars length_tabFP length_tabMot
clearvars p1 p1_L p1_R p2 p2_L p2_R p3 p3_L p3_R p4 p4_L p4_R p5 p6
clearvars t1 t1_L t1_R t2 t2_L t2_R t3 t3_L t3_R t4 t4_L t4_R t5 t6
clearvars best_measure best_shift max_shift measures n_shifts n_cutoff_EMG shifts
clearvars start_ind stop_ind start_stop_indices


%% Find calibration zones
%  ----------------------

% Make narrower version of foot_on_fp_cols
n_safety_margin = 5;
foot_on_fp_cols_narr = narrow_blocks(foot_on_fp_cols, n_safety_margin);

% Get indices of each block = callibration zone
callibration_zones = [];
for i = 1:size(foot_on_fp_cols_narr, 2)
    callibration_zones = [callibration_zones; block_wave_edges(foot_on_fp_cols_narr(:, i))];
end

% Select one zone for the left foot and one for the right foot
% callibration_zone_left = NaN;
% callibration_zone_right = NaN;
% for i = 1:2:size(foot_on_fp_cols_narr, 2)
%     indices = block_wave_edges(foot_on_fp_cols_narr(:, i));
%     if isnan(callibration_zone_left) && ~(isempty(indices))
%         callibration_zone_left = indices(1, :);
%     end
% end
% for i = 2:2:size(foot_on_fp_cols_narr, 2)
%     indices = block_wave_edges(foot_on_fp_cols_narr(:, i));
%     if isnan(callibration_zone_right) && ~(isempty(indices))
%         callibration_zone_right = indices(1, :);
%     end
% end
% callibration_zones = [callibration_zone_left; callibration_zone_right];



%% Check validity: ocluded markers
%  -------------------------------

% Checks if any markers are ocluded
valid_mask = ones(n_down, 1);
for i = 1:n_down
    isnan_arr = isnan(matAllTraj(i, :));
    if sum(isnan_arr) > 0
        valid_mask(i) = 0;
    end
end

% Checks if any foot markers are ocluded
valid_mask_feet = ones(n_down, 1);
for i = 1:n_down
    isnan_arr = isnan([matAllTraj(i, 22:30), matAllTraj(i, 40:48)]);
    if sum(isnan_arr) > 0
        valid_mask_feet(i) = 0;
    end
end


%% Find ankle --> global and global --> ankle transformations
%  ----------------------------------------------------------

% Returns zero matrix instead of tranformation for occluded markers!!!
transfos = ankle_to_global(matAllTraj, valid_mask_feet, AW, mark_diam, mark_base);


%% Plot ankle frame of references for certain range
%  ------------------------------------------------

% Sort of a validation
plot_indices = [];
for i = 1:size(callibration_zones, 1)
    plot_indices = [plot_indices, callibration_zones(i, 1):callibration_zones(i, 2)];
end
plot_window = size(plot_indices, 1);

% Ankle marker coordinates to plot (black)
ANKX = matAllTraj(plot_indices, 22).*sum(foot_on_fp_cols(plot_indices, 1:2:end), 2) + matAllTraj(plot_indices, 40).*sum(foot_on_fp_cols(plot_indices, 2:2:end), 2);
ANKY = matAllTraj(plot_indices, 23).*sum(foot_on_fp_cols(plot_indices, 1:2:end), 2) + matAllTraj(plot_indices, 41).*sum(foot_on_fp_cols(plot_indices, 2:2:end), 2);
ANKZ = matAllTraj(plot_indices, 24).*sum(foot_on_fp_cols(plot_indices, 1:2:end), 2) + matAllTraj(plot_indices, 42).*sum(foot_on_fp_cols(plot_indices, 2:2:end), 2);

% 1, 2, 3, is for x-, y-, z-axis
x = zeros(plot_window, 3);
y = zeros(plot_window, 3);
z = zeros(plot_window, 3);
u = zeros(plot_window, 3);
v = zeros(plot_window, 3);
w = zeros(plot_window, 3);

for i = plot_indices
    % Check which ankle to pick for each time step
    if sum(foot_on_fp_cols(i, 1:2:end)) >= 1
        Tfgl = transfos.Tfgls(4*(i-1)+1:4*i, :);
    else
        Tfgl = transfos.Tfgrs(4*(i-1)+1:4*i, :);
    end

    % Origin is the same for each axis
    x(i, 1) = Tfgl(1, 4);
    y(i, 1) = Tfgl(2, 4);
    z(i, 1) = Tfgl(3, 4);

    x(i, 2) = Tfgl(1, 4);
    y(i, 2) = Tfgl(2, 4);
    z(i, 2) = Tfgl(3, 4);

    x(i, 3) = Tfgl(1, 4);
    y(i, 3) = Tfgl(2, 4);
    z(i, 3) = Tfgl(3, 4);

    % Direction
    u(i, 1) = Tfgl(1, 1);
    v(i, 1) = Tfgl(2, 1);
    w(i, 1) = Tfgl(3, 1);

    u(i, 2) = Tfgl(1, 2);
    v(i, 2) = Tfgl(2, 2);
    w(i, 2) = Tfgl(3, 2);

    u(i, 3) = Tfgl(1, 3);
    v(i, 3) = Tfgl(2, 3);
    w(i, 3) = Tfgl(3, 3);
end

figure
quiver3(x(:, 1), y(:, 1), z(:, 1), u(:, 1), v(:, 1), w(:, 1), "red")
hold on
quiver3(x(:, 2), y(:, 2), z(:, 2), u(:, 2), v(:, 2), w(:, 2), "green")
quiver3(x(:, 3), y(:, 3), z(:, 3), u(:, 3), v(:, 3), w(:, 3), "blue")
scatter3(ANKX, ANKY, ANKZ, 5, "black", "filled")
axis equal
hold off


%% Get the CoP pointclouds to match to find the ankle --> insole transformation
%  ----------------------------------------------------------------------------

% Construct CoP points in global f.o.r. in homogeneous coords
COP_FP_HOM = [cat(2, matFP(:, 4:6), ones(n_down, 1)), ... % FP1
    cat(2, matFP(:, 11:13), ones(n_down, 1)), ... % FP2
    cat(2, matFP(:, 18:20), ones(n_down, 1)), ... % FP3
    cat(2, matFP(:, 25:27), ones(n_down, 1))]; % FP4

% Construct CoP points in local f.o.r. in homogeneous coords
% With y-axis in direction of foot, z-axis upwards and x-axis to the right
COP_INL = cat(2, matMotPrim(:, 3)*insole_length_y, matMotPrim(:, 2)*insole_length_x);
COP_INL = cat(2, COP_INL, zeros(n_down, 1));
COP_INR = cat(2, -matMotPrim(:, 6)*insole_length_y, matMotPrim(:, 5)*insole_length_x);
COP_INR = cat(2, COP_INR, zeros(n_down, 1));

COP_INL_HOM = cat(2, COP_INL, ones(n_down, 1));
COP_INR_HOM = cat(2, COP_INR, ones(n_down, 1));

% Construct arrays with only CoP points when a foot is on a force plate
% Also store index to have proper transformation at your disposal later
% These are the cop points that will be used for finding the insole to
% ankle f.o.r. transformation
COP_INL_ACT = [];
COP_FPL_ACT = [];

COP_INR_ACT = [];
COP_FPR_ACT = [];

% Consider each callibration zone separately
for j = 1:size(callibration_zones, 1)
    for i = callibration_zones(j, 1):callibration_zones(j, 2)
        for k = 1:(size(foot_on_fp_cols, 2)/2)
            % Check that FP is active and double check that insole CoP is not zero (for start/end of steps)
            if (foot_on_fp_cols(i, k*2-1) == 1) && (sum(COP_INL_HOM(i, :) == [0, 0, 0, 1]) < 4) % Left
                COP_INL_ACT = cat(1, COP_INL_ACT, [COP_INL_HOM(i, :), i]);
                COP_FPL_ACT = cat(1, COP_FPL_ACT, [COP_FP_HOM(i, k*4-3:k*4), i]);
            end
            if (foot_on_fp_cols(i, k*2) == 1) && (sum(COP_INR_HOM(i, :) == [0, 0, 0, 1]) < 4) % Right
                COP_INR_ACT = cat(1, COP_INR_ACT, [COP_INR_HOM(i, :), i]);
                COP_FPR_ACT = cat(1, COP_FPR_ACT, [COP_FP_HOM(i, k*4-3:k*4), i]);
            end
        end
    end
end


%% Make some room in the memory/clean up workspace
%  -----------------------------------------------

clearvars COP_FP_HOM COP_INL COP_INL_HOM COP_INR COP_INR_HOM
clearvars plot_indices plot_window isnan_arr
clearvars u v w x y z ANKX ANKY ANKZ


%% Find the ankle --> insole transformation
%  ----------------------------------------

% Transform "active" CoP points from FP/global to ankle f.o.r.
COP_FPL_ACT_FOOT = zeros(length(COP_FPL_ACT), 4);
for i=1:length(COP_FPL_ACT)
    cop_i = transpose(COP_FPL_ACT(i, 1:4));
    ind_i = COP_FPL_ACT(i, 5);
    COP_FPL_ACT_FOOT(i, :) = transpose(transfos.Tgfls(4*(ind_i-1)+1:4*ind_i, :)*cop_i);
end

COP_FPR_ACT_FOOT = zeros(length(COP_FPR_ACT), 4);
for i=1:length(COP_FPR_ACT)
    cop_i = transpose(COP_FPR_ACT(i, 1:4));
    ind_i = COP_FPR_ACT(i, 5);
    COP_FPR_ACT_FOOT(i, :) = transpose(transfos.Tgfrs(4*(ind_i-1)+1:4*ind_i, :)*cop_i);
end

% Plot some stuff to see what you have to match
% Plot CoP points in ankle f.o.r. vs insole f.o.r.
figure
scatter(COP_FPL_ACT_FOOT(:, 1), COP_FPL_ACT_FOOT(:, 2), "green")
hold on
scatter(COP_INL_ACT(:, 1), COP_INL_ACT(:, 2), "red")
axis equal
hold off

figure
scatter(COP_FPR_ACT_FOOT(:, 1), COP_FPR_ACT_FOOT(:, 2), "green")
hold on
scatter(COP_INR_ACT(:, 1), COP_INR_ACT(:, 2), "red")
axis equal
hold off

% Find insole to ankle f.o.r. transformation (no stretching or rotation around x- and y-axes), least squares optimisation
min_fun_L = @(vars) dist_ptclds((map_from_params(vars(1), vars(2), vars(3), vars(4), 1, 1) * COP_INL_ACT(:, 1:4)')', COP_FPL_ACT_FOOT(:, 1:4));
min_fun_R = @(vars) dist_ptclds((map_from_params(vars(1), vars(2), vars(3), vars(4), 1, 1) * COP_INR_ACT(:, 1:4)')', COP_FPR_ACT_FOOT(:, 1:4));
vars0 = [0, -60, 75, 0, 1, 1];

vars_L = fminsearch(min_fun_L, vars0);
vars_R = fminsearch(min_fun_R, vars0);

Tifl = map_from_params(vars_L(1), vars_L(2), vars_L(3), vars_L(4), 1, 1);
Tfil = trans_mat_inv(Tifl);

Tifr = map_from_params(vars_R(1), vars_R(2), vars_R(3), vars_R(4), 1, 1);
Tfir = trans_mat_inv(Tifr);

% Validation of results: transform CoP as measured by FP to insole f.o.r.
% for callibration zones
COP_FPL_ACT_IN = zeros(size(COP_FPL_ACT_FOOT));
COP_FPR_ACT_IN = zeros(size(COP_FPR_ACT_FOOT));

for i=1:size(COP_FPL_ACT_FOOT, 1)
    COP_FPL_ACT_IN(i, :) = transpose(Tfil*transpose(COP_FPL_ACT_FOOT(i, :)));
end

for i=1:size(COP_FPR_ACT_FOOT, 1)
    COP_FPR_ACT_IN(i, :) = transpose(Tfir*transpose(COP_FPR_ACT_FOOT(i, :)));
end

avg_cop_error_L = dist_ptclds(COP_FPL_ACT_IN(:, 1:2), COP_INL_ACT(:, 1:2))/size(COP_INL_ACT, 1);
avg_cop_error_R = dist_ptclds(COP_FPR_ACT_IN(:, 1:2), COP_INR_ACT(:, 1:2))/size(COP_INR_ACT, 1);

% Plot CoP points in insole f.o.r., from FP and insole data
figure
scatter(COP_FPL_ACT_IN(:, 1), COP_FPL_ACT_IN(:, 2), "green")
hold on
scatter(COP_INL_ACT(:, 1), COP_INL_ACT(:, 2), "red")
axis equal
hold off

figure
scatter(COP_FPR_ACT_IN(:, 1), COP_FPR_ACT_IN(:, 2), "green")
hold on
scatter(COP_INR_ACT(:, 1), COP_INR_ACT(:, 2), "red")
axis equal
hold off


%% Calculate RMS value of MVC EMGs over window of 300 milliseconds
%  ---------------------------------------------------------------

% Sampling stuff
rms_window = 0.3;
n_rms_window = int32(rms_window*f_EMG);

% Variables to store RMS
matMVC_TibAnt_act = zeros(size(matMVC_TibAnt));
matMVC_Gast_act = zeros(size(matMVC_Gast));
matMVC_Quad_act = zeros(size(matMVC_Quad));
matMVC_Ham_act = zeros(size(matMVC_Ham));

% Actually calculate RMS
for j = 1:size(matMVC_TibAnt, 1)-n_rms_window
    for k = 1:size(matMVC_TibAnt, 2)
        matMVC_TibAnt_act(j, k) = sqrt(mean(matMVC_TibAnt(j:j+n_rms_window-1, k).^2));
    end
end
for j = 1:size(matMVC_Gast, 1)-n_rms_window
    for k = 1:size(matMVC_Gast, 2)
        matMVC_Gast_act(j, k) = sqrt(mean(matMVC_Gast(j:j+n_rms_window-1, k).^2));
    end
end
for j = 1:size(matMVC_Quad, 1)-n_rms_window
    for k = 1:size(matMVC_Quad, 2)
        matMVC_Quad_act(j, k) = sqrt(mean(matMVC_Quad(j:j+n_rms_window-1, k).^2));
    end
end
for j = 1:size(matMVC_Ham, 1)-n_rms_window
    for k = 1:size(matMVC_Ham, 2)
        matMVC_Ham_act(j, k) = sqrt(mean(matMVC_Ham(j:j+n_rms_window-1, k).^2));
    end
end


%% Plot EMG signals and their RMS for possible inspection
%  ------------------------------------------------------

% Tibialis anterior
figure
plot(matMVC_TibAnt(:, 1), "blue");
hold on
plot(matMVC_TibAnt_act(:, 1), "red");
hold off

% Gastrocnemius medialis, gastrocnemius lateralis
figure
subplot(2, 1, 1);
plot(matMVC_Gast(:, 2), "blue");
hold on
plot(matMVC_Gast_act(:, 2), "red");
hold off
subplot(2, 1, 2);
plot(matMVC_Gast(:, 3), "blue");
hold on
plot(matMVC_Gast_act(:, 3), "red");
hold off

% Vastus medialis, vastus lateralis, rectus femoris
figure
subplot(3, 1, 1);
plot(matMVC_Quad(:, 4), "blue");
hold on
plot(matMVC_Quad_act(:, 4), "red");
hold off
subplot(3, 1, 2);
plot(matMVC_Quad(:, 5), "blue");
hold on
plot(matMVC_Quad_act(:, 5), "red");
hold off
subplot(3, 1, 3);
plot(matMVC_Quad(:, 6), "blue");
hold on
plot(matMVC_Quad_act(:, 6), "red");
hold off

% Semitendinosus, biceps femoris
figure
subplot(2, 1, 1);
plot(matMVC_Ham(:, 7), "blue");
hold on
plot(matMVC_Ham_act(:, 7), "red");
hold off
subplot(2, 1, 2);
plot(matMVC_Ham(:, 8), "blue");
hold on
plot(matMVC_Ham_act(:, 8), "red");
hold off


%% Define contraction zones
%  ------------------------

% Get "amplitudes"
TibAnt_max = max(matMVC_TibAnt_act(:, 1));
GastMed_max = max(matMVC_Gast_act(:, 2));
GastLat_max = max(matMVC_Gast_act(:, 3));
VastMed_max = max(matMVC_Quad_act(:, 4));
VastLat_max = max(matMVC_Quad_act(:, 5));
RectFem_max = max(matMVC_Quad_act(:, 6));
SemTen_max = max(matMVC_Ham_act(:, 7));
BicFem_max = max(matMVC_Ham_act(:, 8));

% Get masks that are 1 for strong contraction
TibAnt_mask = matMVC_TibAnt_act(:, 1) >= 0.5*TibAnt_max;
GastMed_mask = matMVC_Gast_act(:, 2) >= 0.5*GastMed_max;
GastLat_mask = matMVC_Gast_act(:, 3) >= 0.5*GastLat_max;
VastMed_mask = matMVC_Quad_act(:, 4) >= 0.5*VastMed_max;
VastLat_mask = matMVC_Quad_act(:, 5) >= 0.5*VastLat_max;
RectFem_mask = matMVC_Quad_act(:, 6) >= 0.5*RectFem_max;
SemTen_mask = matMVC_Ham_act(:, 7) >= 0.5*SemTen_max;
BicFem_mask = matMVC_Ham_act(:, 8) >= 0.5*BicFem_max;


%% Calculate mean MVC RMS per muscle
%  ---------------------------------

MVC_norms = zeros(1, 8);
MVC_norms(1) = mean(matMVC_TibAnt_act(:, 1).*TibAnt_mask);
MVC_norms(2) = mean(matMVC_Gast_act(:, 2).*GastMed_mask);
MVC_norms(3) = mean(matMVC_Gast_act(:, 3).*GastLat_mask);
MVC_norms(4) = mean(matMVC_Quad_act(:, 4).*VastMed_mask);
MVC_norms(5) = mean(matMVC_Quad_act(:, 5).*VastLat_mask);
MVC_norms(6) = mean(matMVC_Quad_act(:, 6).*RectFem_mask);
MVC_norms(7) = mean(matMVC_Ham_act(:, 7).*SemTen_mask);
MVC_norms(8) = mean(matMVC_Ham_act(:, 8).*BicFem_mask);


%% Calculate ACT value of MVC EMGs
%  -------------------------------

% Filter stuff
fc = 6;
[but_b, but_a] = butter(2, fc/(f_EMG/2), 'low');

% Variables to store ACT
matMVC_TibAnt_act = zeros(size(matMVC_TibAnt));
matMVC_Gast_act = zeros(size(matMVC_Gast));
matMVC_Quad_act = zeros(size(matMVC_Quad));
matMVC_Ham_act = zeros(size(matMVC_Ham));

% Actually calculate ACT
for k = 1:size(matMVC_TibAnt, 2)
    matMVC_TibAnt_act(:, k) = filtfilt(but_b, but_a, abs(matMVC_TibAnt(:, k)));
end
for k = 1:size(matMVC_Gast, 2)
    matMVC_Gast_act(:, k) = filtfilt(but_b, but_a, abs(matMVC_Gast(:, k)));
end
for k = 1:size(matMVC_Quad, 2)
    matMVC_Quad_act(:, k) = filtfilt(but_b, but_a, abs(matMVC_Quad(:, k)));
end
for k = 1:size(matMVC_Ham, 2)
    matMVC_Ham_act(:, k) = filtfilt(but_b, but_a, abs(matMVC_Ham(:, k)));
end


%% Plot EMG signals and their RMS for possible inspection
%  ------------------------------------------------------

% Tibialis anterior
figure
plot(matMVC_TibAnt(:, 1), "blue");
hold on
plot(matMVC_TibAnt_act(:, 1), "red");
hold off

% Gastrocnemius medialis, gastrocnemius lateralis
figure
subplot(2, 1, 1);
plot(matMVC_Gast(:, 2), "blue");
hold on
plot(matMVC_Gast_act(:, 2), "red");
hold off
subplot(2, 1, 2);
plot(matMVC_Gast(:, 3), "blue");
hold on
plot(matMVC_Gast_act(:, 3), "red");
hold off

% Vastus medialis, vastus lateralis, rectus femoris
figure
subplot(3, 1, 1);
plot(matMVC_Quad(:, 4), "blue");
hold on
plot(matMVC_Quad_act(:, 4), "red");
hold off
subplot(3, 1, 2);
plot(matMVC_Quad(:, 5), "blue");
hold on
plot(matMVC_Quad_act(:, 5), "red");
hold off
subplot(3, 1, 3);
plot(matMVC_Quad(:, 6), "blue");
hold on
plot(matMVC_Quad_act(:, 6), "red");
hold off

% Semitendinosus, biceps femoris
figure
subplot(2, 1, 1);
plot(matMVC_Ham(:, 7), "blue");
hold on
plot(matMVC_Ham_act(:, 7), "red");
hold off
subplot(2, 1, 2);
plot(matMVC_Ham(:, 8), "blue");
hold on
plot(matMVC_Ham_act(:, 8), "red");
hold off


%% Define contraction zones
%  ------------------------

% Get "amplitudes"
TibAnt_max = max(matMVC_TibAnt_act(:, 1));
GastMed_max = max(matMVC_Gast_act(:, 2));
GastLat_max = max(matMVC_Gast_act(:, 3));
VastMed_max = max(matMVC_Quad_act(:, 4));
VastLat_max = max(matMVC_Quad_act(:, 5));
RectFem_max = max(matMVC_Quad_act(:, 6));
SemTen_max = max(matMVC_Ham_act(:, 7));
BicFem_max = max(matMVC_Ham_act(:, 8));

% Get masks that are 1 for strong contraction
TibAnt_mask = matMVC_TibAnt_act(:, 1) >= 0.5*TibAnt_max;
GastMed_mask = matMVC_Gast_act(:, 2) >= 0.5*GastMed_max;
GastLat_mask = matMVC_Gast_act(:, 3) >= 0.5*GastLat_max;
VastMed_mask = matMVC_Quad_act(:, 4) >= 0.5*VastMed_max;
VastLat_mask = matMVC_Quad_act(:, 5) >= 0.5*VastLat_max;
RectFem_mask = matMVC_Quad_act(:, 6) >= 0.5*RectFem_max;
SemTen_mask = matMVC_Ham_act(:, 7) >= 0.5*SemTen_max;
BicFem_mask = matMVC_Ham_act(:, 8) >= 0.5*BicFem_max;


%% Calculate mean MVC ACT per muscle
%  ---------------------------------

ACT_norms = zeros(1, 8);
ACT_norms(1) = mean(matMVC_TibAnt_act(:, 1).*TibAnt_mask);
ACT_norms(2) = mean(matMVC_Gast_act(:, 2).*GastMed_mask);
ACT_norms(3) = mean(matMVC_Gast_act(:, 3).*GastLat_mask);
ACT_norms(4) = mean(matMVC_Quad_act(:, 4).*VastMed_mask);
ACT_norms(5) = mean(matMVC_Quad_act(:, 5).*VastLat_mask);
ACT_norms(6) = mean(matMVC_Quad_act(:, 6).*RectFem_mask);
ACT_norms(7) = mean(matMVC_Ham_act(:, 7).*SemTen_mask);
ACT_norms(8) = mean(matMVC_Ham_act(:, 8).*BicFem_mask);

%% save output
%  -----------

% synchronised measurements
writematrix(matMotPrim, syncDataFolder + "\sync_Moticon_prim.csv")
writematrix(matMotSecu, syncDataFolder + "\sync_Moticon_secu.csv")
writematrix(matFP, syncDataFolder + "\sync_FP.csv")
writematrix(matAllTraj, syncDataFolder + "\sync_AllTraj.csv")
writematrix(matEMG, syncDataFolder + "\sync_EMG.csv")
writematrix(matBut, syncDataFolder + "\sync_But.csv")

% transformations
writematrix(Tifl, syncDataFolder + "\Tifl.csv")
writematrix(Tifr, syncDataFolder + "\Tifr.csv")
writematrix(Tfil, syncDataFolder + "\Tfil.csv")
writematrix(Tfir, syncDataFolder + "\Tfir.csv")

% MVC normalisation factors
writematrix(MVC_norms, syncDataFolder + "\MVC_norms.csv")
writematrix(ACT_norms, syncDataFolder + "\ACT_norms.csv")



% The end! %