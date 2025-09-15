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
n_recordings = 5;

dataAllTraj = {};
dataEMG = {};
dataFP = {};
dataBut = {};
dataMot = {};
dataModOut = {};
dataGloAng = {};

for i = 1:n_recordings
    dataAllTraj = [dataAllTraj; {dataFolder + "\Recording" + string(i) + "_AllTraj.csv"}];
    dataEMG = [dataEMG; {dataFolder + "\Recording" + string(i) + "_EMG.csv"}];
    dataFP = [dataFP; {dataFolder + "\Recording" + string(i) + "_ForcePlates.csv"}];
    dataBut = [dataBut; {dataFolder + "\Recording" + string(i) + "_Buttons.csv"}];
    dataMot = [dataMot; {dataFolder + "\Recording" + string(i) + ".txt"}];
    dataModOut = [dataModOut; {dataFolder + "\Recording" + string(i) + "_ModelOutputs.csv"}];
    dataGloAng = [dataGloAng; {dataFolder + "\Recording" + string(i) + "_GlobalAngles.csv"}];
end

% Writing location for processed data
syncDataFolder = "..\..\Gait Classification and References\Data\Subject20\Session1_Sync";
procDataFolderTrn = "..\..\Gait Classification and References\Data\Subject20\Session1_Proc_Trn";
procDataFolderRef = "..\..\Gait Classification and References\Data\Subject20\Session1_Proc_Ref";

% Subject information
subject_tag = "Subj20";
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
trials = get_trial_names();

% Sampling information
f_Mot = 100;

% Define buffer length
T_buf = 1;
n_buf = int32(T_buf*f_Mot);

% Define headers
colnamesAllTraj = get_colnamesAllTraj();
colnamesBut = get_colnamesBut();
colnamesEMG = get_colnamesEMG();
colnamesFP = get_colnamesFP();
colnamesModOut = get_colnamesModOut(marker_mode);
colnamesGloAng = get_colnamesGloAng(marker_mode);
colnamesMot = get_colnamesMot();
colnamesGRF = get_colnamesGRF();
colnamesMasks = get_colnamesMasks();
colnamesFOF = get_colnamesFOF();
colnamesFootMarkers = get_colnamesFootMarkers();

% Read Tifl, Tifr, Tfil, Tfir and MVC_norms
Tifl = readmatrix(syncDataFolder + "\Tifl.csv");
Tifr = readmatrix(syncDataFolder + "\Tifr.csv");
Tfil = readmatrix(syncDataFolder + "\Tfil.csv");
Tfir = readmatrix(syncDataFolder + "\Tfir.csv");
MVC_norms = readmatrix(syncDataFolder + "\MVC_norms.csv");


%% Loop through all recordings and process
%  ---------------------------------------

for i_rec = 1:n_recordings
    % Read data from csv files
    % ------------------------
    
    % Actual data
    tabAllTraj = readmatrix(dataAllTraj{i_rec, :});
    tabEMG = readmatrix(dataEMG{i_rec, :});
    tabFP = readmatrix(dataFP{i_rec, :});
    tabBut = readmatrix(dataBut{i_rec, :});
    tabMot = readmatrix(dataMot{i_rec, :});
    tabModOut = readmatrix(dataModOut{i_rec, :});
    tabGloAng = readmatrix(dataGloAng{i_rec, :});

    % Sampling frequencies
    f_AllTraj = get_sampling_frequency(char(dataAllTraj{i_rec, :}));
    f_EMG = get_sampling_frequency(char(dataEMG{i_rec, :}));
    f_FP = get_sampling_frequency(char(dataFP{i_rec, :}));
    f_But = get_sampling_frequency(char(dataBut{i_rec, :}));
    f_ModOut = get_sampling_frequency(char(dataModOut{i_rec, :}));
    f_GloAng = get_sampling_frequency(char(dataGloAng{i_rec, :}));

    % % Manual correction for Subject6
    % if i_rec == 1
    %     tabBut = [tabBut; zeros(9, 6)];
    %     tabEMG = [tabEMG; zeros(9, 10)];
    % end
    % % Manual correction for Subject13
    % if i_rec == 5
    %     tabBut = [tabBut; zeros(9, 6)];
    %     tabEMG = [tabEMG; zeros(9, 10)];
    % end
    % % Manual correction for Subject14
    % if i_rec == 1
    %     tabBut = [tabBut; zeros(9, 6)];
    %     tabEMG = [tabEMG; zeros(9, 10)];
    % end

    % Check if all vicon data covers the same time segment
    [tabAllTraj, tabBut, tabEMG, tabFP] = check_measurement_lengths(tabAllTraj, tabBut, tabEMG, tabFP, f_AllTraj, f_But, f_EMG, f_FP);


    % Restructure data
    % ----------------
    
    % Trajectories
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
    
    % Model outputs
    matModOut = tabModOut(4:end, 3:end);

    % Global angles
    matGloAng = tabGloAng(4:end, 3:end);
    
    % Make some room in the memory/clean up workspace
    clearvars tabAllTraj tabEMG tabFP tabBut tabMot tabModOut tabGloAng
    
    
    % Fill in isolated missing data
    % -----------------------------
    
    matAllTraj = simple_impute(matAllTraj);
    matEMG = simple_impute(matEMG);
    matFP = simple_impute(matFP);
    matBut = simple_impute(matBut);
    matMotPrim = simple_impute(matMotPrim);
    matMotSecu = simple_impute(matMotSecu);
    matModOut = simple_impute(matModOut);
    matGloAng = simple_impute(matGloAng);
    
    
    % Resample
    % --------
    
    matAllTraj = simple_resample(matAllTraj, f_AllTraj, f_Mot);
    matFP = simple_resample(matFP, f_FP, f_Mot);
    matBut = simple_resample(matBut, f_But, f_Mot);
    matModOut = simple_resample(matModOut, f_ModOut, f_Mot);
    matGloAng= simple_resample(matGloAng, f_GloAng, f_Mot);
    
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
    
    % Vertical forces from force plates and insoles
    Fv_FP1 = -matFP(:, 3);
    Fv_FP2 = -matFP(:, 10);
    Fv_FP3 = -matFP(:, 17);
    Fv_FP4 = -matFP(:, 24);
    
    Fv_INL = matMotPrim(:, 1);
    Fv_INR = matMotPrim(:, 4);
    
    % Visualise
    unsync_figure = visualize_vertical_forces("Unsynchronized signals", ...
        Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, Fv_INL, Fv_INR, t_down, t_Mot);
    
    % % Optionally: manually correct faults in button presses
    % % Subject1 Session2
    % if i_rec == 2
    %     matBut(12000:12049, 2) = 2*ones(50, 1);
    %     matBut(13800:13849, 2) = 2*ones(50, 1);
    %     matBut(49000:49999, 3) = zeros(1000, 1);
    %     matBut(94500:94549, 3) = 2*ones(50, 1);
    % end
    % % Subject1 Session3
    % if i_rec == 1
    %     matBut(39100:39119, 3) = 2*ones(20, 1);
    % end
    % % Subject2 Session1
    % if i_rec == 3
    %     matBut(42650:42799, 2) = zeros(150, 1);
    %     matBut(42650:42799, 3) = zeros(150, 1);
    %     matBut(46950:47149, 2) = zeros(200, 1);
    %     matBut(46950:47149, 3) = zeros(200, 1);
    %     matBut(55300:55499, 2) = zeros(200, 1);
    % end
    % if i_rec == 4
    %     matBut(11800:12049, 2) = zeros(250, 1);
    %     matBut(11800:12049, 3) = zeros(250, 1);
    %     matBut(18200:18249, 3) = 2*ones(50, 1);
    %     matBut(27400:27599, 2) = zeros(200, 1);
    %     matBut(27400:27599, 3) = zeros(200, 1);
    %     matBut(58500:58749, 2) = zeros(250, 1);
    %     matBut(58500:58749, 3) = zeros(250, 1);
    %     matBut(65200:65499, 2) = zeros(300, 1);
    %     matBut(65200:65499, 3) = zeros(300, 1);
    %     matBut(74000:74149, 2) = zeros(150, 1);
    %     matBut(76000:76199, 2) = zeros(200, 1);
    %     matBut(79800:80599, 2) = zeros(800, 1);
    %     matBut(79800:80199, 3) = zeros(400, 1);
    % end
    % % Subject3 Session1
    % if i_rec == 2
    %     matBut(5200:5599, 3) = zeros(400, 1);
    % end
    % % Subject4 Session1
    % if i_rec == 1
    %     matBut(24000:24399, 3) = zeros(400, 1);
    % end
    % if i_rec == 2
    %     matBut(19200:19299, 3) = 2*ones(100, 1);
    % end
    % if i_rec == 3
    %     matBut(9600:9799, 3) = zeros(200, 1);
    % end
    % % Subject6 Session1
    % if i_rec == 1
    %     matBut(43000:43099, 1) = 2*ones(100, 1);
    %     matBut(44700:45099, 1) = zeros(400, 1);
    %     matBut(134000:134199, 1) = 2*ones(200, 1);
    % end
    % % Subject7 Session1
    % if i_rec == 3
    %     matBut(32300:32499, 2) = zeros(200, 1);
    % end
    % if i_rec == 4
    %     matBut(68100:68499, 3) = zeros(400, 1);
    % end
    % % Subject8 Session1
    % if i_rec == 3
    %     matBut(10600:10699, 3) = 2*ones(100, 1);
    %     matBut(10900:11099, 3) = zeros(200, 1);
    % end
    % % Subject9 Session1
    % if i_rec == 1
    %     matBut(80800:80999, 2) = zeros(200, 1);
    %     matBut(171100:171499, 3) = zeros(400, 1);
    % end
    % if i_rec == 2
    %     matBut(54100:54299, 3) = zeros(200, 1);
    %     matBut(116400:116599, 2) = zeros(200, 1);
    %     matBut(116400:116599, 3) = zeros(200, 1);
    %     matBut(120200:121099, 2) = zeros(900, 1);
    %     matBut(120200:121099, 3) = zeros(900, 1);
    %     matBut(144800:144999, 2) = zeros(200, 1);
    %     matBut(144800:144999, 3) = zeros(200, 1);
    % end
    % % Subject10 Session1
    % if i_rec == 2
    %     matBut(13500:13899, 3) = zeros(400, 1);
    % end
    % % Subject13 Session1
    % if i_rec == 2
    %     matBut(2400:2699, 2) = zeros(300, 1);
    % end
    % if i_rec == 4
    %     matBut(35800:36799, 3) = zeros(1000, 1);
    % end
    % % Subject15 Session1
    % if i_rec == 2
    %     matBut(2500:2549, 3) = 2*ones(50, 1);
    %     matBut(65000:65299, 3) = zeros(300, 1);
    % end
    % % Subject16 Session1
    % if i_rec == 1
    %     matBut(15000:15999, 3) = zeros(1000, 1);
    %     matBut(45500:45599, 3) = 2*ones(100, 1);
    %     matBut(63000:63199, 3) = zeros(200, 1);
    %     matBut(143200:143249, 3) = 2*ones(50, 1);
    % end
    % if i_rec == 3
    %     matBut(46200:46249, 3) = 2*ones(50, 1);
    % end
    % % Subject 18
    % if i_rec == 1
    %     matBut(62500:62549, 3) = 2*ones(50, 1);
    % end
    % if i_rec == 5
    %     matBut(18000:18399, 2) = zeros(400, 1);
    %     matBut(23000:23499, 2) = zeros(500, 1);
    % end
    % % Subject 19
    % if i_rec == 1
    %     matBut(18500:19000, 2) = zeros(501, 1);
    %     matBut(32000:32500, 2) = zeros(501, 1);
    %     matBut(42000:42500, 2) = zeros(501, 1);
    %     matBut(61500:62000, 2) = zeros(501, 1);
    %     matBut(84000:84500, 2) = zeros(501, 1);
    %     matBut(101500:102000, 2) = zeros(501, 1);
    %     matBut(127000:127500, 2) = zeros(501, 1);
    % end
    % if i_rec == 2
    %     matBut(6500:7000, 2) = zeros(501, 1);
    %     matBut(8700:9200, 2) = zeros(501, 1);
    %     matBut(39200:39700, 2) = zeros(501, 1);
    %     matBut(68000:68500, 2) = zeros(501, 1);
    %     matBut(87000:87500, 2) = zeros(501, 1);
    %     matBut(95500:96000, 2) = zeros(501, 1);
    %     matBut(98000:98500, 2) = zeros(501, 1);
    %     matBut(102400:102900, 2) = zeros(501, 1);
    % end
    % if i_rec == 3
    %     matBut(7500:8000, 2) = zeros(501, 1);
    %     matBut(9500:10000, 2) = zeros(501, 1);
    %     matBut(19000:19500, 2) = zeros(501, 1);
    %     matBut(24500:25000, 2) = zeros(501, 1);
    %     matBut(28800:29300, 2) = zeros(501, 1);
    %     matBut(34500:35000, 2) = zeros(501, 1);
    %     matBut(37500:38000, 2) = zeros(501, 1);
    %     matBut(41500:42000, 2) = zeros(501, 1);
    %     matBut(44000:44500, 2) = zeros(501, 1);
    %     matBut(49000:49500, 2) = zeros(501, 1);
    %     matBut(51800:52300, 2) = zeros(501, 1);
    %     matBut(54500:55000, 2) = zeros(501, 1);
    %     matBut(57000:57500, 2) = zeros(501, 1);
    % end
    % % Subject 20
    % if i_rec == 2
    %     matBut(19400:19900, 2) = zeros(501, 1);
    %     matBut(25000:25500, 2) = zeros(501, 1);
    %     matBut(58000:58500, 2) = zeros(501, 1);
    %     matBut(76000:76500, 2) = zeros(501, 1);
    %     matBut(80000:80500, 2) = zeros(501, 1);
    %     matBut(86000:87000, 2) = zeros(1001, 1);
    %     matBut(94500:95000, 2) = zeros(501, 1);
    %     matBut(100000:100200, 2) = zeros(201, 1);
    % end
    % if i_rec == 5
    %     matBut(10000:10500, 2) = zeros(501, 1);
    %     matBut(16000:16500, 2) = zeros(501, 1);
    %     matBut(24000:24500, 2) = zeros(501, 1);
    % end

    % Plot error button, rough segmentation button, fine segmentation button and moticon sync signals
    buttons_figure = figure;
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
    
    
    % Crude synchronisation based on synchronisation signal
    % -----------------------------------------------------
    
    % Find index of start and end of block-'wave'
    start_stop_indices = block_wave_edges((matBut(:, 4) > 1));
    start_ind = start_stop_indices(1, 1);
    stop_ind = start_stop_indices(1, 2);
    
    % Cut
    matAllTraj = matAllTraj(start_ind:stop_ind, :);
    matFP = matFP(start_ind:stop_ind, :);
    matBut = matBut(start_ind:stop_ind, :);
    matModOut = matModOut(start_ind:stop_ind, :);
    matGloAng = matGloAng(start_ind:stop_ind, :);
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
        matModOut = matModOut(1:n_Mot, :);
        matGloAng = matGloAng(1:n_Mot, :);
        t_down = t_down(1:n_Mot);
    
        matEMG = matEMG(1:int32(n_Mot*f_EMG/f_Mot), :);
    end
    
    % Vertical forces from force plates and insoles
    Fv_FP1 = -matFP(:, 3);
    Fv_FP2 = -matFP(:, 10);
    Fv_FP3 = -matFP(:, 17);
    Fv_FP4 = -matFP(:, 24);
    
    Fv_INL = matMotPrim(:, 1);
    Fv_INR = matMotPrim(:, 4);
    
    % Visualise
    crude_sync_figure = visualize_vertical_forces("Crudely synchronized signals", ...
        Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, Fv_INL, Fv_INR, t_down, t_Mot);

    % % Manual correction for Subject11
    % if i_rec == 1
    %     Fv_FP3(3500:10999) = zeros(7500, 1);
    %     Fv_FP3(14000:16999) = zeros(3000, 1);
    % end
    % % Manual correction for Subject13
    % if i_rec == 1
    %     Fv_FP3(90000:94749) = zeros(4750, 1);
    %     Fv_FP3(122000:123999) = zeros(2000, 1);
    %     Fv_FP2(90000:94749) = zeros(4750, 1);
    %     Fv_FP2(122000:123999) = zeros(2000, 1);
    %     Fv_FP3(69000:69999) = zeros(1000, 1);
    %     Fv_FP3(75500:76999) = zeros(1500, 1);
    % end
    % if i_rec == 3
    %     Fv_FP3(19500:20499) = zeros(1000, 1);
    % end
    % % Manual correction for Subject14
    % if i_rec == 1
    %     Fv_FP4(17000:18999) = zeros(2000, 1);
    %     Fv_FP4(21500:22999) = zeros(1500, 1);
    %     Fv_FP4(76000:76999) = zeros(1000, 1);
    %     Fv_FP4(117000:119399) = zeros(2400, 1);
    %     Fv_FP4(143000:145999) = zeros(3000, 1);
    %     Fv_FP4(156000:159999) = zeros(4000, 1);
    % end
    % if i_rec == 3
    %     Fv_FP4(12800:13399) = zeros(600, 1);
    %     Fv_FP4(42600:43199) = zeros(600, 1);
    %     Fv_FP4(50000:50499) = zeros(500, 1);
    %     Fv_FP4(52500:52999) = zeros(500, 1);
    % end
    % if i_rec == 4
    %     Fv_FP4(40000:42999) = zeros(3000, 1);
    %     Fv_FP4(62500:63499) = zeros(1000, 1);
    % end
    % % Manual correction for Subject15
    % if i_rec == 2
    %     Fv_FP1(137500:138499) = zeros(1000, 1);
    %     Fv_FP1(139500:140999) = zeros(1500, 1);
    %     Fv_FP1(146200:146999) = zeros(800, 1);
    %     Fv_FP1(148800:149399) = zeros(600, 1);
    %     Fv_FP1(152000:152999) = zeros(1000, 1);
    % end
    % % Manual correction for Subject16
    % if i_rec == 1
    %     Fv_FP3(60400:60699) = zeros(300, 1);
    % end
    % if i_rec == 3
    %     Fv_FP3(47200:47499) = zeros(300, 1);
    % end
    % if i_rec == 5
    %     Fv_FP3(2200:2799) = zeros(600, 1);
    %     Fv_FP4(2200:2799) = zeros(600, 1);
    %     Fv_FP3(4000:4599) = zeros(600, 1);
    %     Fv_FP4(4000:4599) = zeros(600, 1);
    %     Fv_FP3(5900:6499) = zeros(600, 1);
    %     Fv_FP4(5900:6499) = zeros(600, 1);
    %     Fv_FP3(8000:8599) = zeros(600, 1);
    %     Fv_FP4(8000:8599) = zeros(600, 1);
    %     Fv_FP3(10000:11499) = zeros(1500, 1);
    %     Fv_FP4(10000:11499) = zeros(1500, 1);
    % end
    % % Subject 17
    % if i_rec == 1
    %     Fv_FP3(1700:2099) = zeros(400, 1);
    %     Fv_FP3(15800:16099) = zeros(300, 1);
    %     Fv_FP3(18200:18399) = zeros(200, 1);
    %     Fv_FP3(22800:23099) = zeros(300, 1);
    %     Fv_FP3(40000:40299) = zeros(300, 1);
    %     Fv_FP3(43300:43499) = zeros(200, 1);
    %     Fv_FP3(167000:179999) = zeros(13000, 1);
    % end
    % % Subject 18
    % if i_rec == 1
    %     Fv_FP3(69000:82999) = zeros(14000, 1);
    % end
    % if i_rec == 3
    %     Fv_FP3(11500:13999) = zeros(2500, 1);
    % end
    % if i_rec == 4
    %     Fv_FP3(16000:20499) = zeros(4500, 1);
    %     Fv_FP4(16000:20499) = zeros(4500, 1);
    % end
    % % Subject 19
    % if i_rec == 3
    %     Fv_FP2(14000:18000) = zeros(4001, 1);
    % end
    % if i_rec == 5
    %     Fv_FP2(1:3500) = zeros(3500, 1);
    %     Fv_FP2(5000:5500) = zeros(501, 1);
    %     Fv_FP2(6500:7000) = zeros(501, 1);
    %     Fv_FP2(9000:9500) = zeros(501, 1);
    %     Fv_FP2(12000:12500) = zeros(501, 1);
    %     Fv_FP3(1:3500) = zeros(3500, 1);
    %     Fv_FP3(5000:5500) = zeros(501, 1);
    %     Fv_FP3(6500:7000) = zeros(501, 1);
    %     Fv_FP3(9000:9500) = zeros(501, 1);
    %     Fv_FP3(12000:12500) = zeros(501, 1);
    % end
    % % Subject 20
    % if i_rec == 1
    %     Fv_FP2(16000:28000) = zeros(12001, 1);
    % end
    % if i_rec == 2
    %     Fv_FP4(67000:74000) = zeros(7001, 1);
    % end

    
    % Label: foot over force plate
    % ----------------------------
    
    % Force plate activation threshold
    ThLOW = 0.05*BW;
    
    % Calculate which foot is on which FP at any given time
    foot_on_fp_cols = foot_on_fp(-[matFP(:, 3), matFP(:, 10), matFP(:, 17), matFP(:, 24)], ... % Fz data
        [matFP(:, 4:6), matFP(:, 11:13), matFP(:, 18:20), matFP(:, 25:27)], ... % CoP data
        matAllTraj(:, 25:27), matAllTraj(:, 28:30), matAllTraj(:, 43:45), matAllTraj(:, 46:48), ... % LHEE, LTOE, RHEE, RTOE
        ThLOW);
    
    % Plot foot on FP signals
    fof_figure = figure;
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
    
    
    % Slightly finer synchronisation based on correlation
    % ---------------------------------------------------
    
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
        matModOut = matModOut(1:end-best_shift, :);
        matGloAng = matGloAng(1:end-best_shift, :);
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
        matModOut = matModOut(1-best_shift:end, :);
        matGloAng = matGloAng(1-best_shift:end, :);
        matEMG = matEMG(1-int32(best_shift*f_EMG/f_Mot):end, :);
    end
    
    n_Mot = size(matMotPrim, 1);
    n_down = n_Mot;
    
    % Visualise
    finer_sync_figure = visualize_vertical_forces("Finely synchronized signals", ...
        Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, Fv_INL, Fv_INR, t_down, t_Mot);
    
    
    % Piecewise synchronisation per insole
    % ------------------------------------
    
    % Find presses in rough segmentation button signal
    ON_OFF_SEGM_ROUGH = block_wave_edges((matBut(:, 2) > 1));
    n_trials = size(ON_OFF_SEGM_ROUGH, 1)-1;
    best_shifts = zeros(n_trials, 2);
    
    % Look at each individual trial and find best shift for that trial
    for i = 1:n_trials
        start_ind = ON_OFF_SEGM_ROUGH(i, 1);
        stop_ind = ON_OFF_SEGM_ROUGH(i+1, 2);
        max_shift = 50;
        
        % The best shift is NaN if there is no activity on the force plates
        % It is also NaN if both feet were on the same force plate a lot
        [best_shift_L, best_shift_R] = find_best_shifts(Fv_FP1(start_ind:stop_ind), Fv_FP2(start_ind:stop_ind), Fv_FP3(start_ind:stop_ind), ...
            Fv_FP4(start_ind:stop_ind), Fv_INL(start_ind:stop_ind), Fv_INR(start_ind:stop_ind), ...
            matFP(start_ind:stop_ind, :), matAllTraj(start_ind:stop_ind, :), ThLOW, max_shift);
        best_shifts(i, :) = [best_shift_L, best_shift_R];
    end
    
    % Impute NaN values as accurately as possible
    best_shifts(:, 1) = shift_impute(best_shifts(:, 1));
    best_shifts(:, 2) = shift_impute(best_shifts(:, 2));
    
    % Split matMot in left and right part for separate synchronisation
    matMotPrim_L = matMotPrim(:, 1:3);
    matMotPrim_R = matMotPrim(:, 4:6);
    matMotSecu_L = matMotSecu(:, 1:22);
    matMotSecu_R = matMotSecu(:, 23:44);
    
    % Delete some rows of data to apply the shift
    for i = 1:n_trials
        ON_OFF_SEGM_ROUGH = block_wave_edges((matBut(:, 2) > 1));
        start_ind = ON_OFF_SEGM_ROUGH(i, 1);
        stop_ind = ON_OFF_SEGM_ROUGH(i+1, 2);
        best_shift_L = best_shifts(i, 1);
        best_shift_R = best_shifts(i, 2);
    
        if (sign(best_shift_L) == sign(best_shift_R)) && (best_shift_L > 0) && (best_shift_L >= best_shift_R)
            Fv_INL(start_ind:start_ind+best_shift_L-1) = [];
            Fv_INR(start_ind:start_ind+best_shift_R-1) = [];
            Fv_INR(stop_ind-(best_shift_L-best_shift_R)+1:stop_ind) = [];
            matMotPrim_L(start_ind:start_ind+best_shift_L-1, :) = [];
            matMotPrim_R(start_ind:start_ind+best_shift_R-1, :) = [];
            matMotPrim_R(stop_ind-(best_shift_L-best_shift_R)+1:stop_ind, :) = [];
            matMotSecu_L(start_ind:start_ind+best_shift_L-1, :) = [];
            matMotSecu_R(start_ind:start_ind+best_shift_R-1, :) = [];
            matMotSecu_R(stop_ind-(best_shift_L-best_shift_R)+1:stop_ind, :) = [];
            
            Fv_FP1(stop_ind-best_shift_L+1:stop_ind) = [];
            Fv_FP2(stop_ind-best_shift_L+1:stop_ind) = [];
            Fv_FP3(stop_ind-best_shift_L+1:stop_ind) = [];
            Fv_FP4(stop_ind-best_shift_L+1:stop_ind) = [];
            matAllTraj(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matFP(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matBut(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matModOut(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matGloAng(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matEMG(stop_ind - int32(best_shift_L*f_EMG/f_Mot) + 1:stop_ind, :) = [];
            
        elseif (sign(best_shift_L) == sign(best_shift_R)) && (best_shift_L > 0) && (best_shift_L < best_shift_R)
            Fv_INL(start_ind:start_ind+best_shift_L-1) = [];
            Fv_INL(stop_ind-(best_shift_R-best_shift_L)+1:stop_ind) = [];
            Fv_INR(start_ind:start_ind+best_shift_R-1) = [];
            matMotPrim_L(start_ind:start_ind+best_shift_L-1, :) = [];
            matMotPrim_L(stop_ind-(best_shift_R-best_shift_L)+1:stop_ind, :) = [];
            matMotPrim_R(start_ind:start_ind+best_shift_R-1, :) = [];
            matMotSecu_L(start_ind:start_ind+best_shift_L-1, :) = [];
            matMotSecu_L(stop_ind-(best_shift_R-best_shift_L)+1:stop_ind, :) = [];
            matMotSecu_R(start_ind:start_ind+best_shift_R-1, :) = [];
            
            Fv_FP1(stop_ind-best_shift_R+1:stop_ind) = [];
            Fv_FP2(stop_ind-best_shift_R+1:stop_ind) = [];
            Fv_FP3(stop_ind-best_shift_R+1:stop_ind) = [];
            Fv_FP4(stop_ind-best_shift_R+1:stop_ind) = [];
            matAllTraj(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matFP(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matBut(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matModOut(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matGloAng(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matEMG(stop_ind - int32(best_shift_R*f_EMG/f_Mot) + 1:stop_ind, :) = [];
            
        elseif (sign(best_shift_L) == sign(best_shift_R)) && (best_shift_L < 0) && (best_shift_L < best_shift_R)
            Fv_INL(stop_ind+best_shift_L+1:stop_ind) = [];
            Fv_INR(stop_ind+best_shift_R+1:stop_ind) = [];
            Fv_INR(start_ind+(best_shift_L-best_shift_R)+1:start_ind) = [];
            matMotPrim_L(stop_ind+best_shift_L+1:stop_ind, :) = [];
            matMotPrim_R(stop_ind+best_shift_R+1:stop_ind, :) = [];
            matMotPrim_R(start_ind+(best_shift_L-best_shift_R)+1:start_ind, :) = [];
            matMotSecu_L(stop_ind+best_shift_L+1:stop_ind, :) = [];
            matMotSecu_R(stop_ind+best_shift_R+1:stop_ind, :) = [];
            matMotSecu_R(start_ind+(best_shift_L-best_shift_R)+1:start_ind, :) = [];
            
            Fv_FP1(start_ind+best_shift_L+1:start_ind) = [];
            Fv_FP2(start_ind+best_shift_L+1:start_ind) = [];
            Fv_FP3(start_ind+best_shift_L+1:start_ind) = [];
            Fv_FP4(start_ind+best_shift_L+1:start_ind) = [];
            matAllTraj(start_ind+best_shift_L+1:start_ind, :) = [];
            matFP(start_ind+best_shift_L+1:start_ind, :) = [];
            matBut(start_ind+best_shift_L+1:start_ind, :) = [];
            matModOut(start_ind+best_shift_L+1:start_ind, :) = [];
            matGloAng(start_ind+best_shift_L+1:start_ind, :) = [];
            matEMG(start_ind + int32(best_shift_L*f_EMG/f_Mot) + 1:start_ind, :) = [];
            
        elseif (sign(best_shift_L) == sign(best_shift_R)) && (best_shift_L < 0) && (best_shift_L >= best_shift_R)
            Fv_INL(stop_ind+best_shift_L+1:stop_ind) = [];
            Fv_INL(start_ind+(best_shift_R-best_shift_L)+1:start_ind) = [];
            Fv_INR(stop_ind+best_shift_R+1:stop_ind) = [];
            matMotPrim_L(stop_ind+best_shift_L+1:stop_ind, :) = [];
            matMotPrim_L(start_ind+(best_shift_R-best_shift_L)+1:start_ind, :) = [];
            matMotPrim_R(stop_ind+best_shift_R+1:stop_ind, :) = [];
            matMotSecu_L(stop_ind+best_shift_L+1:stop_ind, :) = [];
            matMotSecu_L(start_ind+(best_shift_R-best_shift_L)+1:start_ind, :) = [];
            matMotSecu_R(stop_ind+best_shift_R+1:stop_ind, :) = [];
            
            Fv_FP1(start_ind+best_shift_R+1:start_ind) = [];
            Fv_FP2(start_ind+best_shift_R+1:start_ind) = [];
            Fv_FP3(start_ind+best_shift_R+1:start_ind) = [];
            Fv_FP4(start_ind+best_shift_R+1:start_ind) = [];
            matAllTraj(start_ind+best_shift_R+1:start_ind, :) = [];
            matFP(start_ind+best_shift_R+1:start_ind, :) = [];
            matBut(start_ind+best_shift_R+1:start_ind, :) = [];
            matModOut(start_ind+best_shift_R+1:start_ind, :) = [];
            matGloAng(start_ind+best_shift_R+1:start_ind, :) = [];
            matEMG(start_ind + int32(best_shift_R*f_EMG/f_Mot) + 1:start_ind, :) = [];
            
        elseif (sign(best_shift_L) ~= sign(best_shift_R)) && ((best_shift_L < 0) || (best_shift_R > 0))
            Fv_INL(stop_ind-(best_shift_R-best_shift_L)+1:stop_ind) = [];
            Fv_INR(start_ind-(best_shift_R-best_shift_L)+1:start_ind) = [];
            matMotPrim_L(stop_ind-(best_shift_R-best_shift_L)+1:stop_ind, :) = [];
            matMotPrim_R(start_ind-(best_shift_R-best_shift_L)+1:start_ind, :) = [];
            matMotSecu_L(stop_ind-(best_shift_R-best_shift_L)+1:stop_ind, :) = [];
            matMotSecu_R(start_ind-(best_shift_R-best_shift_L)+1:start_ind, :) = [];
            
            Fv_FP1(start_ind+best_shift_L+1:start_ind) = [];
            Fv_FP1(stop_ind-best_shift_R+1:stop_ind) = [];
            Fv_FP2(start_ind+best_shift_L+1:start_ind) = [];
            Fv_FP2(stop_ind-best_shift_R+1:stop_ind) = [];
            Fv_FP3(start_ind+best_shift_L+1:start_ind) = [];
            Fv_FP3(stop_ind-best_shift_R+1:stop_ind) = [];
            Fv_FP4(start_ind+best_shift_L+1:start_ind) = [];
            Fv_FP4(stop_ind-best_shift_R+1:stop_ind) = [];
            matAllTraj(start_ind+best_shift_L+1:start_ind, :) = [];
            matAllTraj(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matFP(start_ind+best_shift_L+1:start_ind, :) = [];
            matFP(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matBut(start_ind+best_shift_L+1:start_ind, :) = [];
            matBut(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matModOut(start_ind+best_shift_L+1:start_ind, :) = [];
            matModOut(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matGloAng(start_ind+best_shift_L+1:start_ind, :) = [];
            matGloAng(stop_ind-best_shift_R+1:stop_ind, :) = [];
            matEMG(start_ind + int32(best_shift_L*f_EMG/f_Mot) + 1:start_ind, :) = [];
            matEMG(stop_ind - int32(best_shift_R*f_EMG/f_Mot) + 1:stop_ind, :) = [];
    
        elseif (sign(best_shift_L) ~= sign(best_shift_R)) && ((best_shift_R < 0) || (best_shift_L > 0))
            Fv_INL(start_ind-(best_shift_L-best_shift_R)+1:start_ind) = [];
            Fv_INR(stop_ind-(best_shift_L-best_shift_R)+1:stop_ind) = [];
            matMotPrim_L(start_ind-(best_shift_L-best_shift_R)+1:start_ind, :) = [];
            matMotPrim_R(stop_ind-(best_shift_L-best_shift_R)+1:stop_ind, :) = [];
            matMotSecu_L(start_ind-(best_shift_L-best_shift_R)+1:start_ind, :) = [];
            matMotSecu_R(stop_ind-(best_shift_L-best_shift_R)+1:stop_ind, :) = [];
            
            Fv_FP1(start_ind+best_shift_R+1:start_ind) = [];
            Fv_FP1(stop_ind-best_shift_L+1:stop_ind) = [];
            Fv_FP2(start_ind+best_shift_R+1:start_ind) = [];
            Fv_FP2(stop_ind-best_shift_L+1:stop_ind) = [];
            Fv_FP3(start_ind+best_shift_R+1:start_ind) = [];
            Fv_FP3(stop_ind-best_shift_L+1:stop_ind) = [];
            Fv_FP4(start_ind+best_shift_R+1:start_ind) = [];
            Fv_FP4(stop_ind-best_shift_L+1:stop_ind) = [];
            matAllTraj(start_ind+best_shift_R+1:start_ind, :) = [];
            matAllTraj(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matFP(start_ind+best_shift_R+1:start_ind, :) = [];
            matFP(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matBut(start_ind+best_shift_R+1:start_ind, :) = [];
            matBut(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matModOut(start_ind+best_shift_R+1:start_ind, :) = [];
            matModOut(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matGloAng(start_ind+best_shift_R+1:start_ind, :) = [];
            matGloAng(stop_ind-best_shift_L+1:stop_ind, :) = [];
            matEMG(start_ind + int32(best_shift_R*f_EMG/f_Mot) + 1:start_ind, :) = [];
            matEMG(stop_ind - int32(best_shift_L*f_EMG/f_Mot) + 1:stop_ind, :) = [];
        end
    end
    
    % Rejoin left and right matMot
    matMot = [matMotPrim_L, matMotPrim_R, matMotSecu_L, matMotSecu_R];
    clearvars matMotPrim matMotSecu
    clearvars matMotPrim_L matMotPrim_R matMotSecu_L matMotSecu_R
    clearvars t_down T_down t_Mot T_Mot dt_down dt_Mot
    
    n_Mot = size(matMot, 1);
    n_down = n_Mot;
    
    % Visualise
    final_sync_figure = visualize_vertical_forces("Piecewise synchronized signals", ...
        Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, Fv_INL, Fv_INR, NaN, NaN);
    
    
    % Get validity masks
    % ------------------
    
    % Check validity: ocluded markers
    bad_feet_mask = ones(size(matAllTraj, 1), 1);
    for i = 1:size(matAllTraj, 1)
        if sum(isnan([matAllTraj(i, 22:30), matAllTraj(i, 40:48)])) > 0
            bad_feet_mask(i) = 0;
        end
    end
    
    % Find ankle --> global and global --> ankle transformations
    % Returns zero matrix instead of tranformation for occluded markers!!!
    transfos = ankle_to_global(matAllTraj, bad_feet_mask, AW, mark_diam, mark_base);
    
    % Calculate which foot is on which FP at any given time, for final synchronisation
    foot_on_fp_cols = foot_on_fp(-[matFP(:, 3), matFP(:, 10), matFP(:, 17), matFP(:, 24)], ... % Fz data
        [matFP(:, 4:6), matFP(:, 11:13), matFP(:, 18:20), matFP(:, 25:27)], ... % CoP data
        matAllTraj(:, 25:27), matAllTraj(:, 28:30), matAllTraj(:, 43:45), matAllTraj(:, 46:48), ... % LHEE, LTOE, RHEE, RTOE
        ThLOW);
    
    % Define positions of force plates
    pos_FP1 = matFP(1, 4:6);
    pos_FP2 = matFP(1, 11:13);
    pos_FP3 = matFP(1, 18:20);
    pos_FP4 = matFP(1, 25:27);
    
    % % Optionally set them by hand
    % if mat_FP(1, 3) ~= 0
    %     pos_FP1 = [250, 300, 0];
    % end
    % if mat_FP(1, 10) ~= 0
    %     pos_FP2 = [750, 600, 0];
    % end
    % if mat_FP(1, 17) ~= 0
    %     pos_FP3 = [750, 1200, 0];
    % end
    % if mat_FP(1, 24) ~= 0
    %     pos_FP4 = [250, 1200, 0];
    % end

    % Check validity: find bad steps on the force plates
    bad_step_mask = find_bad_steps(matMot, foot_on_fp_cols, insole_length_x, insole_length_y, transfos, Tifl, Tifr, pos_FP1, pos_FP2, pos_FP3, pos_FP4, matAllTraj, BW);

    % Check validity: force plate does/does not give impossible values
    bad_fp_mask = find_bad_fp(matFP, foot_on_fp_cols, pos_FP1, pos_FP2, pos_FP3, pos_FP4);
    
    
    % Transform GRFs
    % --------------
    
    % This means to the original insole frame of reference (as in the tech data
    % sheet), with the assumption that the z-axis is parallel to the global z-axis
    GRF_LR = transform_grfs(matFP, transfos, Tifl, Tifr, foot_on_fp_cols);
    
    % Visualise
    grf_figure = figure;
    p1 = plot(GRF_LR(:, 3));
    t1 = "Left Fz";
    hold on
    p2 = plot(GRF_LR(:, 7));
    t2 = "Right Fz";
    legend([p1, p2], [t1, t2]);
    title("Vertical forces in AP frame of reference");
    hold off

    % Transform COPs
    % --------------

    COPS = transform_cops(matFP, transfos, Tfil, Tfir, foot_on_fp_cols, insole_length_x, insole_length_y);
    
    % plot_slice = [1802:1880, 8555:8636];
    % 
    % figure
    % scatter(COPS(plot_slice, 1), COPS(plot_slice, 2), "green")
    % hold on
    % scatter(matMot(plot_slice, 2), matMot(plot_slice, 3), "red")
    % axis equal
    % hold off
    % 
    % figure
    % scatter(COPS(plot_slice, 3), COPS(plot_slice, 4), "green")
    % hold on
    % scatter(matMot(plot_slice, 5), matMot(plot_slice, 6), "red")
    % axis equal
    % hold off

    
    % Normalise data for GRF references
    % ---------------------------------
    
    % Normalise GRF_LR with repect to body weight and shoe size
    GRF_LR_norm = GRF_LR/BW;
    GRF_LR_norm(:, 4) = -GRF_LR_norm(:, 4)/(insole_length_x); % Original torque expressed in N*mm, sign reversed because expressed in left hand base
    GRF_LR_norm(:, 8) = GRF_LR_norm(:, 8)/(insole_length_x);
    
    % Normalise Moticon Fv estimations
    matMot_norm = matMot;
    matMot_norm(:, 1) = matMot(:, 1)/BW;
    matMot_norm(:, 4) = matMot(:, 4)/BW;
    
    % Normalise Moticon pressure readings
    matMot_norm(:, 7:22) = matMot_norm(:, 7:22)/(BW/(0.1*insole_length_x*0.1*insole_length_y)); % Divide by pressure of full body weight statically on one foot (in N/cm2)
    matMot_norm(:, 29:44) = matMot_norm(:, 29:44)/(BW/(0.1*insole_length_x*0.1*insole_length_y));
    
    % Express Moticon gyroscope readings in rad/s instead of deg/s as then every signal is of same order of magnitude
    matMot_norm(:, 26:28) = matMot_norm(:, 26:28)/(360/(2*pi));
    matMot_norm(:, 48:50) = matMot_norm(:, 48:50)/(360/(2*pi));
       
    % Wait for user input
    % -------------------
    
    txt_input = input("Do you want to continue processing this recording? (y/n)", "s");
    
    if txt_input == "y"
        % Close some windows
        figs2keep = [final_sync_figure, buttons_figure];
        all_figs = findobj(0, 'type', 'figure');
        delete(setdiff(all_figs, figs2keep));
        
        % Find presses in rough and fine segmentation button signal
        ON_OFF_SEGM_ROUGH = block_wave_edges((matBut(:, 2) > 1));
        ON_OFF_SEGM_FINE = block_wave_edges((matBut(:, 3) > 1));
        ON_OFF_SEGM_ROUGH = update_segm_rough(ON_OFF_SEGM_ROUGH, ON_OFF_SEGM_FINE);
        
        % Find presses in error button signal
        ON_OFF_ERROR = block_wave_edges((matBut(:, 1) > 1));
        
        % See if number of executed trials corresponds to number of trial labels
        n_labels = size(trials{i_rec, :}, 2);
        n_trials = size(ON_OFF_SEGM_ROUGH, 1)-1; % -1 because press before and after each task
        n_errors = size(ON_OFF_ERROR, 1);
        
        if n_labels ~= n_trials
            error("Number of executed trials does not match number of trial labels.")
        end
        
        % Set up iteration
        i_trial = 1;
        
        % Loop over all trials
        while (i_trial <= n_trials)
            % Isolate a single trial
            start_ind = ON_OFF_SEGM_ROUGH(i_trial, 1);
            stop_ind = ON_OFF_SEGM_ROUGH(i_trial+1, 2);
        
            % Label of current trial
            trial = trials{i_rec, :};
            trial = trial{:, i_trial};
            
            % Vertical forces from force plates and insoles
            Fv_FP1 = -matFP(start_ind:stop_ind, 3);
            Fv_FP2 = -matFP(start_ind:stop_ind, 10);
            Fv_FP3 = -matFP(start_ind:stop_ind, 17);
            Fv_FP4 = -matFP(start_ind:stop_ind, 24);
            
            Fv_INL = matMot(start_ind:stop_ind, 1);
            Fv_INR = matMot(start_ind:stop_ind, 4);
            
            % Visualise before further processing
            single_task_figure = visualize_vertical_forces(trial, ...
                Fv_FP1, Fv_FP2, Fv_FP3, Fv_FP4, Fv_INL, Fv_INR, NaN, NaN);
            
            % Wait for user input to continue
            txt_input = input("Do you want to continue processing this task? (y/n)", "s");
            %txt_input = "y";
    
            % Cleared to continue
            if txt_input == "y"
                % Save synchronised data per trial execution, per data 'type'
                filename_base_Trn = procDataFolderTrn + "\" + subject_tag + " - Rec" + string(i_rec) + " - " + trial + " - ";
                filename_base_Ref = procDataFolderRef + "\" + subject_tag + " - Rec" + string(i_rec) + " - " + trial + " - ";
    
                % Isolate single trial executions
                start_stop_reps = [];
                for i = 1:2:size(ON_OFF_SEGM_FINE, 1)
                    start_ind_rep = ON_OFF_SEGM_FINE(i, 2);
                    stop_ind_rep = ON_OFF_SEGM_FINE(i+1, 1);
                    if (start_ind_rep >= start_ind) && (stop_ind_rep <= stop_ind)
                        start_stop_reps = [start_stop_reps; [start_ind_rep, stop_ind_rep]];
                    end
                end

                % Filter out trial executions where error button was pressed
                start_stop_reps = update_start_stop_reps(start_stop_reps, ON_OFF_ERROR);
                
                % Save data for each execution
                for i_rep = 1:size(start_stop_reps, 1)
                    start_ind_rep = start_stop_reps(i_rep, 1);
                    stop_ind_rep = start_stop_reps(i_rep, 2);
                    n_samp_rep = stop_ind_rep-start_ind_rep+1;
                    
                    
                    % Training data
                    % -------------

                    % Trajectory data
                    tabAllTraj = array2table(matAllTraj(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesAllTraj);
                    writetable(tabAllTraj, filename_base_Trn + "Rep" + string(i_rep) + " - " + "AllTraj.csv");
    
                    % Button data
                    tabBut = array2table(matBut(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesBut);
                    writetable(tabBut, filename_base_Trn + "Rep" + string(i_rep) + " - " + "But.csv");
    
                    % EMG data
                    tabEMG = array2table(matEMG(int32((start_ind_rep-1)*f_EMG/f_Mot)+1:int32(stop_ind_rep*f_EMG/f_Mot), :), 'VariableNames', colnamesEMG);
                    writetable(tabEMG, filename_base_Trn + "Rep" + string(i_rep) + " - " + "EMG.csv");
    
                    % Force plate data
                    tabFP = array2table(matFP(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesFP);
                    writetable(tabFP, filename_base_Trn + "Rep" + string(i_rep) + " - " + "FP.csv");
    
                    % Model output data
                    tabModOut = array2table(matModOut(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesModOut);
                    writetable(tabModOut, filename_base_Trn + "Rep" + string(i_rep) + " - " + "ModOut.csv");
                    
                    % Global segment angle data
                    tabGloAng = array2table(matGloAng(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesGloAng);
                    writetable(tabGloAng, filename_base_Trn + "Rep" + string(i_rep) + " - " + "GloAng.csv");

                    % Moticon data
                    tabMot = array2table(matMot(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesMot);
                    writetable(tabMot, filename_base_Trn + "Rep" + string(i_rep) + " - " + "Mot.csv");
                    tabMot_norm = array2table(matMot_norm(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesMot);
                    writetable(tabMot_norm, filename_base_Trn + "Rep" + string(i_rep) + " - " + "MotNorm.csv");

                    % Masks
                    tabMasks = array2table([bad_feet_mask, bad_fp_mask, bad_step_mask], 'VariableNames', colnamesMasks);
                    writetable(tabMasks, filename_base_Trn + "Rep" + string(i_rep) + " - " + "Masks.csv");

                    % Foot on force plate data
                    tabFOF = array2table(foot_on_fp_cols(start_ind_rep:stop_ind_rep, :), 'VariableNames', colnamesFOF);
                    writetable(tabFOF, filename_base_Trn + "Rep" + string(i_rep) + " - " + "FootOnFP.csv");


                    % Reference data
                    % --------------
                    
                    % Find parts usable as reference data
                    segm_vec = cell(0, 3);
                    j = start_ind_rep;
                    task_number = 0;
                    while (j <= stop_ind_rep)
                        if check_active([GRF_LR(:, 3), GRF_LR(:, 7)], j, stop_ind_rep, n_buf, 1, ThLOW)
                            task_number = task_number+1;
                            single_start = j;
                            k = 1;
                            while ((check_active([GRF_LR(:, 3), GRF_LR(:, 7)], j+k, stop_ind_rep, n_buf, 1, ThLOW) || check_active([GRF_LR(:, 3), GRF_LR(:, 7)], stop_ind_rep-(j+k), stop_ind_rep, n_buf, 0, ThLOW))) && ((j+k) < stop_ind_rep)
                                k = k+1;
                            end
                            j = j+k;
                            single_stop = j;
                            
                            % single_grf_figure = figure;
                            % plot(GRF_LR(single_start:single_stop, 3))
                            % hold on
                            % plot(GRF_LR(single_start:single_stop, 7))
                            % hold off
                
                            filename = filename_base_Ref + "Rep" + string(i_rep) + "_" + string(task_number);
                            segm_vec = [segm_vec; {single_start, single_stop, filename}];
                        end
                        j = j+1;
                    end
                    
                    % Save reference data for each of those parts
                    n_segm = size(segm_vec, 1);
                    if n_segm ~= 0
                        for j=1:n_segm
                            single_start = segm_vec{j, 1};
                            single_stop = segm_vec{j, 2};
            
                            crit1 = sum(not(bad_feet_mask(single_start:single_stop)));
                            crit2 = sum(not(bad_fp_mask(single_start:single_stop)));
                            crit3 = sum(not(bad_step_mask(single_start:single_stop)));
            
                            filename = segm_vec{j, 3};
                            filename = append_flags(filename, crit1, crit2, crit3);
                            filename = filename + ".csv";
            
                            % actually save data
                            save_matrix = [GRF_LR_norm(single_start:single_stop, :), COPS(single_start:single_stop, :), matMot_norm(single_start:single_stop, :), foot_on_fp_cols(single_start:single_stop, :), bad_feet_mask(single_start:single_stop), bad_fp_mask(single_start:single_stop), bad_step_mask(single_start:single_stop)];
                            tabGRF = array2table(save_matrix, 'VariableNames', colnamesGRF);
                            writetable(tabGRF, filename);

                            % also save model output data for this section
                            filename = segm_vec{j, 3};
                            filename = filename + " - ModelOutputs";
                            filename = append_flags(filename, crit1, crit2, crit3);
                            filename = filename + ".csv";
                            tabModOut = array2table(matModOut(single_start:single_stop, :), 'VariableNames', colnamesModOut);
                            writetable(tabModOut, filename);

                            % also save global segment angle data for this section
                            filename = segm_vec{j, 3};
                            filename = filename + " - GlobalAngles";
                            filename = append_flags(filename, crit1, crit2, crit3);
                            filename = filename + ".csv";
                            tabGloAng = array2table(matGloAng(single_start:single_stop, :), 'VariableNames', colnamesGloAng);
                            writetable(tabGloAng, filename);

                            % also save the foot marker data for this section
                            filename = segm_vec{j, 3};
                            filename = filename + " - FootMarkers";
                            filename = append_flags(filename, crit1, crit2, crit3);
                            filename = filename + ".csv";
                            tabFootMarkers = array2table([matAllTraj(single_start:single_stop, 22:30), matAllTraj(single_start:single_stop, 40:48)], 'VariableNames', colnamesFootMarkers);
                            writetable(tabFootMarkers, filename);

                            % also save EMG data for this section
                            filename = segm_vec{j, 3};
                            filename = filename + " - EMG";
                            filename = append_flags(filename, crit1, crit2, crit3);
                            filename = filename + ".csv";
                            tabEMG = array2table(matEMG(int32((single_start-1)*f_EMG/f_Mot)+1:int32(single_stop*f_EMG/f_Mot), :), 'VariableNames', colnamesEMG);
                            writetable(tabEMG, filename);
                        end
                    end
                end
                
            % Willingly disregard trial
            elseif txt_input == "n"
                close(single_task_figure)
                fprintf("Trial %d (%s) disregarded.\n", i_trial, trial{1, 1})
            
            % Invalid input --> disregard trial
            else
                close(single_task_figure)
                fprintf("Invalid user input, trial %d (%s) disregarded.\n", i_trial, trial{1, 1})
            end
            
            % Close some windows
            figs2keep = [final_sync_figure, buttons_figure];
            all_figs = findobj(0, 'type', 'figure');
            delete(setdiff(all_figs, figs2keep));
            
            % Increment
            i_trial = i_trial+1;
        end
    
        % Close all windows after processing
        close all
    
    % Willingly disregard recording
    elseif txt_input == "n"
        close all
        fprintf("Recording %d disregarded.\n", i_rec)
    
    % Invalid input --> disregard recording
    else
        close all
        fprintf("Invalid user input, recording %d disregarded.\n", i_rec)
    end
end